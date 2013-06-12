module selection

  use detypes

  implicit none

#ifdef USEMPI
  include 'mpif.h'
#endif

  private
  public selector, replace_generation, roundvector

  logical, parameter :: debug_replace_gen=.false.

contains 

  subroutine selector(X, Xnew, U, trialF, trialCr, m, n, lowerbounds, upperbounds, &
                       run_params, fcall, func, quit, accept)

    type(population), intent(in) :: X
    type(population), intent(inout) :: Xnew
    integer, intent(inout) :: fcall, accept
    real(dp), dimension(:), intent(in) :: U
    real(dp), intent(in) :: trialF, trialCr
    integer, intent(in) :: m, n              !current index for population chunk (m) and full population (n)
    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds 
    type(codeparams), intent(in) :: run_params
    logical, intent(inout) :: quit
    real(dp), external :: func

    real(dp) :: trialvalue
    real(dp), dimension(size(U)) :: trialvector
    real(dp), dimension(size(X%vectors_and_derived(1,:))) :: trialderived
    logical :: validvector


    trialderived = 0.

    if (any(U(:) .gt. upperbounds) .or. any(U(:) .lt. lowerbounds)) then 
                                             !trial vector exceeds parameter space bounds: apply boundary constraints
       select case (run_params%DE%bconstrain)
          case (1)                           !'brick wall'
             trialvector(:) = X%vectors(n,:)
             validvector = .false.
          case (2)                           !randomly re-initialize
             call random_number(trialvector(:))
             trialvector(:) = trialvector(:)*(upperbounds - lowerbounds) + lowerbounds
             validvector = .true.
          case (3)                           !reflection
             trialvector = U
             where (U .gt. upperbounds) trialvector = upperbounds - (U - upperbounds)
             where (U .lt. lowerbounds) trialvector = lowerbounds + (lowerbounds - U)
             validvector = .true.
          case default                       !boundary constraints not enforced
             trialvector = U  
             validvector = .true.
          end select
    else                                     !trial vector is within parameter space bounds, so use it
       trialvector = U 
       validvector = .true.
    end if

    trialderived(:run_params%D) = roundvector(trialvector, run_params)

#ifdef USEMPI 
    !call func even if not a valid vector so that all processes can respond to MPI calls inside likelihood function 
    trialvalue = func(trialderived, fcall, quit, validvector)
    if (.not. validvector) trialvalue = huge(1.0_dp)
#else   
    !save time by only evaluating likelihood when necessary
    if (validvector) then 
       trialvalue = func(trialderived, fcall, quit, validvector)
    else
       trialvalue = huge(1.0_dp)
    end if
#endif

    !when the trial vector is at least as good as the current member  
    !of the population, use the trial vector for the next generation
    if (trialvalue .le. X%values(n)) then
       Xnew%vectors(m,:) = trialvector 
       Xnew%vectors_and_derived(m,:) = trialderived
       Xnew%values(m) = trialvalue
       if (run_params%DE%jDE) then            !in jDE, also keep F and Cr
          Xnew%FjDE(m) = trialF
          Xnew%CrjDE(m) = trialCr
       end if
       accept = accept + 1
    else
       Xnew%vectors(m,:) = X%vectors(n,:) 
       Xnew%vectors_and_derived(m,:) = X%vectors_and_derived(n,:)
       Xnew%values(m) = X%values(n)
       if (run_params%DE%jDE) then
          Xnew%FjDE(m) = X%FjDE(n)
          Xnew%CrjDE(m) = X%CrjDE(n)
       end if
    end if

  end subroutine selector


  !rounds vectors to nearest discrete values for all dimensions listed in run_params%discrete
  !all other dimensions are kept the same
  function roundvector(trialvector, run_params)
    real(dp), dimension(:), intent(in) :: trialvector
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D) :: roundvector

    roundvector = trialvector
    roundvector(run_params%discrete) = anint(roundvector(run_params%discrete))

  end function roundvector


  !replaces old generation (X) with the new generation (Xnew) calculated during population loop
  subroutine replace_generation(X, Xnew, run_params, accept, init)
    type(population), intent(inout) :: X              !old population, will be replaced
    type(population), intent(inout) :: Xnew           !recently calculated population chunk
    type(codeparams), intent(in) :: run_params
    integer, intent(inout) :: accept
    logical, intent(in) :: init
    real(dp), dimension(run_params%DE%NP, run_params%D) :: allvecs   !new vector population. For checking for duplicates
    real(dp), dimension(run_params%D, run_params%DE%NP) :: trallvecs !transposed allvecs, to make MPI_Allgather happy
    real(dp), dimension(run_params%DE%NP) :: allvals                 !new values corresponding to allvecs. For checking for duplicates
    real(dp), dimension(run_params%D+run_params%D_derived, run_params%DE%NP) :: trderived !transposed derived
    integer :: k, kmatch                                             !indices for vector compared, possible matching vector
    integer :: ierror, mpi_dp  
    
    !with MPI enabled, Xnew will only contain some elements of the new population. Create allvecs, allvals for duplicate-hunting
#ifdef USEMPI
    !create mpi double precsision real which will be compatible with dp kind specified in detypes.f90
    call MPI_Type_create_f90_real(precision(1.0_dp), range(1.0_dp), mpi_dp, ierror) 

    call MPI_Allgather(transpose(Xnew%vectors), run_params%mpipopchunk*run_params%D, mpi_dp, trallvecs, &
                       run_params%mpipopchunk*run_params%D, mpi_dp, MPI_COMM_WORLD, ierror)
    allvecs = transpose(trallvecs)

    call MPI_Allgather(Xnew%values, run_params%mpipopchunk, mpi_dp, allvals, & 
                       run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)  
#else
    allvecs = Xnew%vectors
    allvals = Xnew%values
#endif

    !weed out any duplicate vectors to maintain population diversity. One duplicate will be kept and the other will revert to 
    !its value in the previous generation (NB for discrete dimensions, we are comparing the underlying non-discrete vectors)
    if (run_params%DE%removeDuplicates .and. .not. init) then
       checkpop: do k=1, run_params%DE%NP-1                                        !look for matches in 1st dim of higher-indexed Xnew%vectors  
          if ( any(allvecs(k,1) .eq. allvecs(k+1:run_params%DE%NP,1)) ) then       !there is at least one possible match

             findmatch: do kmatch=k+1, run_params%DE%NP                            !loop over subpopulation to find the matching vector(s)
                if ( all(allvecs(k,:) .eq. allvecs(kmatch,:)) ) then               !we've found a duplicate vector
                   if (verbose) write (*,*) '  Duplicate vectors:', k, kmatch

                   !Now, compare their counterparts in the previous generation to decide which vector will be kept, which will be reverted
                   picksurvivor: if (all(allvecs(k,:) .eq. X%vectors(k,:)) ) then  !vector at k was inherited, so keep it & revert kmatch
                      if (verbose) write (*,*) '  Reverting vector ', kmatch
                      allvecs(kmatch,:) = X%vectors(kmatch,:)
                      allvals(kmatch) = X%values(kmatch)                     
                      call replace_vector(Xnew, X, run_params, kmatch, accept)

                   else if (all(allvecs(kmatch,:) .eq. X%vectors(kmatch,:))) then  !vector at kmatch was inherited. Keep it
                      if (verbose) write (*,*) '  Reverting vector ', k
                      allvecs(k,:) = X%vectors(k,:)
                      allvals(k) = X%values(k)
                      call replace_vector(Xnew, X, run_params, k, accept) 

                   else if (X%values(k) .lt. X%values(kmatch)) then                !kmatch improved more, so keep it
                      if (verbose) write (*,*) '  Reverting vector ', k
                      allvecs(k,:) = X%vectors(k,:)
                      allvals(k) = X%values(k)
                      call replace_vector(Xnew, X, run_params, k, accept) 

                   else                                                            !k improved more (or the same), so keep it
                      if (verbose) write (*,*) '  Reverting vector ', kmatch
                      allvecs(kmatch,:) = X%vectors(kmatch,:)
                      allvals(kmatch) = X%values(kmatch)
                      call replace_vector(Xnew, X, run_params, kmatch, accept) 

                   end if picksurvivor
                end if

             end do findmatch

          end if
       end do checkpop
    end if


    !replace old population members with those calculated in Xnew
#ifdef USEMPI
    if (debug_replace_gen) then !this just compares the replaced Xnew%vectors & Xnew%values with allvecs and allvals
       call MPI_Allgather(transpose(Xnew%vectors), run_params%mpipopchunk*run_params%D, mpi_dp, trallvecs, &
                          run_params%mpipopchunk*run_params%D, mpi_dp, MPI_COMM_WORLD, ierror)
       X%vectors = transpose(trallvecs)
       if (any(X%vectors .ne. allvecs)) write (*,*) 'ERROR: vectors not transferred properly'
       
       call MPI_Allgather(Xnew%values, run_params%mpipopchunk, mpi_dp, X%values, & 
                          run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
       if (any(X%values .ne. allvals)) write (*,*) 'ERROR: values not transferred properly'

    else                        !vectors and values have already been gathered
       X%vectors = allvecs
       X%values = allvals
    end if
    	
    call MPI_Allgather(transpose(Xnew%vectors_and_derived), run_params%mpipopchunk*&
                       (run_params%D+run_params%D_derived), mpi_dp, trderived, &
                       run_params%mpipopchunk*(run_params%D+run_params%D_derived), &
                       mpi_dp, MPI_COMM_WORLD, ierror)
    X%vectors_and_derived = transpose(trderived)

    if (run_params%DE%jDE) then
       call MPI_Allgather(Xnew%FjDE, run_params%mpipopchunk, mpi_dp, X%FjDE, & 
                          run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
       call MPI_Allgather(Xnew%CrjDE, run_params%mpipopchunk, mpi_dp, X%CrjDE, & 
                          run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
    end if
#else
    !Xnew and X are the same size, so just equate population members
    X%vectors = allvecs
    X%values = allvals
    X%vectors_and_derived = Xnew%vectors_and_derived
    if (run_params%DE%jDE) then
       X%FjDE = Xnew%FjDE
       X%CrjDE = Xnew%CrjDE
    end if
#endif

  end subroutine replace_generation


!replace a vector in Xnew by its counterpart in the previous generation (X)
  subroutine replace_vector(Xnew, X, run_params, n, accept)
    type(population), intent(inout) :: Xnew
    type(population), intent(in) :: X
    type(codeparams), intent(in) :: run_params
    integer, intent(in) :: n                                     !index of vector X to replace
    integer, intent(inout) :: accept
    integer :: m                                                 !index of vector in Xnew (equal to n if no MPI)
    
    m = n - run_params%mpipopchunk*run_params%mpirank

    if ( (m .gt. 0) .and. (m .le. run_params%mpipopchunk) ) then !vector belongs to population chunk in this process
       
       if (debug_replace_gen) then  !checking that this transfer works correctly
          Xnew%vectors(m,:) = X%vectors(n,:)
          Xnew%values(m) = X%values(n)      
       end if

       Xnew%vectors_and_derived(m,:) = X%vectors_and_derived(n,:)

       if (run_params%DE%jDE) then
          Xnew%FjDE(m) = X%FjDE(n)
          Xnew%CrjDE(m) = X%CrjDE(n)
       end if

       if (verbose) write (*,*) n, roundvector(Xnew%vectors(m, :), run_params), '->', Xnew%values(m)

       accept = accept - 1                                       !vector has been 'de-accepted'
    end if
    
  end subroutine replace_vector


end module selection
