module selection

use detypes
use deutils
use mutation, only: init_FjDE, init_lambdajDE
use crossover, only: init_CrjDE

implicit none

#ifdef MPI
  include 'mpif.h'
#endif

private
public selector, replace_generation, roundvector

contains

  subroutine selector(X, Xnew, U, trialF, triallambda, trialCr, m, n, run_params, func, fcall, quit, &
       accept, acceptable_trial_vector)

    type(population), intent(in) :: X
    type(population), intent(inout) :: Xnew
    integer, intent(inout) :: fcall, accept
    real(dp), dimension(:), intent(in) :: U
    real(dp), intent(in) :: trialF, triallambda, trialCr
    integer, intent(in) :: m, n              !current index for population chunk (m) and full population (n)
    type(codeparams), intent(inout) :: run_params
    logical, intent(out) :: acceptable_trial_vector
    logical, intent(inout) :: quit
    procedure(MinusLogLikeFunc) :: func

    real(dp) :: trialvalue
    real(dp), dimension(size(U)) :: trialvector
    real(dp), dimension(size(X%vectors_and_derived(1,:))) :: trialderived
    logical :: validvector, acceptable_fitness


    trialderived = 0.

    if (any(U(:) .gt. run_params%upperbounds) .or. any(U(:) .lt. run_params%lowerbounds)) then
                                             !trial vector exceeds parameter space bounds: apply boundary constraints
       select case (run_params%DE%bconstrain)
          case (1)                           !'brick wall'
             trialvector(:) = X%vectors(n,:)
             validvector = .false.
          case (2)                           !randomly re-initialize
             call random_number(trialvector(:))
             trialvector(:) = trialvector(:)*(run_params%upperbounds - run_params%lowerbounds) + run_params%lowerbounds
             if (run_params%partitionDiscrete) then   !when partitioning, keep old values for discrete variables
                trialvector(run_params%discrete) = U(run_params%discrete)
             end if
             validvector = .true.
          case (3)                           !reflection
             trialvector = U
             where (U .gt. run_params%upperbounds) trialvector = run_params%upperbounds - (U - run_params%upperbounds)
             where (U .lt. run_params%lowerbounds) trialvector = run_params%lowerbounds + (run_params%lowerbounds - U)
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

#ifdef MPI
    !call func even if not a valid vector so that all processes can respond to MPI calls inside likelihood function
    trialvalue = func(trialderived, fcall, quit, validvector, run_params%context)
    if (.not. validvector) trialvalue = huge(1.0_dp)
#else
    !save time by only evaluating likelihood when necessary
    if (validvector) then
       trialvalue = func(trialderived, fcall, quit, validvector, run_params%context)
    else
       trialvalue = huge(1.0_dp)
    end if
#endif

    acceptable_fitness = (trialvalue .le. run_params%max_acceptable_value) !true if point is acceptably fit
    acceptable_trial_vector = (acceptable_fitness .or. .not. run_params%discard_unfit_points) !true if trial vector is acceptable

    if (acceptable_trial_vector) then
       !when the trial vector is at least as good as the current member
       !of the population, use the trial vector for the next generation
       if (trialvalue .le. X%values(n)) then
          Xnew%vectors(m,:) = trialvector
          Xnew%vectors_and_derived(m,:) = trialderived
          Xnew%values(m) = trialvalue
          if (run_params%DE%jDE) then  !in jDE, also keep F and Cr (and lambda)
             Xnew%FjDE(m) = trialF
             Xnew%CrjDE(m) = trialCr
             if (run_params%DE%lambdajDE) then
                Xnew%lambdajDE(m) = triallambda
             end if
          end if
          accept = accept + 1
       else
          Xnew%vectors(m,:) = X%vectors(n,:)
          Xnew%vectors_and_derived(m,:) = X%vectors_and_derived(n,:)
          Xnew%values(m) = X%values(n)
          if (run_params%DE%jDE) then
             Xnew%FjDE(m) = X%FjDE(n)
             Xnew%CrjDE(m) = X%CrjDE(n)
             if (run_params%DE%lambdajDE) then
                Xnew%lambdajDE(m) = X%lambdajDE(n)
             end if
          end if
       end if
    end if

  end subroutine selector



  !replaces old generation (X) with the new generation (Xnew) calculated during population loop
  subroutine replace_generation(X, Xnew, run_params, func, fcall, quit, accept, init)
    type(population), intent(inout) :: X              !old population, will be replaced
    type(population), intent(inout) :: Xnew           !recently calculated population chunk
    type(codeparams), intent(inout) :: run_params
    procedure(MinusLogLikeFunc) :: func
    integer, intent(inout) :: fcall, accept
    logical, intent(inout) :: quit
    logical, intent(in) :: init
    real(dp), dimension(run_params%DE%NP, run_params%D) :: allvecs   !new vector population. For checking for duplicates
    real(dp), dimension(run_params%D, run_params%DE%NP) :: trallvecs !transposed allvecs, to make MPI_Allgather happy
    real(dp), dimension(run_params%D+run_params%D_derived, run_params%DE%NP) :: trderived !transposed derived
    integer :: ierror

    !with MPI enabled, Xnew will only contain some elements of the new population. Create allvecs, allvals for duplicate-hunting
#ifdef MPI
    !TODO
    !Create mpi double precision real, which will be compatible with dp kind specified in detypes.f90,
    !with some sort of configure switch in future.  This should be
    !call MPI_Type_create_f90_real(precision(1.0_dp), range(1.0_dp), mpi_dp, ierror)
    !call MPI_Allgather(transpose(Xnew%vectors), run_params%mpipopchunk*run_params%D, mpi_dp, trallvecs, &
    !                   run_params%mpipopchunk*run_params%D, mpi_dp, MPI_COMM_WORLD, ierror)
    !This cannot be done currently, due to an MPICH bug: http://trac.mpich.org/projects/mpich/ticket/1769

    call MPI_Allgather(transpose(Xnew%vectors), run_params%mpipopchunk*run_params%D, mpi_double_precision, trallvecs, &
                       run_params%mpipopchunk*run_params%D, mpi_double_precision, MPI_COMM_WORLD, ierror)
    allvecs = transpose(trallvecs)

    !weed out duplicate vectors if desired
    if (run_params%DE%removeDuplicates) then
       call remove_duplicate_vectors(X, Xnew, run_params, allvecs, func, fcall, quit, accept, init)
    end if


    !replace old population members with those calculated in Xnew
    X%vectors = allvecs

    !call MPI_Allgather(Xnew%values, run_params%mpipopchunk, mpi_dp, X%values, &
     !                  run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
    call MPI_Allgather(Xnew%values, run_params%mpipopchunk, mpi_double_precision, X%values, &
                       run_params%mpipopchunk, mpi_double_precision, MPI_COMM_WORLD, ierror)

    !call MPI_Allgather(transpose(Xnew%vectors_and_derived), run_params%mpipopchunk*&
    !                   (run_params%D+run_params%D_derived), mpi_dp, trderived, &
    !                   run_params%mpipopchunk*(run_params%D+run_params%D_derived), &
    !                   mpi_dp, MPI_COMM_WORLD, ierror)
    call MPI_Allgather(transpose(Xnew%vectors_and_derived), run_params%mpipopchunk*&
                       (run_params%D+run_params%D_derived), mpi_double_precision, trderived, &
                       run_params%mpipopchunk*(run_params%D+run_params%D_derived), &
                       mpi_double_precision, MPI_COMM_WORLD, ierror)
    X%vectors_and_derived = transpose(trderived)

    if (run_params%DE%jDE) then
       !call MPI_Allgather(Xnew%FjDE, run_params%mpipopchunk, mpi_dp, X%FjDE, &
       !                   run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
       call MPI_Allgather(Xnew%FjDE, run_params%mpipopchunk, mpi_double_precision, X%FjDE, &
                        run_params%mpipopchunk, mpi_double_precision, MPI_COMM_WORLD, ierror)
       !call MPI_Allgather(Xnew%CrjDE, run_params%mpipopchunk, mpi_dp, X%CrjDE, &
       !                   run_params%mpipopchunk, mpi_dp, MPI_COMM_WORLD, ierror)
       call MPI_Allgather(Xnew%CrjDE, run_params%mpipopchunk, mpi_double_precision, X%CrjDE, &
                          run_params%mpipopchunk, mpi_double_precision, MPI_COMM_WORLD, ierror)
       if (run_params%DE%lambdajDE) then
          call MPI_Allgather(Xnew%lambdajDE, run_params%mpipopchunk, mpi_double_precision, X%lambdajDE, &
                          run_params%mpipopchunk, mpi_double_precision, MPI_COMM_WORLD, ierror)
       end if
    end if
#else
    allvecs = Xnew%vectors

    !weed out duplicate vectors if desired
    if (run_params%DE%removeDuplicates) then
       call remove_duplicate_vectors(X, Xnew, run_params, allvecs, func, fcall, quit, accept, init)
    end if

    !Xnew and X are the same size, so just equate population members
    X%vectors = allvecs
    X%values = Xnew%values
    X%vectors_and_derived = Xnew%vectors_and_derived
    if (run_params%DE%jDE) then
       X%FjDE = Xnew%FjDE
       X%CrjDE = Xnew%CrjDE
       if (run_params%DE%lambdajDE) then
          X%lambdajDE = Xnew%lambdajDE
       end if
    end if
#endif

  end subroutine replace_generation


 !weed out any duplicate vectors to maintain population diversity. One duplicate will be kept and the other will revert to
 !its value in the previous generation (NB for discrete dimensions, we are comparing the underlying non-discrete vectors)
  subroutine remove_duplicate_vectors(X, Xnew, run_params, allvecs, func, fcall, quit, accept, init)

    type(population), intent(inout) :: X              !old population, will be replaced
    type(population), intent(inout) :: Xnew           !recently calculated population chunk
    type(codeparams), intent(inout) :: run_params
    procedure(MinusLogLikeFunc) :: func
    integer, intent(inout) :: fcall, accept
    logical, intent(inout) :: quit
    logical, intent(in) :: init
    real(dp), intent(inout), dimension(run_params%DE%NP, run_params%D) :: allvecs

    integer :: k, kmatch                              !indices for vector compared, possible matching vector

    checkpop: do k=1, run_params%DE%NP-1                                        !look for matches in 1st dim of higher-indexed Xnew%vectors
       if ( any(allvecs(k,1) .eq. allvecs(k+1:run_params%DE%NP,1)) ) then       !there is at least one possible match
          !if  (run_params%mpirank .eq. 0) write (*,*) 'WARNING: Possible match for', k

          findmatch: do kmatch=k+1, run_params%DE%NP                            !loop over subpopulation to find the matching vector(s)
             if ( all(allvecs(k,:) .eq. allvecs(kmatch,:)) ) then               !found a duplicate vector (do all()-->any() to prevent dupes in single dimensions)
                if (run_params%verbose .ge. 3) write (*,*) '  Duplicate vectors:', k, kmatch

                !Now, compare their counterparts in the previous generation to decide which vector will be kept, which will be reverted
                picksurvivor: if (init) then
                   if (run_params%verbose .ge. 1) write (*,*) 'WARNING: Duplicate vectors in initial generation'  !This should never happen.
                   if (run_params%verbose .ge. 3) write (*,*) '  Generating new vector ', kmatch
                   !replacing second vector with a random new one
                   call replace_vector(Xnew, allvecs, X, run_params, func, kmatch, fcall, quit, accept, revert=.false.)

                else if (all(allvecs(k,:) .eq. X%vectors(k,:)) .and. &                         !both vectors were inherited, so keep k & randomly re-pick kmatch
                           all(allvecs(kmatch,:) .eq. X%vectors(kmatch,:))) then
                    if (run_params%verbose .ge. 1) write (*,*) 'WARNING: Duplicate vectors inherited from previous generation'
                        !Shouldn't happen, but...
                        !can occur if one vector inherited, other reverted in previous generation to a vector which matches the first
                        !This can be a sign that single-dimension duplicates are polluting the population, or that you're just unlucky
                    if (run_params%verbose .ge. 3) write (*,*) '  Generating new vector ', kmatch              !replacing second vector with a random new one
                   call replace_vector(Xnew, allvecs, X, run_params, func, kmatch, fcall, quit, accept, revert=.false.)

                else if (all(allvecs(k,:) .eq. X%vectors(k,:)) ) then                          !vector at k was inherited, so keep it & revert kmatch
                   if (run_params%verbose .ge. 3) then
                      write (*,*) '    Vector ', k, ' inherited, reverting vector ', kmatch
                   end if
                   call replace_vector(Xnew, allvecs, X, run_params, func, kmatch, fcall, quit, accept, revert=.true.)

                else if (all(allvecs(kmatch,:) .eq. X%vectors(kmatch,:))) then                 !vector at kmatch was inherited. Keep it
                   if (run_params%verbose .ge. 3) then
                      write (*,*) '    Vector ', kmatch, ' inherited, reverting vector ', k
                   end if
                   call replace_vector(Xnew, allvecs, X, run_params, func, k, fcall, quit, accept, revert=.true.)

                else if (X%values(k) .lt. X%values(kmatch)) then                               !kmatch improved more (or the same), so keep it
                   if (run_params%verbose .ge. 3) then
                      write (*,*) '    Vector ', kmatch, ' improved more, reverting vector', k
                   end if
                   call replace_vector(Xnew, allvecs, X, run_params, func, k, fcall, quit, accept, revert=.true.)

                else                                                                           !k improved more, so keep it
                   if (run_params%verbose .ge. 3) then
                      write (*,*) '    Vector ', k, ' improved more, reverting vector', kmatch
                   endif
                   call replace_vector(Xnew, allvecs, X, run_params, func, kmatch, fcall, quit, accept, revert=.true.)

                end if picksurvivor
             end if

          end do findmatch

       end if
    end do checkpop

  end subroutine remove_duplicate_vectors


!replace a duplicate vector in Xnew and allvecs by:
!its counterpart in the previous generation (X) (if revert=.true.) or a new randomly generated vector
  subroutine replace_vector(Xnew, allvecs, X, run_params, func, n, fcall, quit, accept, revert)
    type(population), intent(inout) :: Xnew
    type(codeparams), intent(inout) :: run_params
    real(dp), intent(inout), dimension(run_params%DE%NP, run_params%D) :: allvecs
    type(population), intent(in) :: X
    integer, intent(in) :: n                                     !index of vector X to replace
    procedure(MinusLogLikeFunc) :: func
    integer, intent(inout) :: fcall, accept
    logical, intent(in) :: revert
    logical, intent(inout) :: quit
    integer :: m, root
    real(dp), dimension(run_params%D) :: newvector               !alias for Xnew(m,:) for sharing between processes
    real(dp), dimension(1) :: Fnew, lambdanew, Crnew
    real(dp) :: rand
    integer :: ierror, i !mpi_dp !TODO does not work atm due to an MPICH bug: http://trac.mpich.org/projects/mpich/ticket/1769


    m = n - run_params%mpipopchunk*run_params%mpirank            !index of vector in Xnew (equal to n if no MPI)

    root = (n-1)/run_params%mpipopchunk                          !process which 'owns' the vector being replaced

    !vector belongs to population chunk in this process, so change vector and all associated quantities
    if (run_params%mpirank .eq. root) then

       if (revert) then                                          !replace with previous values
          Xnew%vectors(m,:) = X%vectors(n,:)
          Xnew%values(m) = X%values(n)
          Xnew%vectors_and_derived(m,:) = X%vectors_and_derived(n,:)

          if (run_params%DE%jDE) then
             Xnew%FjDE(m) = X%FjDE(n)
             Xnew%CrjDE(m) = X%CrjDE(n)
             if (run_params%DE%lambdajDE) then
                Xnew%lambdajDE(m) = X%lambdajDE(n)
             end if
          end if

       else                                                      !randomly generate a new vector
          do i = 1, run_params%D
             !Keep the old values for partitioned discrete parameters
             if (run_params%partitionDiscrete .and. any(run_params%discrete .eq. i)) cycle
             call random_number(rand)
             Xnew%vectors(m,i) = rand*(run_params%upperbounds(i) - run_params%lowerbounds(i)) + run_params%lowerbounds(i)
          enddo
          Xnew%vectors_and_derived(m,:run_params%D) = roundvector(Xnew%vectors(m,:), run_params)
          Xnew%values(m) = func(Xnew%vectors_and_derived(m,:), fcall, quit, .true., run_params%context)
          if (run_params%DE%jDE) then
             Fnew = init_FjDE(1)
             Crnew = init_CrjDE(1)
             Xnew%FjDE(m) = Fnew(1)
             Xnew%CrjDE(m) = Crnew(1)
             if (run_params%DE%lambdajDE) then
                lambdanew = init_lambdajDE(1)
                Xnew%lambdajDE(m) = lambdanew(1)
             end if
          end if
          newvector = Xnew%vectors(m,:)
       end if

       if (abs(run_params%verbose) .ge. 3) then
          write (*,*) '    Replacement vector:', n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m)
       end if

    end if

    !All processes: fix values in allvecs corresponding to Xnew%vectors(m,:)
    if (revert) then
       allvecs(n,:) = X%vectors(n,:)                          !all processes switch back to previous value
       accept = accept - 1                                    !vector is 'de-accepted' since reverting to the value in the previous generation
    else
#ifdef MPI
       !root process shares newly-created vector with other processes

       !TODO Implement some sort of configure switch in future.  This should be
       !call MPI_Type_create_f90_real(precision(1.0_dp), range(1.0_dp), mpi_dp, ierror)
       !call MPI_Bcast(newvector, run_params%D, mpi_dp, root, MPI_COMM_WORLD, ierror)
       !This cannot be done at present, due to an MPICH bug: http://trac.mpich.org/projects/mpich/ticket/1769

       call MPI_Bcast(newvector, run_params%D, mpi_double_precision, root, MPI_COMM_WORLD, ierror)
       call MPI_Barrier(MPI_COMM_WORLD,ierror)
       allvecs(n,:) = newvector
#else
       !only one process, so just replace value in allvecs with newly-created vector
       allvecs(n,:) = Xnew%vectors(m,:)
#endif
    end if


  end subroutine replace_vector


end module selection
