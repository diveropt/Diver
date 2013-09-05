module init

use detypes
use deutils
use mutation, only: init_FjDE    !initializes scale factors for jDE
use crossover, only: init_CrjDE  !initializes crossovers for jDE
use selection, only: roundvector, replace_generation

implicit none

#ifdef USEMPI
  include 'mpif.h'
#endif

private
public param_assign, initialize, init_all_random_seeds

contains 

  !Assign parameter values (defaults if not specified) to run_params and print DE parameter values to screen

  subroutine param_assign(run_params, lowerbounds, upperbounds, nDerived, discrete, partitionDiscrete, maxciv, maxgen, NP, F, Cr, lambda, &
                          current, expon, bndry, jDE, removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount)

    type(codeparams), intent(out) :: run_params 
    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    integer, intent(in), optional  :: nDerived	 		!input number of derived quantities to output
    integer, dimension(:), intent(in), optional :: discrete     !lists all discrete dimensions of parameter space
    logical, intent(in), optional  :: partitionDiscrete         !split the population evenly amongst discrete parameters and evolve separately
    integer, intent(in), optional  :: maxciv 			!maximum number of civilisations
    integer, intent(in), optional  :: maxgen 			!maximum number of generations per civilisation
    integer, intent(in), optional  :: NP 			!population size (individuals per generation)
    real(dp), dimension(:), intent(in), optional :: F		!scale factor(s).  Note that this must be entered as an array.
    real(dp), intent(in), optional :: Cr 			!crossover factor
    real(dp), intent(in), optional :: lambda 			!mixing factor between best and rand/current
    logical, intent(in), optional  :: current 			!use current vector for mutation
    logical, intent(in), optional  :: expon 			!use exponential crossover
    integer, intent(in), optional  :: bndry                     !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional  :: jDE                       !use self-adaptive DE 
    logical, intent(in), optional  :: removeDuplicates          !weed out duplicate vectors within a single generation
    logical, intent(in), optional  :: doBayesian                !calculate log evidence and posterior weightings
    real(dp), intent(in), optional :: maxNodePop                !population at which node is partitioned in binary space partitioning for posterior
    real(dp), intent(in), optional :: Ztolerance		!input tolerance in log-evidence
    integer, intent(in), optional  :: savecount			!save progress every savecount generations

    integer :: mpiprocs, mpirank, ierror                        !number of processes running, rank of current process, error code
    character (len=30) :: DEstrategy                    	!for printing mutation/crossover DE strategy

    integer, allocatable :: num_discrete_vals(:)
    integer :: i, discrete_index

#ifdef USEMPI
    call MPI_Comm_size(MPI_COMM_WORLD, mpiprocs, ierror)  !gives the total num of processes. If no MPI, set to 1.
    call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierror)   !rank of the current process. If no MPI, set to 0.
#else
    mpiprocs = 1
    mpirank = 0
#endif

    if (mpirank .eq. 0) then
       write (*,*) 
       write (*,*) '============================='
       write (*,*) ' ******** Begin DE ********* '
       write (*,*) 'DEIS v0.1'
       write (*,*) 'Copyright Elinore Roebber and Pat Scott, 2013'
       write (*,*) 
    end if

    run_params%mpirank = mpirank
    run_params%D=size(lowerbounds)

    if (size(upperbounds) .ne. run_params%D) call quit_de('ERROR: parameter space bounds do not have the same dimensionality')
    if (any(lowerbounds .ge. upperbounds)) call quit_de('ERROR: invalid parameter space bounds.')

    allocate(run_params%lowerbounds(run_params%D), run_params%upperbounds(run_params%D))
    run_params%lowerbounds = lowerbounds
    run_params%upperbounds = upperbounds

    if (present(nDerived)) then 
      call setIfPositive_int(nDerived, run_params%D_derived, 'nDerived') !FIXME: this prevents setting nDerived=0
    else
       run_params%D_derived = 0					!default is no derived quantities
    end if

    if (present(maxciv)) then
      call setIfPositive_int(maxciv, run_params%numciv, 'maxciv')
    else
       run_params%numciv = 2000					!arbitrary default value for numciv
    end if

    if (present(maxgen)) then
      call setIfPositive_int(maxgen, run_params%numgen, 'maxgen')
    else 
       run_params%numgen = 300					!arbitrary default value for numgen
    end if

    if (present(doBayesian)) then 
       run_params%calcZ = doBayesian
    else
       run_params%calcZ = .false.				!default is not to do Bayesian stuff
    end if

    if (present(maxNodePop)) then 
       call setIfPositive_real(maxNodePop, run_params%maxNodePop, 'maxNodePop')
    else
       run_params%maxNodePop = 1.9_dp				!default for maxNodePop
    end if

    if (present(Ztolerance)) then 
       call setIfPositive_real(Ztolerance, run_params%tol, 'Ztolerance')
    else
       run_params%tol = 0.01_dp					!default for tolerance
    end if

    if (present(savecount)) then 
      call setIfPositive_int(savecount, run_params%savefreq, 'savecount')
    else
       run_params%savefreq = 1					!default for tolerance counter
    end if

    if (present(jDE)) then
       run_params%DE%jDE = jDE
    else
       run_params%DE%jDE = .false.                              !default is original DE
    endif

    if (present(bndry)) then
       run_params%DE%bconstrain = bndry
    else
       run_params%DE%bconstrain = 1				!default brick wall boundary constraints
    end if

    !set values for F, Cr, lambda, current, expon for self-adaptive rand/1/bin DE or regular DE
    jDEset: if (run_params%DE%jDE) then

       if (present(F)) write (*,*) 'WARNING: value set for F not used during jDE run'
       if (present(Cr)) write (*,*) 'WARNING: value set for Cr not used during jDE run'
       if (present(lambda)) write (*,*) 'WARNING: lambda not used during jDE run'
       if (present(current)) write (*,*) 'WARNING: current not used during jDE run'
       if (present(expon)) write (*,*) 'WARNING: jDE uses binary crossover'

       !the {F_i} and {Cr_i} are kept as part of the population (X%FjDE and X%CrjDE) with single variables (trialF and trialCr) 
       !in the main subroutine of the program for the trial parameters.
 
       if (present(NP)) then
          if (NP .ge. 4) then 					!required for picking unique vectors during mutation
             run_params%DE%NP = NP
          else
             write (*,*) 'WARNING: NP specified is too small. Using smallest permitted NP.'
             run_params%DE%NP = 4
          end if
       else
          run_params%DE%NP = maxval( [10*run_params%D, 4] )	!conservative rule-of-thumb choice 
       end if

       if (mod(run_params%DE%NP, mpiprocs) .ne. 0) then         !population chunks must be equally sized
          write (*,*) 'WARNING: increasing population size to be a multiple of the number of processes.' 
          run_params%DE%NP = run_params%DE%NP - mod(run_params%DE%NP, mpiprocs) + mpiprocs
       end if


       if (present(removeDuplicates)) then
          run_params%DE%removeDuplicates = removeDuplicates
       else if (mpiprocs .gt. 1) then 
          run_params%DE%removeDuplicates = .true.               !FIXME: MPI seems to lead to duplicates in first generation
       else
          run_params%DE%removeDuplicates = .false.              !with jDE mutation, duplicates are rare (CHECK THIS)
       end if

       run_params%DE%Fsize = 0
       run_params%DE%lambda = 0.
       run_params%DE%current = .false.
       run_params%DE%expon = .false.
       
       DEstrategy = 'self-adaptive rand/1/bin (jDE)'            !for printing to the screen

    else                                                        !not using jDE.  Initialize for normal DE

       if (present(F)) then
          if (any(F .le. 0.0_dp)) write (*,*) 'WARNING: some elements of F are 0 or negative. DE may not converge properly.'
          if (any(F .ge. 1.0_dp)) write (*,*) 'WARNING: some elements of F are 1 or greater. DE may not converge properly.'
          run_params%DE%Fsize = size(F)
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = F
       else
          run_params%DE%Fsize = 1
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = (/0.7_dp/) 				!rule of thumb: 0.4<F<1.0
       end if

       if (present(NP)) then                                    
          if (NP .ge. (2*run_params%DE%Fsize + 2)) then 	!required for picking unique vectors during mutation
                run_params%DE%NP = NP
          else !nb: if current=true and/or lambda=0, NP could be smaller, but it's a bad idea, so not implemented
             write (*,*) 'WARNING: NP specified is too small. Using smallest permitted NP.'
             run_params%DE%NP = 2*run_params%DE%Fsize + 3
          end if
       else
          run_params%DE%NP = maxval( [10*run_params%D, 2*run_params%DE%Fsize + 2] )	!conservative rule-of-thumb choice 
       end if

       if (mod(run_params%DE%NP, mpiprocs) .ne. 0) then
          write (*,*) 'WARNING: increasing population size to be a multiple of the number of processes.' 
          run_params%DE%NP = run_params%DE%NP - mod(run_params%DE%NP, mpiprocs) + mpiprocs
       end if

       if (present(Cr)) then  
          if (Cr .lt. 0.0_dp) then
             write (*,*) 'WARNING: Cr < 0. Using Cr = 0.' 	!although Cr<0 is functionally equivalent to Cr=0
             run_params%DE%Cr = 0.0_dp
          else if (Cr .gt. 1.0_dp) then
             write (*,*) 'WARNING: Cr > 1. Using Cr = 1.' 	!although Cr>1 is functionally equivalent to Cr=1
             run_params%DE%Cr = 1.0_dp
          else
             run_params%DE%Cr = Cr
          end if
       else
          run_params%DE%Cr = 0.9_dp 
       end if

       if (present(lambda)) then
          if (lambda .lt. 0.0_dp) write (*,*) 'WARNING: lambda < 0. DE may not converge properly.'
          if (lambda .gt. 1.0_dp) write (*,*) 'WARNING: lambda > 1. DE may not converge properly.'
          run_params%DE%lambda = lambda
       else
          run_params%DE%lambda = 0.0_dp     			!default rand/1/bin
       end if

       if (present(current)) then 
          run_params%DE%current = current
       else
          run_params%DE%current = .false. 			!default rand/1/bin
       end if
       
       if (present(expon)) then
          run_params%DE%expon = expon
       else
          run_params%DE%expon = .false.     			!default rand/1/bin
       end if

       if (present(removeDuplicates)) then
          run_params%DE%removeDuplicates = removeDuplicates
       else if (run_params%DE%current) then
          run_params%DE%removeDuplicates = .false.              !with current mutation, duplicate are rare (!FIXME CHECK THIS)
       else
          run_params%DE%removeDuplicates = .true.
       end if

       !for printing the parameter choice, DE mutation/crossover strategy, and boundary constraints to screen
       if (run_params%DE%lambda .eq. 0.0_dp) then  		!mutation strategy
          if (run_params%DE%current) then
             DEstrategy = 'current/'
          else
             DEstrategy = 'rand/'
          end if
       else if (run_params%DE%lambda .eq. 1.0_dp) then
          DEstrategy = 'best/'
       else 
          if (run_params%DE%current) then
             DEstrategy = 'current-to-best/'
          else
             DEstrategy = 'rand-to-best/'
          end if
       end if

       DEstrategy = trim(DEstrategy)//int_to_string(run_params%DE%Fsize)
       
       if(run_params%DE%expon) then                  		!crossover strategy
          DEstrategy = trim(DEstrategy)//'/exp'
       else
          DEstrategy = trim(DEstrategy)//'/bin'
       end if
    endif jDEset

    if (present(discrete)) then
       if (any(discrete .gt. run_params%D) .or. any(discrete .lt. 1)) then     
          call quit_de('ERROR: Discrete dimensions specified must not be < 1 or > '//trim(int_to_string(run_params%D)))
       end if     
       run_params%D_discrete = size(discrete)
       !Also check that discrete dimensions are not doubly-specified (will crash the partioned case if so, just sloppy otherwise.)
       do i = 1, run_params%D_discrete
          if (count(discrete .eq. discrete(i)) .ne. 1) then
             call quit_de('ERROR: Discrete dimension '//trim(int_to_string(discrete(i)))//'listed multiple times in call to run_de.')
          endif
       enddo
       allocate(run_params%discrete(run_params%D_discrete))
       run_params%discrete = discrete
    else
       run_params%D_discrete = 0
       allocate(run_params%discrete(0))
    end if

    if (present(partitionDiscrete)) then
       run_params%partitionDiscrete = partitionDiscrete
       if (partitionDiscrete) then
          if (run_params%D_discrete .eq. 0) then
             write(*,*) 'WARNING: keyword partitionDiscrete ignored because no discrete parameters were indicated.'
             run_params%partitionDiscrete = .false.
          else
             !Determine on a scale of how many individuals each discrete parameter should repeat
             allocate(run_params%repeat_scales(run_params%D_discrete), num_discrete_vals(run_params%D_discrete))
             do i = 1, run_params%D_discrete
                num_discrete_vals(i) = nint(run_params%upperbounds(run_params%discrete(i)) - &
                                            run_params%lowerbounds(run_params%discrete(i)) + 1 )
                run_params%repeat_scales(i) = run_params%DE%NP
                discrete_index = i
                do while (discrete_index .ne. 1)
                   discrete_index = discrete_index - 1
                   run_params%repeat_scales(i) = nint( dble(run_params%repeat_scales(i)) / dble(num_discrete_vals(discrete_index)) )
                enddo
             enddo
             if ( mod(run_params%DE%NP, product(num_discrete_vals)) .ne. 0) then
                call quit_de('ERROR: partitionDiscrete = true requires that NP must divide up evenly into    the implied number of sub-populations.')
             else
                !Work out how many individuals should be in each discrete partition of the population (ie each subpopulation)
                run_params%subpopNP = run_params%DE%NP / product(num_discrete_vals)
             endif
             deallocate(num_discrete_vals)
          endif
       endif
    else
      run_params%partitionDiscrete = .false.
      run_params%subpopNP = run_params%DE%NP
    endif

    !split up work over multiple processes for MPI (if no mpi, single process does all the work)
    run_params%mpipopchunk = run_params%DE%NP/mpiprocs

    !print feedback about strategy choice, values of parameters to the screen
    if (run_params%mpirank .eq. 0) then
       write (*,*) DEstrategy
       if (size(run_params%discrete) .gt. 0) write (*,*) 'Discrete dimensions:', run_params%discrete
       write (*,*) 'Parameters:'
       write (*,*) ' NP = ', trim(int_to_string(run_params%DE%NP))
       if ((run_params%DE%lambda .ne. 1.0_dp) .and. (run_params%DE%lambda .ne. 0.0_dp)) then
          write (*,'(A10, F6.4)') ' lambda = ', run_params%DE%lambda
       endif
       if (.not. run_params%DE%jDE) then                   
          write (*,'(A5, F6.4)') ' F = ', run_params%DE%F
          write (*,'(A6, F6.4)') ' Cr = ', run_params%DE%Cr 
       endif

       write (*,*) 'Number of processes: ', trim(int_to_string(mpiprocs))

       select case (run_params%DE%bconstrain)                       !boundary constraint choice
       case (1) 
          write (*,*) 'Brick wall boundary constraints'
       case (2)
          write (*,*) 'Random re-initialization boundary constraints'
       case (3)
          write (*,*) 'Reflective boundary constraints'
       case default
          write (*,*) 'WARNING: Invalid value entered for bndry. Boundary constraints not enforced.'
       end select
    end if

  end subroutine param_assign


  !set parameter only if value is positive (real parameter)
  subroutine setIfPositive_real(invar, outvar, string)

    real(dp) :: invar, outvar
    character(LEN=*) :: string

    if (invar .le. 0.0_dp) then
       call quit_de('ERROR: '//string//' cannot be negative.')
    else
      outvar = invar
    endif

  end subroutine setIfPositive_real


  !set parameter only if value is positive (integer parameter)
  subroutine setIfPositive_int(invar, outvar, string)

    integer :: invar, outvar
    character(LEN=*) :: string

    if (invar .le. 0) then
      call quit_de('ERROR: '//string//' cannot be negative.')
    else
      outvar = invar
    endif

  end subroutine setIfPositive_int


  !initializes first generation of target vectors
  subroutine initialize(X, Xnew, run_params, fcall, func, quit) 

    type(population), intent(inout) :: X
    type(population), intent(inout) :: Xnew
    type(codeparams), intent(in) :: run_params
    integer, intent(inout) :: fcall
    logical, intent(inout) :: quit
    real(dp), external :: func
    integer :: n, m, i, discrete_index, accept=0

    X%multiplicities = 1.0_dp !Initialise to 1 in case posteriors are not calculated

    if (run_params%DE%jDE) then !initialize population of F and Cr parameters
       Xnew%FjDE = init_FjDE(run_params, run_params%mpipopchunk)
       Xnew%CrjDE = init_CrjDE(run_params, run_params%mpipopchunk)
    end if

    !loop over the vectors belonging to each population chunk
       do m=1,run_params%mpipopchunk

          if (verbose .or. run_params%partitionDiscrete) then 
             n = run_params%mpipopchunk*run_params%mpirank + m !true population index (equal to m if no mpi)
          endif

          call random_number(Xnew%vectors(m,:))

          do i = 1, run_params%D

             if (run_params%partitionDiscrete .and. any(run_params%discrete .eq. i)) then
                !This is a discrete parameter that needs to be partitioned, and therefore set rather than chosen randomly

                !Determine the index of the discrete parameter (first discrete param, second, etc)
                do discrete_index = 1, run_params%D_discrete
                   if (run_params%discrete(discrete_index) .eq. i) exit
                enddo

                !Set the value of this index for this individual
                Xnew%vectors(m,i) = anint( dble(mod(n-1,run_params%repeat_scales(discrete_index))) / &
                                           dble(run_params%repeat_scales(discrete_index)) * &
                                           (run_params%upperbounds(i) - run_params%lowerbounds(i) + 1) &
                                           + run_params%lowerbounds(i) - 0.5_dp + 100._dp*epsilon(0.5_dp) ) 

             else
                !This is a normal parameter that needs to be chosen randomly
                Xnew%vectors(m,i) = Xnew%vectors(m,i)*(run_params%upperbounds(i) - run_params%lowerbounds(i)) + run_params%lowerbounds(i)

             endif

          enddo

          Xnew%vectors_and_derived(m,:run_params%D) = roundvector(Xnew%vectors(m,:), run_params)
          Xnew%values(m) = func(Xnew%vectors_and_derived(m,:), fcall, quit, .true.)

          if (verbose) then
             if (run_params%DE%jDE) then
                write (*,*) n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m), '|', Xnew%FjDE(m), Xnew%CrjDE(m)
             else
                write (*,*) n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m)
             end if
          end if

       end do

       call replace_generation(X, Xnew, run_params, func, fcall, quit, accept, init=.true.)
    
  end subroutine initialize


  !based on init_random_seed below. Calls init_random_seed for master process,
  !then (if using MPI) generates new seeds for secondary processes and distributes them
  subroutine init_all_random_seeds(nprocs, mpirank)
    integer, intent(in) :: nprocs         !number of processes that need seeds
    integer, intent(in) :: mpirank 
    integer :: i, n, clock, ierror
    real, dimension(:), allocatable :: rand
    integer, dimension(:,:), allocatable :: allseeds
    integer, dimension(:), allocatable :: seed

    !initialize master seed
    if (mpirank .eq. 0) call init_random_seed()

    !initialize random seeds for secondary processes
#ifdef USEMPI
    call random_seed(size = n)
    allocate(allseeds(n, nprocs))
    allocate(seed(n))

    if (mpirank .eq. 0) then                     !master process
       allocate(rand(nprocs))
       call random_number(rand)
       call system_clock(count=clock)
       rand = clock*(1 + rand)

       allseeds = spread(37*(/ (i - 1, i = 1, n) /), dim=2, ncopies=nprocs)
       !equivalent to seed = 37 * (/ (i - 1, i = 1, n) /), but with nprocs copies (n x nprocs)

       allseeds = allseeds + spread(int(rand), dim=1, ncopies=n)
       !equivalent to seed = seed + rand, but rand is stretched to n x nprocs

       deallocate(rand)
    end if

    call MPI_SCATTER(allseeds, n, MPI_INTEGER, seed, n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    if (mpirank .ne. 0) call random_seed(put = seed)

    deallocate(seed, allseeds)
#endif

  end subroutine init_all_random_seeds


  !Yanked from the gfortran documentation
  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)

    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)
  END SUBROUTINE


end module init
