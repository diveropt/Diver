module init

use detypes
use deutils
use mutation, only: init_FjDE, init_lambdajDE    !initializes scale factors for jDE
use crossover, only: init_CrjDE                  !initializes crossovers for jDE
use selection, only: replace_generation

implicit none

#ifdef MPI
  include 'mpif.h'
#endif

private
public param_assign, initialize, init_all_random_seeds

character (len=*), parameter :: version_number = "1.0.5"

contains

  !Assign parameter values (defaults if not specified) to run_params and print DE parameter values to screen

  subroutine param_assign(run_params, lowerbounds, upperbounds, nDerived, discrete, partitionDiscrete, maxciv,     &
                          maxgen, NP, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, convthresh, convsteps, &
                          removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount, outputSamples,          &
                          init_population_strategy, discard_unfit_points, max_initialisation_attempts,             &
                          max_acceptable_value, seed, context, verbose)

    use iso_c_binding, only: C_NULL_PTR, c_ptr

    type(codeparams), intent(out) :: run_params
    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds      !boundaries of parameter space
    integer, intent(in), optional  :: nDerived                          !input number of derived quantities to output
    integer, dimension(:), intent(in), optional :: discrete             !lists all discrete dimensions of parameter space
    logical, intent(in), optional  :: partitionDiscrete                 !split the population evenly amongst discrete parameters and evolve separately
    integer, intent(in), optional  :: maxciv                            !maximum number of civilisations
    integer, intent(in), optional  :: maxgen                            !maximum number of generations per civilisation
    integer, intent(in), optional  :: NP                                !population size (individuals per generation)
    real(dp), dimension(:), intent(in), optional :: F                   !scale factor(s).  Note that this must be entered as an array.
    real(dp), intent(in), optional :: Cr                                !crossover factor
    real(dp), intent(in), optional :: lambda                            !mixing factor between best and rand/current
    logical, intent(in), optional  :: current                           !use current vector for mutation
    logical, intent(in), optional  :: expon                             !use exponential crossover
    integer, intent(in), optional  :: bndry                             !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional  :: jDE                               !use self-adaptive DE
    logical, intent(in), optional  :: lambdajDE                         !use self-adaptive DE with rand-to-best mutation strategy
    real(dp), intent(in), optional :: convthresh                        !threshold for generation-level convergence
    integer, intent(in), optional  :: convsteps                         !number of steps to smooth over when checking convergence
    logical, intent(in), optional  :: removeDuplicates                  !weed out duplicate vectors within a single generation
    logical, intent(in), optional  :: doBayesian                        !calculate log evidence and posterior weightings
    real(dp), intent(in), optional :: maxNodePop                        !population at which node is partitioned in binary space partitioning for posterior
    real(dp), intent(in), optional :: Ztolerance                        !input tolerance in log-evidence
    integer, intent(in), optional  :: savecount                         !save progress every savecount generations
    logical, intent(in), optional  :: outputSamples                     !write samples as output
    integer, intent(in), optional  :: init_population_strategy          !initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
    logical, intent(in), optional  :: discard_unfit_points              !recalculate any trial vector whose fitness is above max_acceptable_value
    integer, intent(in), optional  :: max_initialisation_attempts       !maximum number of times to try to find a valid vector for each slot in the initial population.
    real(dp), intent(in), optional :: max_acceptable_value              !maximum fitness to accept for the initial generation if init_population_strategy > 0. Also applies to later generations if discard_unfit_points = .true.
    integer, intent(in), optional  :: seed                              !base seed for random number generation; non-positive or absent means seed from the system clock
    type(c_ptr), intent(inout), optional  :: context                    !context pointer/integer, used for passing info from the caller to likelihood/prior
    integer, intent(in), optional  :: verbose                           !how much info to print to screen: 0-quiet, 1-basic info, 2-civ info, 3+ everything

    integer :: mpiprocs, mpirank, ierror                                !number of processes running, rank of current process, error code
    character (len=70) :: DEstrategy                                    !for printing mutation/crossover DE strategy

    integer, allocatable :: num_discrete_vals(:)
    integer :: i, discrete_index

#ifdef MPI
    call MPI_Comm_size(MPI_COMM_WORLD, mpiprocs, ierror)  !gives the total num of processes. If no MPI, set to 1.
    call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierror)   !rank of the current process. If no MPI, set to 0.
#else
    mpiprocs = 1
    mpirank = 0
#endif

    run_params%mpiprocs = mpiprocs
    run_params%mpirank = mpirank

    !save the seed for the random number generator
    if (present(seed)) then
       if (seed .le. 0) then
         run_params%seed = -1
       else
         run_params%seed = seed
       endif
    else
       run_params%seed = -1
    endif

    !default level of output is to print errors, most warnings, info about the program & final results
    call setIfNonNegative_int('verbose', run_params%verbose, 1, invar=verbose)

    !to ensure that lines are only printed by one process at a time, the non-root processes are set to negative integers
    !(most printing is done by the root process alone, but information about specific population members comes from individual processes)
    if (run_params%mpirank .ne. 0) then
       run_params%verbose = -run_params%verbose
    end if

    if (run_params%verbose .ge. 1) then
       write (*,*)
       write (*,*) '============================================='
       write (*,*) '************ Begin DE Sampling **************'
       write (*,*) 'Diver v'//trim(version_number)
       write (*,*) 'Copyright Elinore Roebber and Pat Scott, 2013'
       write (*,*) '============================================='
       write (*,*)
    end if

    run_params%D=size(lowerbounds)

    if (size(upperbounds) .ne. run_params%D) call quit_de('ERROR: parameter space bounds do not have the same dimensionality')
    if (any(lowerbounds .ge. upperbounds)) call quit_de('ERROR: invalid parameter space bounds.')

    allocate(run_params%lowerbounds(run_params%D), run_params%upperbounds(run_params%D))
    run_params%lowerbounds = lowerbounds
    run_params%upperbounds = upperbounds

    call setIfNonNegative_int('nDerived', run_params%D_derived, 0, invar=nDerived) !default is no derived quantities
    call setIfPositive_int('maxgen', run_params%numgen, 300, invar=maxgen)
    call setIfPositive_real('convthresh', run_params%convthresh, 1e-3_dp, invar=convthresh)

    if (run_params%convthresh .ge. 1_dp) then                                      !fractional improvement always le 1
       call quit_de('ERROR: threshold for convergence (convthresh) must be < 1')
    endif

    call setIfPositive_int('convsteps', run_params%convsteps, 10, invar=convsteps)
    allocate(run_params%improvements(run_params%convsteps))

    call set_logical(run_params%calcZ, .false., invar=doBayesian)                  !default is not to do Bayesian stuff

    if (run_params%calcZ) then
       call setIfPositive_real('maxNodePop', run_params%maxNodePop, 1.9_dp, invar=maxNodePop)
       call setIfPositive_real('Ztolerance', run_params%tol, 0.01_dp, invar=Ztolerance)
       call setIfPositive_int('maxciv', run_params%numciv, 2000, invar=maxciv)
    else
       !when not doing evidence calculations, no need for many civilizations
       call setIfPositive_int('maxciv', run_params%numciv, 1, invar=maxciv)
    endif

    call setIfPositive_int('savecount', run_params%savefreq, 1, invar=savecount)   !tolerance counter

    call set_logical(run_params%DE%lambdajDE, .true., invar=lambdajDE)             !default is to turn on lambdajDE
    call set_logical(run_params%DE%jDE, .true., invar=jDE)                         !default is to turn on jDE
    if (run_params%DE%lambdajDE) then
       if (.not. run_params%DE%jDE) call quit_de('ERROR: cannot have lambdajDE=true with jDE=false')
    end if

    call setIfNonNegative_int('bndry', run_params%DE%bconstrain, 1, invar=bndry)   !default brick wall boundary constraints
    if (run_params%DE%bconstrain .gt. 3) then
          if (run_params%verbose .ge. 1) then
             write (*,*) 'Legal values for bndry (enforces boundary constraints):'
             write (*,*) ' 0: Not enforced'
             write (*,*) ' 1: Brick wall'
             write (*,*) ' 2: Random re-initialization'
             write (*,*) ' 3: Reflection'
          end if
          call quit_de('ERROR: Invalid value entered for bndry.')
       end if


    !set values for F, Cr, lambda, current, expon for self-adaptive rand/1/bin DE or regular DE
    jDEset: if (run_params%DE%jDE) then

       if (run_params%verbose .ge. 1) then
          if (present(F)) write (*,*) 'WARNING: value set for F not used during jDE run'
          if (present(Cr)) write (*,*) 'WARNING: value set for Cr not used during jDE run'
          !the {F_i} and {Cr_i} are kept as part of the population (X%FjDE and X%CrjDE) with single variables (trialF and trialCr)
          !in the main subroutine of the program for the trial parameters.
          if (present(current)) then
             if (current) write (*,*) 'WARNING: current not used during jDE run'
          end if
          if (present(expon)) then
             if (expon) write (*,*) 'WARNING: jDE uses binary crossover. Value set for expon will be ignored.'
          end if
       end if

       run_params%DE%current = .false.
       run_params%DE%expon = .false.

       !With jDE, Fsize becomes size(population%FjDE)/NP, not size(run_params%DE%F)
       run_params%DE%Fsize = 1
       allocate(run_params%DE%F(run_params%DE%Fsize))
       run_params%DE%F = 0.

       !Options for lambda:
       if (run_params%DE%lambdajDE) then
          run_params%DE%lambda = 0.0_dp
          if (present(lambda) .and. (run_params%verbose .ge. 1) ) then
             write (*,*) 'WARNING: value set for lambda not used during lambdajDE run'
          end if
       else if (present(lambda)) then
          if (run_params%verbose .ge. 1) then
             if (lambda .lt. 0.0_dp) write (*,*) 'WARNING: lambda < 0. DE may not converge properly.'
             if (lambda .gt. 1.0_dp) write (*,*) 'WARNING: lambda > 1. DE may not converge properly.'
          end if
          run_params%DE%lambda = lambda
       else
          run_params%DE%lambda = 0.0_dp                 !default rand/1/bin
       end if

       ! Summarize mutation strategy: for printing to screen
       if (run_params%DE%lambdajDE) then
          DEstrategy = 'self-adaptive rand-to-best/1/bin (lambdajDE)'
       else if (run_params%DE%lambda .eq. 0.0_dp) then
          DEstrategy = 'self-adaptive rand/1/bin (jDE)'
       else if (run_params%DE%lambda .eq. 1.0_dp) then
          DEstrategy = 'self-adaptive best/1/bin (jDE with additional fixed lambda)'
       else
          DEstrategy = 'self-adaptive rand-to-best/1/bin (jDE with additional fixed lambda)'
       end if

    else
       !not using jDE.  Initialize for normal DE
       if (present(F)) then
          if (any(F .le. 0.0_dp) .and. (run_params%verbose .ge. 1) ) then
             write (*,*) 'WARNING: some elements of F are 0 or negative. DE may not converge properly.'
          else if (any(F .ge. 1.0_dp) .and. (run_params%verbose .ge. 1) ) then
             write (*,*) 'WARNING: some elements of F are 1 or greater. DE may not converge properly.'
          end if
          run_params%DE%Fsize = size(F)
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = F
       else
          run_params%DE%Fsize = 1
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = (/0.7_dp/)                 !rule of thumb: 0.4<F<1.0
       end if

       if (present(Cr)) then
          if (Cr .lt. 0.0_dp) then
             if (run_params%verbose .ge. 1)  write (*,*) 'WARNING: Cr < 0. Using Cr = 0.' !although Cr<0 is functionally equivalent to Cr=0
             run_params%DE%Cr = 0.0_dp
          else if (Cr .gt. 1.0_dp) then
             if (run_params%verbose .ge. 1) write (*,*) 'WARNING: Cr > 1. Using Cr = 1.' !although Cr>1 is functionally equivalent to Cr=1
             run_params%DE%Cr = 1.0_dp
          else
             run_params%DE%Cr = Cr
          end if
       else
          run_params%DE%Cr = 0.9_dp
       end if

       if (present(lambda)) then
          if (run_params%verbose .ge. 1) then
             if (lambda .lt. 0.0_dp) write (*,*) 'WARNING: lambda < 0. DE may not converge properly.'
             if (lambda .gt. 1.0_dp) write (*,*) 'WARNING: lambda > 1. DE may not converge properly.'
          end if
          run_params%DE%lambda = lambda
       else
          run_params%DE%lambda = 0.0_dp                 !default rand/1/bin
       end if


       !default rand/1/bin
       call set_logical(run_params%DE%current, .false., invar=current)
       call set_logical(run_params%DE%expon, .false., invar=expon)


       !for printing the parameter choice, DE mutation/crossover strategy, and boundary constraints to screen
       if (run_params%DE%lambda .eq. 0.0_dp) then          !mutation strategy
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

       if(run_params%DE%expon) then                          !crossover strategy
          DEstrategy = trim(DEstrategy)//'/exp'
       else
          DEstrategy = trim(DEstrategy)//'/bin'
       end if
    endif jDEset

    !Set whether or not to remove duplicates
    if ((run_params%DE%jDE .or. run_params%DE%current) .and. (mpiprocs .eq. 1)) then
       !for current mutation or jDE, natural duplicates are rare, but it can be a useful diagnostic for problems with MPI
       call set_logical(run_params%DE%removeDuplicates, .false., invar=removeDuplicates)
    else
       call set_logical(run_params%DE%removeDuplicates, .true., invar=removeDuplicates)
    end if

    !Set up array of discrete dimensions
    if (present(discrete)) then
       if (any(discrete .gt. run_params%D) .or. any(discrete .lt. 1)) then
          call quit_de('ERROR: Discrete dimensions specified must not be < 1 or > '// &
                       trim(int_to_string(run_params%D)))
       end if
       run_params%D_discrete = size(discrete)
       !Also check that discrete dimensions are not doubly-specified (will crash the partioned case if so, just sloppy otherwise.)
       do i = 1, run_params%D_discrete
          if (count(discrete .eq. discrete(i)) .ne. 1) then
             call quit_de('ERROR: Discrete dimension '//trim(int_to_string(discrete(i)))// &
                          'listed multiple times in call to run_de.')
          endif
       enddo
       allocate(run_params%discrete(run_params%D_discrete))
       run_params%discrete = discrete
    else
       run_params%D_discrete = 0
       allocate(run_params%discrete(0))
    end if

    !Setting NP: default is the larger of a conservative rule-of-thumb value or the minimum required for mutation
    call setIfPositive_int('NP', run_params%DE%NP, maxval((/10*run_params%D, 2*run_params%DE%Fsize + 3/)), invar=NP)

    if (run_params%DE%NP .lt. (2*run_params%DE%Fsize + 3)) then
       if (run_params%verbose .ge. 1) write (*,*) 'WARNING: NP specified is too small. Using smallest permitted NP...'
       run_params%DE%NP = 2*run_params%DE%Fsize + 3  !required to pick unique vectors during mutation
    end if

    !Check that population chunks are equally sized
    if (mod(run_params%DE%NP, mpiprocs) .ne. 0) then
       if (run_params%verbose .ge. 1) then
          write (*,*) 'WARNING: increasing population size to be a multiple of the number of processes...'
       end if
       run_params%DE%NP = run_params%DE%NP - mod(run_params%DE%NP, mpiprocs) + mpiprocs
    end if

    !Max sure we don't saturate the print statements
    if (run_params%DE%NP .gt. 99999999)  &
      call quit_de('ERROR: NP cannot be greater than 99999999.  Your requested run requires NP = ' &
                   //trim(int_to_string(run_params%DE%NP))//'.  Please modify NP in your call to run_de.')

    !split up work over multiple processes for MPI (if no mpi, single process does all the work)
    run_params%mpipopchunk = run_params%DE%NP/mpiprocs

    !Set up subpopulations for discrete partitioning
    call set_logical(run_params%partitionDiscrete, .false., invar=partitionDiscrete)
    if (run_params%partitionDiscrete) then
       if (run_params%D_discrete .eq. 0) then
          if (run_params%verbose .ge. 1) then
             write(*,*) 'WARNING: keyword partitionDiscrete ignored because no discrete parameters were indicated.'
          end if
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
                run_params%repeat_scales(i) = nint( dble(run_params%repeat_scales(i)) / &
                                                 dble(num_discrete_vals(discrete_index)) )
             enddo
          enddo
          if ( mod(run_params%DE%NP, product(num_discrete_vals)) .ne. 0) then
             if (run_params%verbose .ge. 1) then
                write(*,*) 'Population size (NP): '//trim(int_to_string(run_params%DE%NP))
                write(*,*) 'Range(s) of discrete parameter(s):', num_discrete_vals
                write(*,*) 'Implied number of sub-populations (product of the ranges): '// &
                              trim(int_to_string(product(num_discrete_vals)))
             end if
             call quit_de('ERROR: partitionDiscrete = true requires that NP must divide up'// &
                           ' evenly into the implied number of sub-populations.')
          else
             !Work out how many individuals should be in each discrete partition of the population (ie each subpopulation)
             run_params%subpopNP = run_params%DE%NP / product(num_discrete_vals)
             if (run_params%subpopNP .lt. (2*run_params%DE%Fsize + 3) ) then
                if (run_params%verbose .ge. 1) then
                   write(*,*) 'Population size (NP): '//trim(int_to_string(run_params%DE%NP))
                   write(*,*) 'Implied number of sub-populations: '//trim(int_to_string(product(num_discrete_vals)))
                   write(*,*) 'Sub-population size: '//trim(int_to_string(run_params%subpopNP))
                   write(*,*) 'Minimum allowed sub-population size: '//trim(int_to_string(2*run_params%DE%Fsize + 3))
                   write(*,*) 'Please increase NP to at least '// &
                           trim(int_to_string(product(num_discrete_vals)*(2*run_params%DE%Fsize + 3)))
                end if
                call quit_de('ERROR: NP is too small to partition into the desired number of sub-populations.')
             end if
          endif
          deallocate(num_discrete_vals)
       endif
    else
       run_params%subpopNP = run_params%DE%NP !no partitioning, so the subpopulation is the same as the whole population
    endif

    !Default is to output parameter samples.
    call set_logical(run_params%outputSamples, .true., invar=outputSamples)

    !Just set up a dummy null context pointer if it happens to be missing
    if (present(context)) then
      run_params%context = context
    else
      run_params%context = C_NULL_PTR
    endif

    !Default is not to demand valid points in the initial generation.
    call setIfNonNegative_int("init_population_strategy", run_params%init_population_strategy, 0, invar=init_population_strategy)

    !Default is not to demand trial vectors to pass fitness threshold.
    call set_logical(run_params%discard_unfit_points, .false., invar=discard_unfit_points)

    !Default is to allow 10,000 initialisation attempts in cases where valid points are preferred for the initial generation.
    call setIfPositive_int("max_initialisation_attempts", run_params%max_initialisation_attempts, 10000, &
     invar=max_initialisation_attempts)

    !Default is to consider points with fitnesses less than 1e6 valid when preferring valid points for initial generation.
    call setIfPositive_real("max_acceptable_value", run_params%max_acceptable_value, 1d6, invar=max_acceptable_value)

    !Parameters have all been set. Now print feedback about strategy choice, values of parameters to the screen
    if (run_params%verbose .ge. 1) then
       write (*,*) DEstrategy
       if (size(run_params%discrete) .gt. 0) write (*,*) 'Discrete dimensions:', run_params%discrete
       write (*,*) 'Parameters:'
       write (*,*) ' NP = ', trim(int_to_string(run_params%DE%NP))
       if (.not. run_params%DE%lambdajDE) then               !print lambda to the screen when fixed, including to 0 and 1
          write (*,'(A10, F6.4)') ' lambda = ', run_params%DE%lambda
          if (.not. run_params%DE%jDE) then                  !print F and Cr to the screen when fixed
             write (*,'(A5, F6.4)') ' F = ', run_params%DE%F
             write (*,'(A6, F6.4)') ' Cr = ', run_params%DE%Cr
          endif
       end if

       write (*,*) 'Number of processes: ', trim(int_to_string(mpiprocs))

       select case (run_params%DE%bconstrain)                !boundary constraint choice
       case (0)
          write (*,*) 'Boundary constraints not enforced'
       case (1)
          write (*,*) 'Brick wall boundary constraints'
       case (2)
          write (*,*) 'Random re-initialization boundary constraints'
       case (3)
          write (*,*) 'Reflective boundary constraints'
       end select

       select case (run_params%init_population_strategy)     !initialisation strategy
       case (0)
          write (*,*) 'Validity of initial generation will not be enforced.'
       case (1)
          write (*,*) 'Will make', run_params%max_initialisation_attempts, 'attempts to find a point with value below', &
           run_params%max_acceptable_value
          write (*,*) 'when generating the initialial population.  After this, invalid points will be permitted.'
       case (2)
          write (*,*) 'Will get', run_params%max_initialisation_attempts, 'attempts to find a point with value below', &
           run_params%max_acceptable_value
          write (*,*) 'before accepting any individual for the initial generation; will halt if unsuccessful.'
       end select

       select case (run_params%discard_unfit_points) !trial vector strategy
       case (.false.)
          write (*,*) 'Validity of trial vectors will not be enforced.'
       case (.true.)
          write (*,*) 'New trial vectors will be generated until all have values below', &
           run_params%max_acceptable_value
       end select


    end if

  end subroutine param_assign


!TODO: add subroutine to set reals between bounds (inclusive) for Cr, lambda

  !check if parameter is present & set parameter only if value is positive (real parameter)
  subroutine setIfPositive_real(string, outvar, defaultvar, invar)

    character(LEN=*), intent(in) :: string
    real(dp), intent(out) :: outvar
    real(dp), intent(in) :: defaultvar
    real(dp), intent(in), optional :: invar

    if (present(invar)) then
       if (invar .le. 0.0_dp) then
          call quit_de('ERROR: '//string//' must be positive.')
       else
          outvar = invar
       endif
   else
      outvar = defaultvar
   end if

  end subroutine setIfPositive_real


  !set parameter only if value is not negative (integer parameter)
  subroutine setIfNonNegative_int(string, outvar, defaultvar, invar)

    character(LEN=*), intent(in) :: string
    integer, intent(out) :: outvar
    integer, intent(in) :: defaultvar
    integer, intent(in), optional :: invar

    if (present(invar)) then
       if (invar .ne. 0) then
          call setIfPositive_int(string, outvar, defaultvar, invar=invar)
       else
          outvar = invar
       endif
    else
       outvar = defaultvar
    end if

  end subroutine setIfNonNegative_int

  !set parameter only if value is positive; if not present, set to default (integer parameter)
  subroutine setIfPositive_int(string, outvar, defaultvar, invar)
    character(LEN=*), intent(in) :: string
    integer, intent(out) :: outvar
    integer, intent(in) :: defaultvar
    integer, intent(in), optional :: invar

    if (present(invar)) then
       if (invar .le. 0) then
          call quit_de('ERROR: '//string//' cannot be negative.')
       else
          outvar = invar
       endif
    else
       outvar = defaultvar
    endif

  end subroutine setIfPositive_int

  subroutine set_logical(outvar, defaultvar, invar)
    logical, intent(out) :: outvar
    logical, intent(in) :: defaultvar
    logical, intent(in), optional :: invar

    if (present(invar)) then
       outvar = invar
    else
       outvar = defaultvar
    end if
  end subroutine set_logical


  !initializes first generation of target vectors
  subroutine initialize(X, Xnew, run_params, func, fcall, quit, accept)

    type(population), intent(inout) :: X
    type(population), intent(inout) :: Xnew
    type(codeparams), intent(inout) :: run_params
    integer, intent(inout) :: fcall
    logical, intent(inout) :: quit
    procedure(MinusLogLikeFunc) :: func
    integer :: n, m, i, discrete_index, attempt_count, max_attempts, accept, fcall_this_gen

    fcall_this_gen = 0 !Initialise to zero so as not to mess up acceptance in subsequent civilisations.

    X%multiplicities = 1.0_dp !Initialise to 1 in case posteriors are not calculated

    if (run_params%DE%jDE) then                              !initialize population of F and Cr parameters
       Xnew%FjDE = init_FjDE(run_params%mpipopchunk)
       Xnew%CrjDE = init_CrjDE(run_params%mpipopchunk)
       if (run_params%DE%lambdajDE) then
          Xnew%lambdajDE = init_lambdajDE(run_params%mpipopchunk)
       end if
    end if

    !loop over the vectors belonging to each population chunk
    do m=1,run_params%mpipopchunk

       n = run_params%mpipopchunk*run_params%mpirank + m !true population index (equal to m if no mpi)

       !if init_population_strategy is not 0, try to find a valid individual to put in the initial population.
       attempt_count = 0
       max_attempts = merge(1, run_params%max_initialisation_attempts, run_params%init_population_strategy .eq. 0)
       do while (attempt_count .lt. max_attempts)

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
                Xnew%vectors(m,i) = Xnew%vectors(m,i)*(run_params%upperbounds(i) - &
                                    run_params%lowerbounds(i)) + run_params%lowerbounds(i)

             endif

          enddo

          Xnew%vectors_and_derived(m,:run_params%D) = roundvector(Xnew%vectors(m,:), run_params)

          Xnew%values(m) = func(Xnew%vectors_and_derived(m,:), fcall_this_gen, quit, .true., run_params%context)

          if (quit) call quit_de('ERROR: quit flag raised whilst initialising first generation.  Forcing hard quit.')

          if (Xnew%values(m) .lt. run_params%max_acceptable_value .or. run_params%init_population_strategy .eq. 0) exit

          attempt_count = attempt_count + 1

       enddo

       !crash if valid vectors have been demanded in the initial population but could not be found.
       if (run_params%init_population_strategy .ge. 2 .and. attempt_count .eq. max_attempts) then
         call quit_de('ERROR: init_population_strategy = 2 but could not find valid points within max_initialisation_attempts!')
       endif

       if (abs(run_params%verbose) .ge. 3) then
          if (run_params%DE%lambdajDE) then
             write (*,*) n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m), '|', &
                            Xnew%lambdajDE(m), Xnew%FjDE(m), Xnew%CrjDE(m)
          else if (run_params%DE%jDE) then
             write (*,*) n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m), '|', Xnew%FjDE(m), Xnew%CrjDE(m)
          else
             write (*,*) n, Xnew%vectors_and_derived(m,:), '->', Xnew%values(m)
          end if
       end if

    end do

    accept = run_params%mpipopchunk * run_params%mpipopchunk / fcall_this_gen
    fcall = fcall + fcall_this_gen

    call replace_generation(X, Xnew, run_params, func, fcall, quit, accept, init=.true.)

  end subroutine initialize


  !Generates a seed for master process, then (if using MPI) generates new seeds
  !for secondary processes and distributes them
  subroutine init_all_random_seeds(nprocs, mpirank, input_seed)
    integer, intent(in) :: nprocs         !number of processes that need seeds
    integer, intent(in) :: mpirank
    integer, intent(in) :: input_seed
    integer :: i, n, clock, ierror
    real(dp), dimension(:), allocatable :: rand
    integer, dimension(:,:), allocatable :: allseeds
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    !initialize master seed
    if (mpirank .eq. 0) then

      if (input_seed .gt. 0) then
        clock = input_seed
      else
        call system_clock(count=clock)
      endif

      !yanked from the gfortran documentation
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

    endif

#ifdef MPI
    !initialize random seeds for secondary processes
    allocate(allseeds(n, nprocs))

    if (mpirank .eq. 0) then                     !master process
       allocate(rand(nprocs*n))

       if (input_seed .gt. 0) then
         rand = input_seed + 13 * (/ (i, i = 1, nprocs*n) /)
       else
         call random_number(rand)
         call system_clock(count=clock)
         rand = 0.5_dp*clock*(1_dp + rand)
       endif

       allseeds = reshape(int(rand), (/n, nprocs/))
       deallocate(rand)
    end if

    call MPI_Scatter(allseeds, n, MPI_INTEGER, seed, n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    if (mpirank .ne. 0) call random_seed(put = seed)

    deallocate(seed, allseeds)
#endif


  end subroutine init_all_random_seeds


end module init
