module init

use detypes
use mutation, only: init_FjDE    !initializes scale factors for jDE
use crossover, only: init_CrjDE  !initializes crossovers for jDE
use selection, only: roundvector, replace_generation

#ifdef USEMPI
  use MPI
#endif

implicit none

private
public param_assign, initialize, init_random_seed, quit_de, int_to_string

contains 

  !Assign parameter values (defaults if not specified) to run_params and print DE parameter values to screen

  subroutine param_assign(run_params, lowerbounds, upperbounds, nDerived, discrete, maxciv, maxgen, NP, F, Cr, lambda, &
                          current, expon, bndry, jDE, removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount)

    type(codeparams), intent(out) :: run_params 
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    integer, intent(in), optional :: nDerived	 		!input number of derived quantities to output
    integer, dimension(:), intent(in), optional :: discrete     !lists all discrete dimensions of parameter space
    integer, intent(in), optional :: maxciv 			!maximum number of civilisations
    integer, intent(in), optional :: maxgen 			!maximum number of generations per civilisation
    integer, intent(in), optional :: NP 			!population size (individuals per generation)
    real, dimension(:), intent(in), optional :: F		!scale factor(s).  Note that this must be entered as an array.
    real, intent(in), optional :: Cr 				!crossover factor
    real, intent(in), optional :: lambda 			!mixing factor between best and rand/current
    logical, intent(in), optional :: current 			!use current vector for mutation
    logical, intent(in), optional :: expon 			!use exponential crossover
    integer, intent(in), optional :: bndry                      !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional :: jDE                        !use self-adaptive DE 
    logical, intent(in), optional :: removeDuplicates           !weed out duplicate vectors within a single generation
    logical, intent(in), optional  :: doBayesian                !calculate log evidence and posterior weightings
    real, intent(in), optional  :: maxNodePop                   !population at which node is partitioned in binary space partitioning for posterior
    real, intent(in), optional :: Ztolerance			!input tolerance in log-evidence
    integer, intent(in), optional :: savecount			!save progress every savecount generations

    integer :: mpiprocs, mpirank, ierror                        !number of processes running, rank of current process, error code
    character (len=30) :: DEstrategy                    	!for printing mutation/crossover DE strategy

#ifdef USEMPI
    call MPI_Comm_size(MPI_COMM_WORLD, mpiprocs, ierror)  !gives the total num of processes. If no MPI, set to 1.
    call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierror)   !rank of the current process. If no MPI, set to 0.
#else
    mpiprocs = 1
    mpirank = 0
#endif

    if (mpirank .eq. 0) then
       write (*,*) '============================='
       write (*,*) ' ******** Begin DE *********'
    end if

    run_params%mpirank = mpirank
    run_params%D=size(lowerbounds)

    if (size(upperbounds) .ne. run_params%D) call quit_de('ERROR: parameter space bounds do not have the same dimensionality')
    if (any(lowerbounds .ge. upperbounds)) call quit_de('ERROR: invalid parameter space bounds.')

    if (present(nDerived)) then 
      call setIfPositive_int(nDerived, run_params%D_derived, 'nDerived') !FIXME: this prevents setting nDerived=0
    else
       run_params%D_derived = 0					!default is no derived quantities
    end if

    if (present(discrete)) then
       if (any(discrete .gt. run_params%D) .or. any(discrete .lt. 1)) then     
          call quit_de('ERROR: Discrete dimensions specified must not be < 1 or > '//trim(int_to_string(run_params%D)))
       end if
       !also check that discrete dimensions are not doubly-specified? Doesn't seem to crash...
       run_params%D_discrete = size(discrete)
       allocate(run_params%discrete(run_params%D_discrete))
       run_params%discrete = discrete
    else
       run_params%D_discrete = 0
       allocate(run_params%discrete(0))
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
       run_params%maxNodePop = 1.9				!default for maxNodePop
    end if

    if (present(Ztolerance)) then 
       call setIfPositive_real(Ztolerance, run_params%tol, 'Ztolerance')
    else
       run_params%tol = 0.01 					!default for tolerance
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
          if (any(F .le. 0.)) write (*,*) 'WARNING: some elements of F are 0 or negative. DE may not converge properly.'
          if (any(F .ge. 1.)) write (*,*) 'WARNING: some elements of F are 1 or greater. DE may not converge properly.'
          run_params%DE%Fsize = size(F)
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = F
       else
          run_params%DE%Fsize = 1
          allocate(run_params%DE%F(run_params%DE%Fsize))
          run_params%DE%F = (/0.7/) 				!rule of thumb: 0.4<F<1.0
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
          if (Cr .lt. 0.0) then
             write (*,*) 'WARNING: Cr < 0. Using Cr = 0.' 	!although Cr<0 is functionally equivalent to Cr=0
             run_params%DE%Cr = 0.0
          else if (Cr .gt. 1.0) then
             write (*,*) 'WARNING: Cr > 1. Using Cr = 1.' 	!although Cr>1 is functionally equivalent to Cr=1
             run_params%DE%Cr = 1.0
          else
             run_params%DE%Cr = Cr
          end if
       else
          run_params%DE%Cr = 0.9 
       end if

       if (present(lambda)) then
          if (lambda .lt. 0.0) write (*,*) 'WARNING: lambda < 0. DE may not converge properly.'
          if (lambda .gt. 1.0) write (*,*) 'WARNING: lambda > 1. DE may not converge properly.'
          run_params%DE%lambda = lambda
       else
          run_params%DE%lambda = 0.0     			!default rand/1/bin
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
          run_params%DE%removeDuplicates = .false.              !with current mutation, duplicate are rare (CHECK THIS)
       else
          run_params%DE%removeDuplicates = .true.
       end if

       !for printing the parameter choice, DE mutation/crossover strategy, and boundary constraints to screen
       if (run_params%DE%lambda .eq. 0) then  			!mutation strategy
          if (run_params%DE%current) then
             DEstrategy = 'current/'
          else
             DEstrategy = 'rand/'
          end if
       else if (run_params%DE%lambda .eq. 1) then
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

    !split up work over multiple processes for MPI (if no mpi, single process does all the work)
    run_params%mpipopchunk = run_params%DE%NP/mpiprocs

    !print feedback about strategy choice, values of parameters to the screen
    if (run_params%mpirank .eq. 0) then
       write (*,*) DEstrategy
       if (size(run_params%discrete) .gt. 0) write (*,*) 'Discrete dimensions:', run_params%discrete
       write (*,*) 'Parameters:'
       write (*,*) ' NP = ', trim(int_to_string(run_params%DE%NP))
       if ((run_params%DE%lambda .ne. 1.) .and. (run_params%DE%lambda .ne. 0.)) then
          write (*,*) ' lambda =', run_params%DE%lambda
       endif
       if (.not. run_params%DE%jDE) then                   
          write (*,*) ' F =', run_params%DE%F
          write (*,*) ' Cr =', run_params%DE%Cr 
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

    real :: invar, outvar
    character(LEN=*) :: string

    if (invar .le. 0.0) then
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
  subroutine initialize(X, Xnew, run_params, lowerbounds, upperbounds, fcall, func, quit) 

    type(population), intent(inout) :: X
    type(population), intent(inout) :: Xnew
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    logical, intent(inout) :: quit
    real, external :: func
    real, dimension(run_params%D) :: evalvector  !vectors used to evaluate function (with 'discrete' dimensions rounded)
    integer :: n, m, accept

    X%multiplicities = 1.d0 !Initialise to 1 in case posteriors are not calculated

    if (run_params%DE%jDE) then !initialize population of F and Cr parameters
       Xnew%FjDE = init_FjDE(run_params)
       Xnew%CrjDE = init_CrjDE(run_params)
    end if

       !$OMP PARALLEL DO
       do m=1,run_params%mpipopchunk !loop over the vectors belonging to each population chunk

          n = run_params%mpipopchunk*run_params%mpirank + m !true population index (equal to m if no mpi)

          call random_number(Xnew%vectors(m,:))

          Xnew%vectors(m,:) = Xnew%vectors(m,:)*(upperbounds - lowerbounds) + lowerbounds
          evalvector = roundvector(Xnew%vectors(m,:), run_params)

          Xnew%values(m) = func(evalvector, Xnew%derived(m,:), fcall, quit)

          if (verbose .and. run_params%DE%jDE) write (*,*) n, evalvector, '->', Xnew%values(m), '|', Xnew%FjDE(m), Xnew%CrjDE(m)
          if (verbose .and. .not. run_params%DE%jDE) write (*,*) n, evalvector, '->', Xnew%values(m)

       end do
       !$END OMP PARALLEL DO

       call replace_generation(X, Xnew, run_params, accept, init=.true.)
    
  end subroutine initialize


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

  function int_to_string(int)
    integer, intent(in) :: int
    character (len=30) :: string, int_to_string

    write (string, *) int
    string = adjustl(string)
    string = trim(string)

    int_to_string = string

  end function int_to_string


  subroutine quit_de(error_message)
    character(LEN=*), intent(in), optional :: error_message
    integer ierror

    if (present(error_message)) write (*,*) error_message
#ifdef USEMPI
    call MPI_Finalize(ierror)
#endif
    stop

  end subroutine quit_de


end module init
