module init

use detypes
use converge

implicit none

private
public param_assign, initialize, init_random_seed

contains 

  !Assign parameter values (defaults if not specified) to run_params and print DE parameter values to screen

  subroutine param_assign(run_params, lowerbounds, upperbounds, nDerived, maxciv, maxgen, NP, F, Cr, lambda, current, &
                          expon, bndry, jDE, removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount)

    type(codeparams), intent(out) :: run_params 
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    integer, intent(in), optional :: nDerived	 		!input number of derived quantities to output
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

    character (len=30) :: DEstrategy, Fsize			!for printing mutation/crossover DE strategy

    run_params%D=size(lowerbounds)
    if (size(upperbounds) .ne. run_params%D) then
       write (*,*) 'ERROR: parameter space bounds do not have the same dimensionality'
       stop
    endif

    if (any(lowerbounds .ge. upperbounds)) then
       write (*,*) 'ERROR: invalid parameter space bounds.'
       stop
    endif

    if (present(nDerived)) then 
      call setIfPositive_int(nDerived, run_params%D_derived, 'nDerived')
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

       write (Fsize, *) run_params%DE%Fsize			!number of mutation scale factors
       Fsize = adjustl(Fsize)
       DEstrategy = trim(DEstrategy)//trim(Fsize) 	
       
       if(run_params%DE%expon) then                  		!crossover strategy
          DEstrategy = trim(DEstrategy)//'/exp'
       else
          DEstrategy = trim(DEstrategy)//'/bin'
       end if
    endif jDEset

    !print feedback about strategy choice, values of parameters to the screen
    write (*,*) DEstrategy
    write (*,*) 'Parameters:'
    write (*,*) ' NP =', run_params%DE%NP
    if ((run_params%DE%lambda .ne. 1.) .and. (run_params%DE%lambda .ne. 0.)) then
       write (*,*) ' lambda =', run_params%DE%lambda
    endif
    if (.not. run_params%DE%jDE) then                   
       write (*,*) ' F =', run_params%DE%F  
       write (*,*) ' Cr =', run_params%DE%Cr 
    endif

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

  end subroutine param_assign


  !set parameter only if value is positive (real parameter)
  subroutine setIfPositive_real(invar, outvar, string)

    real :: invar, outvar
    character(LEN=*) :: string

    if (invar .le. 0.0) then
      write(*,*) 'ERROR: '//string//' cannot be negative.'
      stop
    else
      outvar = invar
    endif

  end subroutine setIfPositive_real


  !set parameter only if value is positive (integer parameter)
  subroutine setIfPositive_int(invar, outvar, string)

    integer :: invar, outvar
    character(LEN=*) :: string

    if (invar .le. 0) then
      write(*,*) 'ERROR: '//string//' cannot be negative.'
      stop
    else
      outvar = invar
    endif

  end subroutine setIfPositive_int


  !initializes first generation of target vectors
  subroutine initialize(X, run_params, lowerbounds, upperbounds, fcall, func, quit) 

    type(population), intent(inout) :: X
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    logical, intent(inout) :: quit
    real, external :: func
    integer :: i
    real, dimension(run_params%DE%NP) :: rand

    if (verbose) write (*,*) '-----------------------------'
    if (verbose) write (*,*) 'Generation: ', '1'

    X%multiplicities = 1.d0 !Initialise to 1 in case posteriors are not calculated

    !$OMP PARALLEL DO
    do i=1,run_params%DE%NP
       call random_number(X%vectors(i,:))
       X%vectors(i,:) = X%vectors(i,:)*(upperbounds - lowerbounds) + lowerbounds
       X%values(i) = func(X%vectors(i,:), X%derived(i,:), fcall, quit)
       if (run_params%DE%jDE) then !initialize population of F and Cr parameters
          call random_number(rand)
          X%FjDE = rand*0.9 + 0.1
          call random_number(rand)
          X%CrjDE = rand
          if (verbose) write (*,*) i, X%vectors(i, :), '->', X%values(i), '|', X%FjDE(i), X%CrjDE(i)
       else
          if (verbose) write (*,*) i, X%vectors(i, :), '->', X%values(i)
       end if
    end do      
    !$END OMP PARALLEL DO

    if (converged(X, 1)) write (*,*) 'ERROR: initial population converges.' !initializing converged()

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


end module init
