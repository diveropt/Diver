module init

use detypes
use converge

implicit none

private
public param_assign, initialize

contains 

  !Assign parameter values (defaults if not specified) and print values to screen

  subroutine param_assign(run_params, bconstrain, lowerbounds, upperbounds, nDerived, maxciv, maxgen, NP, F, &
                          Cr, lambda, current, expon, bndry, tolerance, tolcount, savecount)

    type(codeparams), intent(out) :: run_params 
    integer, intent(out) :: bconstrain				!boundary constraints for selection
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
    real, intent(in), optional :: tolerance			!input tolerance in log-evidence
    integer, intent(in), optional :: tolcount	 		!input number of times delta ln Z < tol in a row for convergence
    integer, intent(in), optional :: savecount			!save progress every savecount generations

    character (len=22) :: DEstrategy, Fsize			!for printing mutation/crossover DE strategy

    run_params%D=size(lowerbounds)
    if (size(upperbounds) .ne. run_params%D) then
       write (*,*) 'ERROR: parameter space bounds do not have the same dimensionality'
       stop
    endif

    if (present(nDerived)) then 
       run_params%D_derived = nDerived
    else
       run_params%D_derived = 0					!default is no derived quantities
    end if

    if (present(maxciv)) then
       run_params%numciv = maxciv
    else
       run_params%numciv = 2000					!arbitrary default value for numciv
    end if

    if (present(maxgen)) then
       run_params%numgen = maxgen
    else 
       run_params%numgen = 100					!arbitrary default value for numgen
    end if

    if (present(NP)) then
       if (NP .ge. (2*size(run_params%DE%F) + 2)) then 		!required for picking unique vectors during mutation
          run_params%DE%NP = NP
       else !nb: if current=true, NP=2*size(run_params%DE%F)+1 would be ok, but not implemented
          write (*,*) 'WARNING: NP specified is too small. Using smallest permitted NP.'
          run_params%DE%NP = 2*size(run_params%DE%F) + 2
       end if
    else
       run_params%DE%NP = maxval( (/10*run_params%D, 2*size(run_params%DE%F) + 2/) )	!conservative rule-of-thumb choice 
    end if

    if (present(F)) then
       if (any(F .le. 0)) write (*,*) 'WARNING: some elements of F are 0 or negative. DE may not converge properly.'
       if (any(F .ge. 1)) write (*,*) 'WARNING: some elements of F are 1 or greater. DE may not converge properly.'
       allocate(run_params%DE%F(size(F)))
       run_params%DE%F = F
    else
       allocate(run_params%DE%F(1))
       run_params%DE%F = (/0.7/) 				!rule of thumb: 0.4<F<1.0
    end if

    if (present(Cr)) then  
       if (Cr .lt. 0.0) then
          write (*,*) 'WARNING: Cr < 0. Using Cr = 0.' 		!although Cr<0 is functionally equivalent to Cr=0
          run_params%DE%Cr = 0.0
       else if (Cr .gt. 1.0) then
          write (*,*) 'WARNING: Cr > 1. Using Cr = 1.' 		!although Cr>1 is functionally equivalent to Cr=1
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
       run_params%DE%lambda = 0.0     				!default rand/1/bin
    end if

    if (present(current)) then 
       run_params%DE%current = current
    else
       run_params%DE%current = .false. 				!default rand/1/bin
    end if

    if (present(expon)) then
       run_params%DE%expon = expon
    else
       run_params%DE%expon = .false.     			!default rand/1/bin
    end if

    !printing the parameter choice and DE mutation/crossover strategy to screen

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

    write (Fsize, *) size(run_params%DE%F)			!number of mutation scale factors
    Fsize = adjustl(Fsize)
    DEstrategy = trim(DEstrategy)//trim(Fsize) 	

    if(run_params%DE%expon) then                  		!crossover strategy
       DEstrategy = trim(DEstrategy)//'/exp'
    else
       DEstrategy = trim(DEstrategy)//'/bin'
    end if

    write (*,*) DEstrategy
    write (*,*) 'Parameters:'
    write (*,*) ' NP =', run_params%DE%NP
    if ((run_params%DE%lambda .ne. 1.) .and. (run_params%DE%lambda .ne. 0.)) write (*,*) ' lambda =', run_params%DE%lambda
    write (*,*) ' F =', run_params%DE%F  
    write (*,*) ' Cr =', run_params%DE%Cr 

    !assigning specified/default boundary constraints and printing to screen

    if (present(bndry)) then
       bconstrain = bndry
    else
       bconstrain = 1 						!default brick wall boundary constraints
    end if

    select case (bconstrain)
       case (1) 
          write (*,*) 'Brick wall boundary constraints'
       case (2)
          write (*,*) 'Random re-initialization boundary constraints'
       case (3)
          write (*,*) 'Reflective boundary constraints'
       case default
          write (*,*) 'WARNING: Invalid value entered for bndry. Boundary constraints not enforced.'
    end select

    if (present(tolerance)) then 
       run_params%tol = tolerance
    else
       run_params%tol = 1.d-4					!default for tolerance
    end if
    if (run_params%tol .gt. 0.0) run_params%calcZ = .true.	!enable posterior and evidence calculation if tol is +ve

    if (present(tolcount)) then 
       run_params%convcountreq = tolcount
    else
       run_params%convcountreq = 4				!default for tolerance counter
    end if

    if (present(savecount)) then 
       run_params%savefreq = savecount
    else
       run_params%savefreq = 1					!default for tolerance counter
    end if

  end subroutine param_assign




  !initializes first generation of target vectors
  subroutine initialize(X, run_params, lowerbounds, upperbounds, fcall, func) 

    type(population), intent(out) :: X
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    real, external :: func
    integer :: i

    if (verbose) write (*,*) '-----------------------------'
    if (verbose) write (*,*) 'Generation: ', '1'

    !Deallocated at end of run_de
    allocate(X%vectors(run_params%DE%NP, run_params%D))
    allocate(X%derived(run_params%DE%NP, run_params%D_derived))
    allocate(X%values(run_params%DE%NP), X%weights(run_params%DE%NP), X%multiplicities(run_params%DE%NP)) 
    X%multiplicities = 1.d0 !Initialise to 1 in case posteriors are not calculated

    !$OMP PARALLEL DO
    do i=1,run_params%DE%NP
       call random_number(X%vectors(i,:))
       X%vectors(i,:) = X%vectors(i,:)*(upperbounds - lowerbounds) + lowerbounds
       X%values(i) = func(X%vectors(i,:), X%derived(i,:), fcall)
       if (verbose) write (*,*) i, X%vectors(i, :), '->', X%values(i)
    end do      
    !$END OMP PARALLEL DO

    if (converged(X, 1)) write (*,*) 'ERROR: initial population converges.'

  end subroutine initialize



end module init
