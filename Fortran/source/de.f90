module de

use detypes
use mutation
use crossover
use converge
use posterior

implicit none

private
public run_de

!logical, parameter :: verbose = .true.				!moved to detypes

contains 


  !Main differential evolution routine.  
  subroutine run_de(func, prior, lowerbounds, upperbounds, maxciv, maxgen, NP, F, Cr, lambda, current, exp, bndry, tolerance)
    real, external :: func, prior 				!function to be optimized, prior functions
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    integer, intent(in), optional :: maxciv 			!maximum number of civilisations
    integer, intent(in), optional :: maxgen 			!maximum number of generations per civilisation
    integer, intent(in), optional :: NP 			!population size (individuals per generation)
    real, dimension(:), intent(in), optional :: F		!scale factor(s).  Note that this must be entered as an array.
    real, intent(in), optional :: Cr 				!crossover factor
    real, intent(in), optional :: lambda 			!mixing factor between best and rand/current
    logical, intent(in), optional :: current 			!use current vector for mutation
    logical, intent(in), optional :: exp 			!use exponential crossover
    integer, intent(in), optional :: bndry                      !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    real, intent(in), optional :: tolerance			!tolerance in log-evidence for 
     
    integer :: D 						!dimension of parameter space; we know this from the bounds given
    type(deparams) :: params 					!carries the differential evolution parameters 
    integer :: bconstrain					!boundary constraint parameter

    type(population), target :: X, BF           		!population of target vectors, best-fit vector          
    real, dimension(size(lowerbounds)) :: V, U			!donor, trial vectors

    integer :: fcall, accept					!fcall counts function calls, accept counts acceptance rate
    integer :: civ, gen, n					!civ, gen, n for iterating civilisation, generation, population loops
    integer :: numciv, numgen					!maximum number of civilizations, generations
    real :: tol							!tolerance in log-evidence

    real, dimension(size(lowerbounds)) :: avgvector, bestvector	!for calculating final average and best fit
    real :: bestvalue
    integer :: bestloc(1)

    logical :: calcZ = .false.					!whether to bother with posterior and evidence or not
    real :: Zold, Z = 0						!evidence
    integer :: convcount = 0					!number of times delta ln Z < tol in a row so far
    integer, parameter :: convcountreq = 4			!number of times delta ln Z < tol in a row for convergence

    logical testing
    

    write (*,*) '============================='
    write (*,*) ' ******** Begin DE *********'

    if (any(lowerbounds .ge. upperbounds)) then
       write (*,*) 'ERROR: invalid parameter space bounds.'
    else !proceed with program

       D=size(lowerbounds)

       !assign specified or default values to params, bconstrain
       call param_assign(params, bconstrain, D, NP=NP, F=F, Cr=Cr, lambda=lambda, current=current, exp=exp, bndry=bndry)

       if (present(maxciv)) then
          numciv = maxciv
       else
          numciv = 500						!arbitrary default value for numciv
       end if
       if (present(maxgen)) then
          numgen = maxgen
       else 
          numgen = 200						!arbitrary default value for numgen
       end if

       if (present(tolerance)) then 
          tol = tolerance
       else
          tol = 1e-3
       end if

       if (tol .gt. 0.0) calcZ = .true.

       allocate(BF%vectors(1,D), BF%values(1))
    
       fcall = 0
       BF%values(1) = huge(BF%values(1))

       !If required, initialise the linked tree used for estimating the evidence and posterior
       if (calcZ) call initree(lowerbounds,upperbounds)

       !Run a number of sequential DE optimisations, exiting either after a set number of
       !runs through or after the evidence has been calculated to a desired accuracy
       do civ = 1, numciv

          if (verbose) write (*,*) '-----------------------------'
          if (verbose) write (*,*) 'Civilisation: ', civ

          !Initialise the first generation
          call initialize(X, params, lowerbounds, upperbounds, fcall, func)
          if (calcZ) call doBayesian(X, Z, prior, fcall)        

          !Internal (normal) DE loop: calculates population for each generation
          do gen = 2, numgen 

             if (verbose) write (*,*) '  -----------------------------'
             if (verbose) write (*,*) '  Generation: ', gen
   
             accept = 0

             !$OMP PARALLEL DO
             do n=1, params%NP                    	!evolves one member of the population

                V = genmutation(X, n, params)    	!donor vectors
                U = gencrossover(X, V, n, params)  	!trial vectors

                !choose next generation of target vectors
                call selection(X, U, n, lowerbounds, upperbounds, bconstrain, fcall, func, accept)
 
                if (verbose) write (*,*) n, X%vectors(n, :), '->', X%values(n)
             end do
             !$END OMP PARALLEL DO
 
             if (verbose) write (*,*) '  Acceptance rate: ', accept/real(params%NP)

             if (calcZ) call doBayesian(X, Z, prior, fcall)        

             if (converged(X, gen)) exit             !Check generation-level convergence: if satisfied, exit loop

          end do

          avgvector = sum(X%vectors, dim=1)/real(params%NP)
          bestloc = minloc(X%values)
          bestvector = X%vectors(bestloc(1),:)
          bestvalue = minval(X%values)
          
          !Update current best fit
          if (bestvalue .le. BF%values(1)) then
             BF%values(1) = bestvalue
             BF%vectors(1,:) = bestvector
          endif
          
          if (verbose) write (*,*)
          if (verbose) write (*,*) '  ============================='
          if (verbose) write (*,*) '  Number of generations in this civilisation: ', min(gen,numgen)
          if (verbose) write (*,*) '  Average final vector in this civilisation: ', avgvector
          if (verbose) write (*,*) '  Value at average final vector in this civilisation: ', func(avgvector, fcall) 
          if (verbose) write (*,*) '  Best final vector in this civilisation: ', bestvector
          if (verbose) write (*,*) '  Value at best final vector in this civilisation: ', bestvalue
          if (verbose) write (*,*) '  Cumulative function calls: ', fcall
      
          !if (calcZ) write(*,*) abs(log(Z)-log(Zold)), tol
          !Break out if posterior/evidence is converged
          if (calcZ .and. abs(log(Z)-log(Zold)) .lt. tol) then
             if (convcount .eq. convcountreq-1) exit
             convcount = convcount + 1
          else
             convcount = 0
          endif
          Zold = Z
          
       enddo

!    if (verbose) write (*,*)
!    if (verbose) write (*,*) '============================='
!    if (verbose) write (*,*) 'Number of civilisations: ', min(civ,numciv)
!    if (verbose) write (*,*) 'Best final vector: ', BF%vectors(1,:)
!    if (verbose) write (*,*) 'Value at best final vector: ', BF%values(1)
!    if (calcZ) write (*,*)   'ln(Evidence): ', log(Z)
!    if (verbose) write (*,*) 'Total Function calls: ', fcall

       write (*,*) '============================='
       write (*,*) 'Number of civilisations: ', min(civ,numciv)
       write (*,*) 'Best final vector: ', BF%vectors(1,:)
       write (*,*) 'Value at best final vector: ', BF%values(1)
       if (calcZ) write (*,*)   'ln(Evidence): ', log(Z)
       write (*,*) 'Total Function calls: ', fcall

       deallocate(X%vectors, X%values, params%F, BF%vectors, BF%values)

    end if

  end subroutine run_de



  !assign parameter values (defaults if not specified) and print values to screen
  !moved this to its own subroutine for ease-of-reading in the main program
  subroutine param_assign(params, bconstrain, D, NP, F, Cr, lambda, current, exp, bndry)

    type(deparams), intent(out) :: params 
    integer, intent(out) :: bconstrain		!boundary constraints for selection
    integer, intent(in) :: D 			  
    integer, optional, intent(in) :: NP 
    real, dimension(:), optional, intent(in) :: F 	
    real, optional, intent(in) :: Cr 	
    real, optional, intent(in) :: lambda 
    logical, optional, intent(in) :: current 
    logical, optional, intent(in) :: exp 
    integer, optional, intent(in) :: bndry

    character (len=22) :: DEstrategy, Fsize	!for printing mutation/crossover DE strategy

    !assigning specified or default values of the DE parameters to params

    params%D = D

    if (present(F)) then
       if (any(F .le. 0)) write (*,*) 'WARNING: some elements of F are 0 or negative. DE may not converge properly.'
       if (any(F .ge. 1)) write (*,*) 'WARNING: some elements of F are 1 or greater. DE may not converge properly.'
       allocate(params%F(size(F)))
       params%F = F
    else
       allocate(params%F(1))
       params%F = (/0.7/) 			!rule of thumb: 0.4<F<1.0
    end if

    if (present(NP)) then
       if (NP .ge. (2*size(params%F) + 2)) then !required for picking unique vectors during mutation
          params%NP = NP
       else !nb: if current=true, NP=2*size(params%F)+1 would be ok, but not implemented
          write (*,*) 'WARNING: NP specified is too small. Using smallest permitted NP.'
          params%NP = 2*size(params%F) + 2
       end if
    else
       params%NP = maxval( (/10*D, 2*size(params%F) + 2/) )	!conservative rule-of-thumb choice 
    end if

    if (present(Cr)) then  
       if (Cr .lt. 0.0) then
          write (*,*) 'WARNING: Cr < 0. Using Cr = 0.' !although Cr<0 is functionally equivalent to Cr=0
          params%Cr = 0.0
       else if (Cr .gt. 1.0) then
          write (*,*) 'WARNING: Cr > 1. Using Cr = 1.' !although Cr>1 is functionally equivalent to Cr=1
          params%Cr = 1.0
       else
          params%Cr = Cr
       end if
    else
       params%Cr = 0.9 
    end if

    if (present(lambda)) then
       if (lambda .lt. 0.0) write (*,*) 'WARNING: lambda < 0. DE may not converge properly.'
       if (lambda .gt. 1.0) write (*,*) 'WARNING: lambda > 1. DE may not converge properly.'
       params%lambda = lambda
    else
       params%lambda = 0.0     			!default rand/1/bin
    end if

    if (present(current)) then 
       params%current = current
    else
       params%current = .false. 		!default rand/1/bin
    end if

    if (present(exp)) then
       params%exp = exp
    else
       params%exp = .false.     		!default rand/1/bin
    end if

    !printing the parameter choice and DE mutation/crossover strategy to screen

    if (params%lambda .eq. 0) then  		!mutation strategy
       if (params%current) then
          DEstrategy = 'current/'
       else
          DEstrategy = 'rand/'
       end if
    else if (params%lambda .eq. 1) then
       DEstrategy = 'best/'
    else 
       if (params%current) then
          DEstrategy = 'current-to-best/'
       else
          DEstrategy = 'rand-to-best/'
       end if
    end if

    write (Fsize, *) size(params%F)		!number of mutation scale factors
    Fsize = adjustl(Fsize)
    DEstrategy = trim(DEstrategy)//trim(Fsize) 	

    if(params%exp) then                  	!crossover strategy
       DEstrategy = trim(DEstrategy)//'/exp'
    else
       DEstrategy = trim(DEstrategy)//'/bin'
    end if

    write (*,*) DEstrategy
    write (*,*) 'Parameters:'
    write (*,*) ' NP =', params%NP
    if ((params%lambda .ne. 1.) .and. (params%lambda .ne. 0.)) write (*,*) ' lambda =', params%lambda
    write (*,*) ' F =', params%F  
    write (*,*) ' Cr =', params%Cr 

    !assigning specified/default boundary constraints and printing to screen

    if (present(bndry)) then
       bconstrain = bndry
    else
       bconstrain = 1 				!default brick wall boundary constraints
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

  end subroutine param_assign




  !initializes first generation of target vectors
  subroutine initialize(X, params, lowerbounds, upperbounds, fcall, func) 

    type(population), intent(out) :: X
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    real, external :: func
    integer :: i

    if (verbose) write (*,*) '-----------------------------'
    if (verbose) write (*,*) 'Generation: ', '1'

    allocate(X%vectors(params%NP, params%D), X%values(params%NP), X%weights(params%NP)) !deallocated at end of run_de

    !$OMP PARALLEL DO
    do i=1,params%NP
       call random_number(X%vectors(i,:))
       X%vectors(i,:) = X%vectors(i,:)*(upperbounds - lowerbounds) + lowerbounds

       X%values(i) = func(X%vectors(i,:), fcall)
       if (verbose) write (*,*) i, X%vectors(i, :), '->', X%values(i)
    end do      
    !$END OMP PARALLEL DO

    if (converged(X, 1)) write (*,*) 'ERROR: initial population converges.'

  end subroutine initialize



  !select next generation of target vectors
  subroutine selection(X, U, n, lowerbounds, upperbounds, bconstrain, fcall, func, accept) 

    type(population), intent(inout) :: X
    real, dimension(:), intent(in) :: U
    integer, intent(inout) :: fcall, accept
    integer, intent(in) :: n 
    real, dimension(:), intent(in) :: lowerbounds, upperbounds
    integer, intent(in) :: bconstrain      
    real, external :: func
    real :: trialvalue
    real, dimension(size(U)) :: trialvector  

    if (any(U(:) .gt. upperbounds) .or. any(U(:) .lt. lowerbounds)) then 
       !trial vector exceeds parameter space bounds: apply boundary constraints
       select case (bconstrain)
          case (1)                           !'brick wall'
             trialvalue = huge(1.0)
             trialvector(:) = X%vectors(n,:)
          case (2)                           !randomly re-initialize
             call random_number(trialvector(:))
             trialvector(:) = trialvector(:)*(upperbounds - lowerbounds) + lowerbounds
             trialvalue = func(trialvector, fcall)
          case (3)                           !reflection
             trialvector = U
             where (U .gt. upperbounds)  trialvector = upperbounds - (U - upperbounds)
             where (U .lt. lowerbounds)  trialvector = lowerbounds + (lowerbounds - U)
             trialvalue = func(trialvector, fcall)
          case default                       !boundary constraints not enforced
             trialvector = U                
             trialvalue = func(U(:), fcall)  
          end select

    else                                     !trial vector is within parameter space bounds, so use it
       trialvector = U                    
       trialvalue = func(U(:), fcall)  
    end if

    !when the trial vector is at least as good, use it for the next generation
    if (trialvalue .le. X%values(n)) then
       X%vectors(n,:) = trialvector 
       X%values(n) = trialvalue
       accept = accept + 1
    end if

  end subroutine selection



  !Get posterior weights and update evidence
  subroutine doBayesian(X, Z, prior, fcall)
  
    type(population), intent(inout) :: X		!current generation
    real, intent(inout) :: Z				!evidence
    real, external :: prior 				!prior funtion
    integer, intent(in) :: fcall			!running number of samples
    integer, save :: fcall_prev = 0			!last number of samples
    
    !Find weights for posterior pdf / evidence calculation
    call getweights(X,prior)
    
    !FIXME multiplicities for outputting in chains = X%weights/fcall*exp(-X%values)

    !Update evidence
    Z = (Z*dble(fcall_prev) + sum(X%weights*exp(-X%values)))/dble(fcall)

    !Save number of points for next time
    fcall_prev = fcall

    !write(*,*) Z, sum(X%weights/fcall*exp(-X%values)), sum(X%weights), sum(exp(-X%values))

  end subroutine doBayesian



end module de
