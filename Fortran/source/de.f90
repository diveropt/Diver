module de

use detypes
use mutation
use crossover
use posterior

implicit none

private
public run_de

contains 


  !Main differential evolution routine.  
  subroutine run_de(func, prior, lowerbounds, upperbounds, numciv, numgen, NP, F, Cr, tol)
    real, external :: func, prior 				!function to be optimized, prior funtions
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    integer, intent(in) :: numciv 				!maximum number of civilisations
    integer, intent(in) :: numgen 				!maximum number of generations per civilisation
    integer, intent(in) :: NP 					!population size (individuals per generation)
    real, intent(in) :: F 					!scale factor
    real, intent(in) :: Cr 					!crossover factor
    real, intent(in) :: tol					!tolerance in log-evidence for 
      
    integer :: D 		!dimension of parameter space; we know this from the bounds given
    type(deparams) :: params 	!carries the differential evolution parameters 

    type(population) :: X, BF                  			!population of target vectors, best-fit vector          
    real, dimension(size(lowerbounds)) :: V, U			!donor, trial vectors

    integer :: fcall, accept			    !fcall counts function calls, accept counts acceptance rate
    integer :: civ, gen, n			    !civ, gen, n for iterating civilisation, generation, population loops
    
    real, dimension(size(lowerbounds)) :: avgvector, bestvector	!for calculating final average and best fit
    real :: bestvalue
    integer :: bestloc(1)

    logical :: calcZ = .false.					!whether to bother with posterior and evidence or not
    real :: Zold, Z = 0						!evidence
    integer :: convcount = 0					!number of times delta ln Z < tol in a row so far
    integer, parameter :: convcountreq = 4			!number of times delta ln Z < tol in a row for convergence
    
    if (tol .gt. 0.0) calcZ = .true.

    D=size(lowerbounds)
    params%NP = NP
    params%D = D
    params%F = F
    params%Cr = Cr

    allocate(BF%vectors(1,D), BF%values(1))

    fcall = 0
    BF%values(1) = huge(BF%values(1))
    write (*,*) 'Begin DE'
    write (*,*) 'Parameters:'
    write (*,*) ' NP=', params%NP
    write (*,*) ' F=', params%F 
    write (*,*) ' Cr=', params%Cr 

    !Run a number of sequential DE optmisations, exiting either after a set number of
    !runs through or after the evidence has been calculated to a desired accuracy
    do civ = 1, numciv

      write (*,*) '-----------------------------'
      write (*,*) 'Civilisation: ', civ

      !Initialise the first generation
      call initialize(X, params, lowerbounds, upperbounds, fcall, func)

      !Internal (normal) DE loop: calculates population for each generation
      do gen = 2, numgen 

         write (*,*) '  -----------------------------'
         write (*,*) '  Generation: ', gen
   
         accept = 0

         !$OMP PARALLEL DO
         do n=1, NP                            !evolves one member of the population

            V = rand1mutation(X, n, params)    !donor vectors
            U = bincrossover(X, V, n, params)  !trial vectors

            !choose next generation of target vectors
            call selection(X, U, n, lowerbounds, upperbounds, fcall, func, accept)
 
            write (*,*) n, X%vectors(n, :), '->', X%values(n)
         end do
         !$END OMP PARALLEL DO
 
         write (*,*) '  Acceptance rate: ', accept/real(NP)

         if (calcZ) then          
           !Find weights for posterior pdf / evidence calculation
           call getweights(X,prior)
           !FIXME multiplicities for outputting in chains = X%weights*exp(-X%values)
           !Update evidence
           Z = Z + sum(X%weights*exp(-X%values))
         endif

         !Check generation-level convergence: if satisfied, exit loop (!FIXME to be implemented)

      end do

      avgvector = sum(X%vectors, dim=1)/real(NP)
      bestloc = minloc(X%values)
      bestvector = X%vectors(bestloc(1),:)
      bestvalue = minval(X%values)

      !Update current best fit
      if (bestvalue .le. BF%values(1)) then
        BF%values(1) = bestvalue
        BF%vectors(1,:) = bestvector
      endif

      write (*,*)
      write (*,*) '  ============================='
      write (*,*) '  Number of generations in this civilisation: ', min(gen,numgen)
      write (*,*) '  Average final vector in this civilisation: ', avgvector
      write (*,*) '  Value at average final vector in this civilisation: ', func(avgvector, fcall) 
      write (*,*) '  Best final vector in this civilisation: ', bestvector
      write (*,*) '  Value at best final vector in this civilisation: ', bestvalue
      write (*,*) '  Cumulative function calls: ', fcall
      
      !Break out if posterior/evidence is converged
      if (calcZ .and. Z/Zold .lt. exp(tol)) then
        if (convcount .eq. convcountreq-1) exit
        convcount = convcount + 1
      else
        convcount = 0
      endif
      Zold = Z

    enddo

    write (*,*)
    write (*,*) '============================='
    write (*,*) 'Number of civilisations: ', min(civ,numciv)
    write (*,*) 'Best final vector: ', BF%vectors(1,:)
    write (*,*) 'Value at best final vector: ', BF%values(1)
    if (calcZ) write (*,*) 'ln(Evidence): ', log(Z)
    write (*,*) 'Total Function calls: ', fcall

    deallocate(X%vectors, X%values, BF%vectors, BF%values)

  end subroutine run_de




  !initializes first generation of target vectors
  subroutine initialize(X, params, lowerbounds, upperbounds, fcall, func) 

    type(population), intent(out) :: X
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    real, external :: func
    integer :: i

    write (*,*) '-----------------------------'
    write (*,*) 'Generation: ', '1'

    allocate(X%vectors(params%NP, params%D), X%values(params%NP), X%weights(params%NP)) !deallocated at end of run_de

    !$OMP PARALLEL DO
    do i=1,params%NP
       call random_number(X%vectors(i,:))
       X%vectors(i,:) = X%vectors(i,:)*(upperbounds - lowerbounds) + lowerbounds

       X%values(i) = func(X%vectors(i,:), fcall)
       write (*,*) i, X%vectors(i, :), '->', X%values(i)
    end do      
    !$END OMP PARALLEL DO
  end subroutine initialize



  !select next generation of target vectors
  subroutine selection(X,U, n, lowerbounds, upperbounds, fcall, func, accept) 

    type(population), intent(inout) :: X
    real, dimension(:), intent(in) :: U
    integer, intent(inout) :: fcall, accept
    integer, intent(in) :: n 
    real, dimension(:), intent(in) :: lowerbounds, upperbounds
    real, external :: func
    real :: trialvalue

    !check that results stay within bounds. 'Brick wall'
    if (all(U .le. upperbounds) .and. all(U .ge. lowerbounds)) then

       trialvalue = func(U, fcall)

       !when the trial vector is at least as good, use it for the next generation
       if (trialvalue .le. X%values(n)) then
          X%vectors(n,:) = U 
          X%values(n) = trialvalue
          accept = accept + 1
       end if
    end if
  end subroutine selection


end module de
