module de

use detypes
use mutation
use crossover

implicit none

private
public run_de

contains 


  !this is the main differential evolution routine.  
  subroutine run_de(func, lowerbounds, upperbounds, numgen, NP, F, Cr)
    real, external :: func !function to be optimized
    !boundaries of parameter space
    real, dimension(:), intent(in) :: lowerbounds, upperbounds 
    integer, intent(in) :: numgen !maximum number of generations
    integer, intent(in) :: NP !population size
    real, intent(in) :: F !scale factor
    real, intent(in) :: Cr !crossover factor
    
    integer :: D !dimension of parameter space; we know this from the bounds given
    type(deparams) :: params !carries the differential evolution parameters 

    type(population) :: X                              !population of target vectors          
    real, allocatable, dimension(:) :: V, U            !donor, trial vectors

    !fcall counts function calls, accept counts acceptance rate
    !gen, n for iterating generation loops, population loops
    integer :: gen, n, fcall, accept
    
    real, allocatable, dimension(:) :: avgvector       !for calculating final average

    D=size(lowerbounds)
    params%NP = NP
    params%D = D
    params%F = F
    params%Cr = Cr

    allocate(V(D), U(D))

    call initialize(X, params, lowerbounds, upperbounds, fcall, func)

    !main loop of program: calculates population for each generation
    do gen = 2, numgen 

       write (*,*) '-----------------------------'
       write (*,*) 'Generation: ', gen
   
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
 
       write (*,*) 'Acceptance rate: ', accept/real(NP)

       !check convergence criteria: if satisfied, exit loop (& increment generation counter?)

    end do

    deallocate(V, U)

    write (*,*) '============================='
    write (*,*) 'Number of generations: ', gen-1 !this may not always work...

    allocate(avgvector(D))
    avgvector = sum(X%vectors, dim=1)/real(NP)
    write (*,*) 'Average final vector: ', avgvector
    write (*,*) 'Value at average final vector: ', func(avgvector, fcall) 
    deallocate(avgvector)

    write (*,*) 'Best final vector: ', X%vectors(minloc(X%values), :)
    write (*,*) 'Value at best final vector: ', minval(X%values)

    write (*,*) 'Function calls: ', fcall

    deallocate(X%vectors, X%values)
  end subroutine run_de




  !initializes first generation of target vectors
  subroutine initialize(X, params, lowerbounds, upperbounds, fcall, func) 

    type(population), intent(out) :: X
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    real, external :: func
    integer :: i

    fcall = 0

    write (*,*) 'Begin DE'
    write (*,*) 'Parameters:'
    write (*,*) ' NP=', params%NP
    write (*,*) ' F=', params%F 
    write (*,*) ' Cr=', params%Cr 
    write (*,*) '-----------------------------'
    write (*,*) 'Generation: ', '1'

    allocate(X%vectors(params%NP, params%D), X%values(params%NP)) !deallocated at end of run_de

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
    if (all(U(:) .le. upperbounds) .and. all(U(:) .ge. lowerbounds)) then

       trialvalue = func(U(:), fcall)

       !when the trial vector is at least as good, use it for the next generation
       if (trialvalue .le. X%values(n)) then
          X%vectors(n,:) = U(:) 
          X%values(n) = trialvalue
          accept = accept + 1
       end if
    end if
  end subroutine selection


end module de
