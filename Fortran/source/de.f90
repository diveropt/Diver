!Fortran 90 version of DE by Elinore Roebber
!needs convergence criteria, a better random number generator, kinds, and
!inputs for parameters, boundaries
!rearranged so that we pass individual vectors to mutation, crossover, selection

program de

implicit none

integer, parameter :: NP=10, numgen=20             !enforce NP>2
real, parameter :: F=0.7, Cr=0.9                   !enforce 0<F<1, 0<=Cr<=1
real, parameter, dimension(2) :: lowerbounds=-50.0 !boundaries of parameter space
real, parameter, dimension(2) :: upperbounds=50.0
!will eventually want to input NP, F, Cr, lowerbounds, upperbounds
!structure for the parameters?

integer, parameter :: D=size(upperbounds)          !dimension of the problem 
 
real, external :: func                             !function to be optimized.


type population
  real, dimension(NP, D) :: vectors
  real, dimension(NP) :: values
end type population

!fcall counts function calls, accept counts acceptance rate
integer :: gen, n, fcall, accept

real, dimension(D) :: result

type(population) :: X                              !population of target vectors
real, dimension(D) :: V, U                         !donor, trial vectors




call initialize(X, lowerbounds, upperbounds, fcall, func)
 
do gen = 2, numgen

   write (*,*) '-----------------------------'
   write (*,*) 'Generation: ', gen
   
   accept = 0


   do n=1, NP                    !for each member of the population
      V = rand1mutation(X, n)    !rand/1 mutation produces donor vectors
      U = bincrossover(X, V, n)  !bin crossover produces trial vectors
      call selection(X, U, n, fcall, func, accept) !choose next generation of target vectors
      write (*,*) X%vectors(n, :), '->', X%values(n)
   end do

 
   write (*,*) 'Acceptance rate: ', accept/real(NP)

   !check convergence criteria: if satisfied, exit loop

end do


write (*,*) '============================='
write (*,*) 'Number of generations: ', gen-1

result = sum(X%vectors, dim=1)/real(NP)
write (*,*) 'Average final vector: ', result
write (*,*) 'Value at average final vector: ', func(result, fcall) 

write (*,*) 'Best final vector: ', X%vectors(minloc(X%values), :)
write (*,*) 'Value at best final vector: ', minval(X%values)

write (*,*) 'Function calls: ', fcall





contains 

  subroutine initialize(X, lowerbounds, upperbounds, fcall, func) !initializes first generation of target vectors

    type(population), intent(out) :: X
    real, dimension(:), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    integer :: i
    real, external :: func

    fcall = 0

    write (*,*) 'Begin DE'
    write (*,*) 'Parameters:'
    write (*,*) ' NP=', NP
    write (*,*) ' F=', F 
    write (*,*) ' Cr=', Cr 
    write (*,*) '-----------------------------'
    write (*,*) 'Generation: ', '1'

    do i=1,NP
       call random_real(X%vectors(i,:), lowerbounds, upperbounds)
       X%values(i) = func(X%vectors(i,:), fcall)
       write (*,*) X%vectors(i, :), '->', X%values(i)
    end do    
  end subroutine initialize



  function rand1mutation(X, n)
    type(population), intent(in) :: X   !current generation of target vectors
    integer, intent(in) :: n            !index of current vector
    real, dimension(D) :: rand1mutation !donor vector
    integer :: r1, r2, r3

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    call random_int(r1, 1, NP)    !pick 1st vector from population
    do                            !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, NP) 
       if (r2 .ne. r1) exit
    end do
    do                            !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, NP) 
       if ((r3 .ne. r1) .and. (r3 .ne. r2)) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    rand1mutation(:) = X%vectors(r1,:) + F*(X%vectors(r2,:) - X%vectors(r3,:))
  end function rand1mutation



  function bincrossover(X, V, n)         !binomial crossover to create trial vectors

    type(population), intent(in) :: X    !current generation of target vectors
    real, dimension(D), intent(in) :: V  !donor vectors
    integer, intent(in) :: n             !index of current target vector
    real, dimension(D) :: bincrossover   !trial vector created
    integer :: jrand           
    real, dimension(D) :: randj

    call random_int(jrand, 1, D)         !choose a guaranteed crossover for each vector.
    call random_number(randj)
    where (randj .le. Cr)
       bincrossover(:) = V(:)            !use donor vector
    elsewhere
       bincrossover(:) = X%vectors(n, :) !use target vector
    end where
    bincrossover(jrand) = V(jrand)       !guaranteed crossover of donor
    
  end function bincrossover



  subroutine selection(X,U, n, fcall, func, accept) !select next generation of target vectors

    type(population), intent(inout) :: X
    real, dimension(D), intent(in) :: U
    real :: trialvalue
    integer, intent(inout) :: fcall, accept
    integer, intent(in) :: n
    real, external :: func

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


!assorted random number generators.

  subroutine random_real(harvest, min, max)  !choose array of random reals between min and max
    real, dimension(:), intent(out) :: harvest
    real, optional, dimension(:), intent(in) :: min, max !must be same size as harvest
    real, dimension(size(harvest)) :: range

    range = max - min
    call random_number(harvest)
    harvest = range*harvest + min
  end subroutine random_real


  subroutine random_int(harvest, min, max) !choose a random integer between min and max, inclusive
    integer, intent(out) :: harvest
    integer, intent(in) :: min, max
    real :: range
    real :: r
   
    range = max - min
    call random_number(r)
    r = r*range + min
    harvest = nint(r)
  end subroutine random_int


end program de



