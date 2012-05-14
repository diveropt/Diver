!Fortran 90 version of DE by Elinore Roebber
!needs convergence criteria, a better random number generator, kinds, and
!inputs for parameters, boundaries
!has external function

program de

implicit none

integer, parameter :: NP=10, numgen=50     !will need to enforce NP>2
real, parameter :: F=0.7, Cr=0.9    !will need to enforce 0<F<1, 0<=Cr<=1
real, parameter, dimension(2) :: lowerbounds=-50.0  !boundaries of parameter space
real, parameter, dimension(2) :: upperbounds=50.0
!will eventually want to input NP, F, Cr, bounds

integer, parameter :: D=size(upperbounds)
real, external :: func

type population
  integer :: generation
  real, dimension(NP, D) :: vectors
  real, dimension(NP) :: values
end type population

integer :: n, fcall !counts the times the objective function is called
real, dimension(D) :: result

type(population) :: X !target vectors

real, dimension(NP, D) :: V, U !donor, trial vectors





call initialize(X, lowerbounds, upperbounds, fcall, func)

!loop until convergence criteria satisfied. 
do n = 1, numgen-1

  V = rand1mutation(X)    !rand/1 mutation produces donor vectors
  U = bincrossover(X, V)  !bin crossover produces trial vectors
  call selection(X, U, fcall, func)    !choose next generation of target vectors

end do



!average over final generation for the result,
!also give best final vector?

print *, '============================='
print *, 'Number of generations: ', X%generation

result = sum(X%vectors, dim=1)/real(NP)
print *, 'Average final vector: ', result
print *, 'Value at average final vector: ', func(result, fcall) 

do n = 1, NP-1
  if (func(X%vectors(n,:), fcall) .ge. func(X%vectors(n+1,:), fcall)) then 
     result = X%vectors(n+1,:)
  end if
end do


print *, 'Function calls: ', fcall





contains 

  subroutine initialize(X, lowerbounds, upperbounds, fcall, func) !initializes first generation of target vectors
    !by choosing random components within boundaries for all NP D-dimensional vectors

    type(population), intent(out) :: X
    real, dimension(:), intent(in) :: lowerbounds, upperbounds
    integer, intent(inout) :: fcall
    integer i
    real, external :: func

    fcall = 0
    X%generation = 1

    print *, 'Begin DE'
    print *, 'Parameters:'
    print *, ' NP=', NP
    print *, ' F=', F 
    print *, ' Cr=', Cr 
    print *, '-----------------------------'
    print *, 'Generation: ', X%generation

    do i=1,NP
       call random_real(X%vectors(i,:), lowerbounds, upperbounds)
       X%values(i) = func(X%vectors(i,:), fcall)
       print *, X%vectors(i, :), '->', X%values(i)
    end do    
  end subroutine initialize



  function rand1mutation(X)
    type(population), intent(in) :: X  !current generation of target vectors
    real, dimension(NP, D) :: rand1mutation     !donor vectors
    integer :: r1, r2, r3, i

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    do i = 1, NP
       call random_int(r1, 1, NP)    !pick 1st vector from population
       do                            !pick 2nd vector; ensure vectors are distinct
          call random_int(r2, 1, NP) 
          if (r2 .ne. r1) exit
       end do
       do
          call random_int(r3, 1, NP) !pick 3rd vector
          if ((r3 .ne. r1) .and. (r3 .ne. r2)) exit
       end do
       !V = Xr1 + F*(Xr2 - Xr3)
       rand1mutation(i,:) = X%vectors(r1,:) + F*(X%vectors(r2,:) - X%vectors(r3,:))
   end do
  end function rand1mutation



  function bincrossover(X, V)
    !set each D-dimensional trial vector in U by comparing D random numbers with Cr and picking 
    !components from X or V as needed
    type(population), intent(in) :: X !current generation of target vectors
    real, dimension(NP, D), intent(in) :: V !donor vectors
    real, dimension(NP, D) :: bincrossover  !trial vectors created
    integer :: i, jrand            
    real, dimension(D) :: randj

    do i= 1, NP
      call random_int(jrand, 1, D)        !choose a guaranteed crossover for each vector.
      call random_number(randj)
      where (randj .le. Cr)
         bincrossover(i,:) = V(i,:) !use donor vector
      elsewhere
         bincrossover(i,:) = X%vectors(i,:) !use target vector
      end where
      bincrossover(i,jrand) = V(i,jrand) !guaranteed crossover of donor
   end do
  end function bincrossover



  subroutine selection(X,U, fcall, func)
    type(population), intent(inout) :: X
    real, dimension(NP, D), intent(in) :: U
    real :: trialvalue
    integer, intent(inout) :: fcall
    integer i
    real, external :: func
    
    !set new generation of X depending on whether f(X) or f(U) is lower, for each 
    !D-dimensional vector

   X%generation = X%generation + 1
   print *, '-----------------------------'
   print *, 'Generation: ', X%generation

    do i=1, NP
       !check that results stay within bounds. 'Brick wall'
       if (all(U(i,:) .le. upperbounds) .and. all(U(i,:) .ge. lowerbounds)) then
          !when the trial vector is at least as good, use it for the next generation
          trialvalue = func(U(i,:), fcall)
          if (trialvalue .le. X%values(i)) then
             X%vectors(i,:) = U(i,:) 
             X%values(i) = trialvalue
          end if                 
        end if
        print *, X%vectors(i, :), '->', X%values(i)
    end do
 
  end subroutine selection




!following functions, subroutines should eventually be moved elsewhere
!assorted random number generators.

  subroutine random_real(harvest, min, max) !choose array of random reals between min and max
    real, dimension(:), intent(out) :: harvest
    real, optional, dimension(:), intent(in) :: min, max !must be same size as harvest
    real, dimension(size(harvest)) :: range

    range = max - min
    call random_number(harvest)
    harvest = range*harvest + min
  end subroutine random_real


  subroutine random_int(harvest, min, max)        !choose a random integer between min and max, inclusive
    integer, intent(out) :: harvest
    integer, intent(in) :: min, max
    real :: range
    real r
   
    range = max - min
    call random_number(r)
    r = r*range + min
    harvest = nint(r)
  end subroutine random_int


end program de



