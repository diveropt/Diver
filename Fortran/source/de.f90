!Has problems with loss of diversity.

!needs convergence criteria, a better random number generator, kinds, and
!inputs for parameters, boundaries, function (all currently hard-coded in)

program de

implicit none

integer, parameter :: NP=5, D=2     !will need to enforce NP>2, D>0
real, parameter :: F=0.7, Cr=0.9    !will need to enforce 0<F<1, 0<=Cr<=1
real, parameter, dimension(D) :: lowerbounds=-10.0  !boundaries of parameter space
real, parameter, dimension(D) :: upperbounds=10.0
!will eventually want to input NP, F, Cr, bounds and get D from function to optimize
integer n
real, dimension(D) :: result


type population
  integer generation
  real, dimension(NP, D) :: vectors
end type population

type(population) :: X, V, U !target, donor, trial vectors





call initialize(X, lowerbounds, upperbounds)

!loop until convergence criteria satisfied. 
do n = 1, 70
  print *, X
  V = rand1mutation(X)    !rand/1 mutation produces donor vectors
  print *,  V
  U = bincrossover(X, V)  !bin crossover produces trial vectors
  print *,  U
  call selection(X, U)    !choose next generation of target vectors
end do



!average over final generation for the result,
!also give best final vector
print *, X
print *, 'Generation: ', X%generation

result = sum(X%vectors, dim=1)/real(NP)
print *, 'Average final vector: ', result
print *, 'Value at average final vector: ', func(result) 

do n = 1, NP-1
  if (func(X%vectors(n,:)) .ge. func(X%vectors(n+1,:))) then 
     result = X%vectors(n+1,:)
  end if
end do
print *, 'Best final vector: ', result
print *, 'Value at best final vector: ', func(result)





contains 

  subroutine initialize(X, lowerbounds, upperbounds) !initializes first generation of target vectors
    !by choosing random components within boundaries for all NP D-dimensional vectors

    type(population), intent(out) :: X
    real, dimension(:), intent(in) :: lowerbounds, upperbounds
    integer i

    X%generation = 1

    do i=1,NP
       call random_real(X%vectors(i,:), lowerbounds, upperbounds)

    end do    
  end subroutine initialize



  function rand1mutation(X)
    type(population), intent(in) :: X  !current generation of target vectors
    type(population) rand1mutation     !donor vectors
    integer :: r1, r2, r3, i

    rand1mutation%generation = X%generation

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X; use forall?
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
       rand1mutation%vectors(i,:) = X%vectors(r1,:) + F*(X%vectors(r2,:) - X%vectors(r3,:))
   end do
  end function rand1mutation



  function bincrossover(X, V)
    !set each D-dimensional trial vector in U by comparing D random numbers with Cr and picking 
    !components from X or V as needed
    type(population), intent(in) :: X, V !current generation of target, donor vectors
    type(population) bincrossover        !trial vectors created
    integer :: i, jrand            
    real, dimension(D) :: randj

    bincrossover%generation = X%generation

    do i= 1, NP
      call random_int(jrand, 1, D)        !choose a guaranteed crossover for each vector.
      call random_number(randj)
      where (randj .le. Cr)
         bincrossover%vectors(i,:) = V%vectors(i,:) !use donor vector
      elsewhere
         bincrossover%vectors(i,:) = X%vectors(i,:) !use target vector
      end where
      bincrossover%vectors(i,jrand) = V%vectors(i,jrand) !guaranteed crossover of donor
   end do
  end function bincrossover



  subroutine selection(X,U)
    type(population), intent(inout) :: X
    type(population), intent(in) :: U
    integer i
    
    !set new generation of X depending on whether f(X) or f(U) is lower, for each 
    !D-dimensional vector

   X%generation = X%generation + 1

    do i=1, NP
       !check that results stay within bounds. 'Brick wall'
       if (all(U%vectors(i,:) .le. upperbounds) .and. all(U%vectors(i,:) .ge. lowerbounds)) then
          !when the trial vector is at least as good, use it for the next generation
          if (func(U%vectors(i,:)) .le. func(X%vectors(i,:))) then
             X%vectors(i,:) = U%vectors(i,:) 
          end if                 
        end if
    end do
 
  end subroutine selection




!following functions, subroutines should eventually be moved elsewhere


  !function to be minimized.
  function func(X)
    real, dimension(2), intent(in) :: X  !D=2.  Make sure this will work with above; no checks yet. 
    real func
    func = (1.0 - X(1))**2 + (5.0 - X(2))**2
  end function func


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
