module mutation 

use detypes

implicit none

private
public rand1mutation

contains

    function rand1mutation(X, n, params)
    type(population), intent(in) :: X   !current generation of target vectors
    integer, intent(in) :: n            !index of current vector
    type(deparams), intent(in) :: params
    real, dimension(params%D) :: rand1mutation !donor vector
    integer :: r1, r2, r3

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    call random_int(r1, 1, params%NP)    !pick 1st vector from population
    do                            !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, params%NP) 
       if (r2 .ne. r1) exit
    end do
    do                            !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, params%NP) 
       if ((r3 .ne. r1) .and. (r3 .ne. r2)) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    rand1mutation(:) = X%vectors(r1,:) + params%F*(X%vectors(r2,:) - X%vectors(r3,:))
  end function rand1mutation


!add other mutation strategies


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


end module mutation
