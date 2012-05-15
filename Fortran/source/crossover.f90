module crossover 

use detypes

implicit none

private
public bincrossover

contains



  function bincrossover(X, V, n, params)       !binomial crossover to create trial vectors

    type(population), intent(in) :: X          !current generation of target vectors
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(params%D) :: bincrossover  !trial vector created
    integer :: jrand           
    real, dimension(params%D) :: randj

    call random_int(jrand, 1, params%D)        !choose a guaranteed crossover for each vector.
    call random_number(randj)
    where (randj .le. params%Cr)
       bincrossover(:) = V(:)                  !use donor vector
    elsewhere
       bincrossover(:) = X%vectors(n, :)       !use target vector
    end where
    bincrossover(jrand) = V(jrand)             !guaranteed crossover of donor   
  end function bincrossover




  !add expcrossover




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

end module crossover
