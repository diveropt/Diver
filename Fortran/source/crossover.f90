module crossover 

use detypes

implicit none

private
public gencrossover

contains


  function gencrossover(X, V, n, params)

    type(population), intent(in) :: X          !current generation of target vectors
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(params%D) :: gencrossover  !trial vector created

    if (params%expon) then
       gencrossover = expcrossover(X, V, n, params)
    else
       gencrossover = bincrossover(X, V, n, params)
    end if
    
  end function gencrossover



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




  function expcrossover(X, V, n, params)

    type(population), intent(in) :: X          !current generation of target vectors
    type(deparams), intent(in) :: params
    real, dimension(params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(params%D) :: expcrossover  !trial vector created

    integer L, j                               !length of crossover, index for vectors
    real rand

    L=0
    do j=1, params%D                           !determine length of the crossover
       L = L + 1
       call random_number(rand)
       if (rand .gt. params%Cr) exit   
    end do

    call random_int(j, 1, params%D)            !beginning of crossover
    
    !rewrite this in terms of the modulo function?
    if (j+L .gt. params%D) then
       expcrossover(:) = X%vectors(n,:)
       expcrossover(j:params%D) = V(j:params%D)
       expcrossover(1:j+L-params%D-1) = V(1:j+L-params%D-1)
    else
       expcrossover(:) = X%vectors(n,:)
       expcrossover(j:j+L-1) = V(j:j+L-1)
    end if

  end function expcrossover




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
