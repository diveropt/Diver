module crossover 

use detypes

implicit none

private
public gencrossover

contains


  function gencrossover(X, V, n, run_params)

    type(population), intent(in) :: X          !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(run_params%D) :: gencrossover  !trial vector created

    if (run_params%DE%expon) then
       gencrossover = expcrossover(X, V, n, run_params)
    else
       gencrossover = bincrossover(X, V, n, run_params)
    end if
    
  end function gencrossover



  function bincrossover(X, V, n, run_params)       !binomial crossover to create trial vectors

    type(population), intent(in) :: X          !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(run_params%D) :: bincrossover  !trial vector created
    integer :: jrand           
    real, dimension(run_params%D) :: randj

    call random_int(jrand, 1, run_params%D)        !choose a guaranteed crossover for each vector.
    call random_number(randj)
    where (randj .le. run_params%DE%Cr)
       bincrossover(:) = V(:)                  !use donor vector
    elsewhere
       bincrossover(:) = X%vectors(n, :)       !use target vector
    end where
    bincrossover(jrand) = V(jrand)             !guaranteed crossover of donor   
  end function bincrossover




  function expcrossover(X, V, n, run_params)

    type(population), intent(in) :: X          !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                   !index of current target vector
    real, dimension(run_params%D) :: expcrossover  !trial vector created

    integer L, j                               !length of crossover, index for vectors
    real rand

    L=0
    do j=1, run_params%D                           !determine length of the crossover
       L = L + 1
       call random_number(rand)
       if (rand .gt. run_params%DE%Cr) exit   
    end do

    call random_int(j, 1, run_params%D)            !beginning of crossover
    
    !rewrite this in terms of the modulo function?
    if (j+L .gt. run_params%D) then
       expcrossover(:) = X%vectors(n,:)
       expcrossover(j:run_params%D) = V(j:run_params%D)
       expcrossover(1:j+L-run_params%D-1) = V(1:j+L-run_params%D-1)
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
