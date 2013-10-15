module crossover 

use detypes

implicit none

private
public gencrossover, init_CrjDE

real(dp), parameter :: tau=0.1                          !jDE control parameter from Brest et al. 2006

contains


  subroutine gencrossover(X, V, U, n, run_params, trialCr)

    type(population), intent(in) :: X               !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D), intent(in) :: V  !donor vectors
    real(dp), dimension(run_params%D), intent(out) :: U !trial vectors
    integer, intent(in) :: n                        !index of current target vector
    real(dp), intent(out) :: trialCr                    !self-adapting Cr

    if (run_params%DE%jDE) then 
       trialCr = newCr(X, n)
       U = bincrossover(X, V, n, run_params, trialCr)
    elseif (run_params%DE%expon) then
       U = expcrossover(X, V, n, run_params)
    else
       U = bincrossover(X, V, n, run_params)
    end if
    
  end subroutine gencrossover


  function bincrossover(X, V, n, run_params, trialCr)  !binomial crossover to create trial vectors

    type(population), intent(in) :: X              !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D), intent(in) :: V !donor vectors
    integer, intent(in) :: n                       !index of current target vector
    real(dp), intent(in), optional :: trialCr          !self-adapting Cr. Used for jDE
    real(dp), dimension(run_params%D) :: bincrossover  !trial vector created
    integer :: jrand           
    real(dp), dimension(run_params%D) :: randj
    real(dp) :: Cr

    if (present(trialCr)) then
       Cr = trialCr
    else
       Cr = run_params%DE%Cr
    endif

    call random_int(jrand, 1, run_params%D)        !choose a guaranteed crossover for each vector.
    call random_number(randj)

    where (randj .le. Cr)
       bincrossover(:) = V(:)                      !use donor vector
    elsewhere
       bincrossover(:) = X%vectors(n, :)           !use target vector
    end where

    bincrossover(jrand) = V(jrand)                 !guaranteed crossover of donor   
  end function bincrossover


  function expcrossover(X, V, n, run_params)

    type(population), intent(in) :: X   
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D), intent(in) :: V
    integer, intent(in) :: n               
    real(dp), dimension(run_params%D) :: expcrossover 

    integer L, j                                   !length of crossover, index for vectors
    real(dp) rand

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


  function newCr(X, n)
    type(population), intent(in) :: X
    integer, intent(in) :: n 
    real(dp) :: newCr
    real(dp) :: rand

    call random_number(rand)
    if (rand .lt. tau) then
       call random_number(newCr)
    else
       newCr = X%CrjDE(n)                        !use Cr from previous generation
    endif 
  end function newCr


  function init_CrjDE(size)
    integer, intent(in) :: size
    real(dp), dimension(size) :: init_CrjDE
    real(dp), dimension(size) :: rand

    call random_number(rand)
    init_CrjDE = rand

  end function init_CrjDE


  subroutine random_int(harvest, min, max) !choose a random integer between min and max, inclusive
    integer, intent(out) :: harvest
    integer, intent(in) :: min, max
    real(dp) :: range
    real(dp) :: r
   
    range = max - min
    call random_number(r)
    r = r*range + min
    harvest = nint(r)
  end subroutine random_int

end module crossover
