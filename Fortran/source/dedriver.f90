program dedriver !testing rand/1/bin differential evolution

use de

implicit none

integer, parameter :: NP=10, numgen=20             !enforce NP>2
real, parameter :: F=0.7, Cr=0.9                   !enforce 0<F<1, 0<=Cr<=1
real, parameter, dimension(2) :: lowerbounds=-50.0 !boundaries of parameter space
real, parameter, dimension(2) :: upperbounds=50.0
real, external :: func


call run_de(func, lowerbounds, upperbounds, numgen, NP, F, Cr)

end program dedriver




!function to be minimized.
function func(X, fcall)
  implicit none
  real, dimension(2), intent(in) :: X
  integer, intent(inout) :: fcall
  real func
  
  fcall = fcall + 1
  func = (1.0 - X(1))**2 + (5.0 - X(2))**2
end function func
