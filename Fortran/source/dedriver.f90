module driverconsts

use de

implicit none

integer, parameter :: NP=10, numgen=50, numciv=500	!enforce NP>2
real, parameter :: F=0.7, Cr=0.9, tol = 0.01		!enforce 0<F<1, 0<=Cr<=1.  Set tol negative to forget posterior/evidence
real, parameter, dimension(2) :: lowerbounds=-50.0	!boundaries of parameter space
real, parameter, dimension(2) :: upperbounds=50.0
real, parameter, dimension(2) :: ranges = upperbounds - lowerbounds
real, parameter :: dPrior = ranges(1)*ranges(2)

contains

!function to be minimized.
function func(X, fcall)
  implicit none
  real, dimension(2), intent(in) :: X
  integer, intent(inout) :: fcall
  real func
  
  fcall = fcall + 1
  func = (1.0 - X(1))**2 + (5.0 - X(2))**2
end function func

!flat prior distribution for parameters of function
real function prior(X)
  implicit none
  real, dimension(2), intent(in) :: X
  prior = 1.0 / dPrior
end function prior

end module driverconsts



program dedriver !testing rand/1/bin differential evolution

use driverconsts
implicit none
call run_de(func, prior, lowerbounds, upperbounds, numciv, numgen, NP, F, Cr, tol)

end program dedriver


