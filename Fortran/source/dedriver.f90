module driverconsts

use de

implicit none

integer, parameter :: NP=10, numgen=15, numciv=1
real, parameter ::  Cr=0.9, tol = 1e-3			!recommend 0<F<1, 0<=Cr<=1.  Set tol negative to forget posterior/evidence
real, parameter, dimension(1) :: F=0.6
real, parameter, dimension(2) :: lowerbounds=-50.0	!boundaries of parameter space
real, parameter, dimension(2) :: upperbounds=50.0
real, parameter, dimension(2) :: ranges = upperbounds - lowerbounds
real, parameter :: dPrior = ranges(1)*ranges(2)

contains

!function to be minimized.
real function func(X, fcall)
  implicit none
  real, dimension(2), intent(in) :: X
  integer, intent(inout) :: fcall
  
  fcall = fcall + 1
  !if (X(1) .gt. 0.0) then
  !  func = 0. 
  !else 
  !  func = 1.0
  !endif
  !func=X(1)
  !func = X(1)**2+X(2)**2
  func = (1.0 - X(1))**2 + (5.0 - X(2))**2
end function func

!flat prior distribution for parameters of function
real function prior(X)
  implicit none
  real, dimension(2), intent(in) :: X
  prior = 1.0 / dPrior
end function prior

end module driverconsts



program dedriver !testing general differential evolution

use driverconsts
implicit none
call run_de(func, prior, lowerbounds, upperbounds, expon=.true.)

end program dedriver

