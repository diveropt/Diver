module examples

use de

implicit none

 integer, parameter :: NP=10, numgen=15, numciv=1, nDerived=2
 character (len=300) :: path='example_output/example'
 real, parameter ::  Cr=0.9, tol = 1e-3			!recommend 0<F<1, 0<=Cr<=1.  Set tol negative to forget posterior/evidence
 real, parameter, dimension(1) :: F=0.6
 real, parameter, dimension(2) :: lowerbounds=-50.0	!boundaries of parameter space
 real, parameter, dimension(2) :: upperbounds=50.0
 real, parameter, dimension(2) :: ranges = upperbounds - lowerbounds
 real, parameter :: dPrior = ranges(1)*ranges(2)

contains

!Functions to be minimized.  Assumed to be -ln(Likelihood)

real function constant(params, derived, fcall)

  real, dimension(size(lowerbounds)), intent(in) :: params
  real, dimension(nDerived), intent(out) :: derived
  integer, intent(inout) :: fcall 

  fcall = fcall + 1
  !-lnlike
  constant = 0. 
  !derived quantities (other functions of the parameters)
  derived = [2.*params(1),params(1)+params(2)]

end function constant


real function step(params, derived, fcall)

  real, dimension(size(lowerbounds)), intent(in) :: params
  real, dimension(nDerived), intent(out) :: derived
  integer, intent(inout) :: fcall
  
  fcall = fcall + 1
  if (params(1) .gt. 0.0) then
    step = 0. 
  else 
    step = 1.
  endif
  derived = [2.*params(1),params(1)+params(2)]
  
end function step


real function linear(params, derived, fcall)

  real, dimension(size(lowerbounds)), intent(in) :: params
  real, dimension(nDerived), intent(out) :: derived
  integer, intent(inout) :: fcall
  
  fcall = fcall + 1
  if (params(1) .gt. 0.0) then
    linear = params(1) 
  else 
    linear = 0.
  endif
  derived = [2.*params(1),params(1)+params(2)]
  
end function linear


real function gauss1(params, derived, fcall)

  real, dimension(size(lowerbounds)), intent(in) :: params
  real, dimension(nDerived), intent(out) :: derived
  integer, intent(inout) :: fcall
  
  fcall = fcall + 1
  gauss1 = params(1)*params(1) + params(2)*params(2)
  derived = [2.*params(1),params(1)+params(2)]
  
end function gauss1


real function gauss2(params, derived, fcall)

  real, dimension(size(lowerbounds)), intent(in) :: params
  real, dimension(nDerived), intent(out) :: derived
  integer, intent(inout) :: fcall
  
  fcall = fcall + 1
  gauss2 = (1.-params(1))*(1.-params(1)) + (5.-params(2))*(5.-params(2))
  derived = [2.*params(1),params(1)+params(2)]
  
end function gauss2


!Example prior distributions

!Flat prior distribution for all parameters
real function flatprior(X)

  real, dimension(2), intent(in) :: X
  flatprior = 1.0 / dPrior
end function flatprior


end module examples


program dedriver !Tester for general differential evolution

use examples

implicit none
  
  call run_de(gauss2, flatprior, lowerbounds, upperbounds, path, nDerived=nDerived, expon=.true.)

end program dedriver
