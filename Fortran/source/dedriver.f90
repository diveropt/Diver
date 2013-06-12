module examples

use de
use detypes

implicit none

 integer, parameter :: NP=10, numgen=15, numciv=1, nDerived=2
 character (len=300) :: path='example_output/example'
 real(dp), parameter ::  Cr=0.9, tol = 1e-3			!recommend 0<F<1, 0<=Cr<=1
 real(dp), parameter, dimension(1) :: F=0.6
 real(dp), parameter, dimension(4) :: lowerbounds=-50.0	!boundaries of parameter space
 real(dp), parameter, dimension(4) :: upperbounds=50.0
 real(dp), parameter, dimension(4) :: ranges = upperbounds - lowerbounds
 real(dp), parameter :: dPrior = ranges(1)*ranges(2)*ranges(3)*ranges(4)

contains

!Functions to be minimized.  Assumed to be -ln(Likelihood)

real(dp) function constant(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall 
  logical, intent(out) :: quit
  logical, intent(in) :: validvector

  fcall = fcall + 1
  quit = .false.
  !-lnlike
  constant = 0.0_dp 
  if (.not. validvector) constant=huge(1.0_dp)
  !derived quantities (other functions of the parameters)
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]

end function constant


real(dp) function step(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  if (.not. validvector) then 
     step = huge(1.0_dp)
  else if (params(1) .gt. 0.0_dp) then
    step = 0.0_dp 
  else 
    step = 1.0_dp
  endif
  params(size(lowerbounds)+1:) = [2.0_dp*params(1),params(1)+params(2)]
  
end function step


real(dp) function linear(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  if (.not. validvector) then 
     linear = huge(1.0_dp)
  else if (params(1) .gt. 0.0_dp) then
    linear = params(1) 
  else 
    linear = 0.0_dp
  endif
  params(size(lowerbounds)+1:) = [2.0_dp*params(1),params(1)+params(2)]
  
end function linear


real(dp) function gauss1(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  gauss1 = params(1)*params(1) + params(2)*params(2)
  if (.not. validvector) gauss1=huge(1.0_dp)
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  
end function gauss1


real(dp) function gauss2(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  !if (fcall .lt. 30000) then 
    quit = .false.
  !else
  !  quit = .true.
  !endif
  gauss2 = (1.-params(1))*(1.-params(1)) + (5.-params(2))*(5.-params(2)) - 8.14897_dp
  params(size(lowerbounds)+1:) = [2.0_dp*params(1),params(1)+params(2)]
  if (.not. validvector) gauss2=huge(1.0_dp)
  
end function gauss2

real(dp) function gauss8(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  integer :: i  

  fcall = fcall + 1
  !if (fcall .lt. 30000) then 
    quit = .false.
  !else
  !  quit = .true.
  !endif
  gauss8=0.0_dp
  do i = 1,8
    gauss8 = gauss8+params(i)*params(i)
  enddo
  params(size(lowerbounds)+1:) = [2.0_dp*params(1),params(1)+params(2)]
  if (.not. validvector) gauss8=huge(1.0_dp)
  
end function gauss8

real(dp) function gauss4(params, fcall, quit, validvector)

  real(dp), dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  integer :: i  

  fcall = fcall + 1
  !if (fcall .lt. 30000) then 
    quit = .false.
  !else
  !  quit = .true.
  !endif
  gauss4=0.0_dp
  do i = 1,4
    gauss4 = gauss4+params(i)*params(i)
  enddo
  params(size(lowerbounds)+1:) = [2.0_dp*params(1),params(1)+params(2)]

  if (.not. validvector) gauss4=huge(1.0_dp)
  
end function gauss4

!Example prior distributions

!Flat prior distribution for all parameters
real(dp) function flatprior(X)

  real(dp), dimension(size(lowerbounds)), intent(in) :: X
  flatprior = 1.0_dp / dPrior

end function flatprior

end module examples



program dedriver !Tester for general differential evolution

  use examples

  implicit none

  call run_de(gauss4, flatprior, lowerbounds, upperbounds, path, nDerived=nDerived, doBayesian=.true., &
  resume=.false., Ztolerance=0.1_dp, lambda=0.8_dp, bndry=1)

end program dedriver
