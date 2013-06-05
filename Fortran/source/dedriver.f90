module examples

use de

implicit none

 integer, parameter :: NP=10, numgen=15, numciv=1, nDerived=2
 character (len=300) :: path='example_output/example'
 real, parameter ::  Cr=0.9, tol = 1e-3			!recommend 0<F<1, 0<=Cr<=1
 real, parameter, dimension(1) :: F=0.6
 real, parameter, dimension(4) :: lowerbounds=-50.0	!boundaries of parameter space
 real, parameter, dimension(4) :: upperbounds=50.0
 real, parameter, dimension(4) :: ranges = upperbounds - lowerbounds
 real, parameter :: dPrior = ranges(1)*ranges(2)*ranges(3)*ranges(4)

contains

!Functions to be minimized.  Assumed to be -ln(Likelihood)

real function constant(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall 
  logical, intent(out) :: quit
  logical, intent(in) :: validvector

  fcall = fcall + 1
  quit = .false.
  !-lnlike
  constant = 0. 
  if (.not. validvector) constant=huge(1.0)
  !derived quantities (other functions of the parameters)
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]

end function constant


real function step(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  if (.not. validvector) then 
     step = huge(1.0)
  else if (params(1) .gt. 0.0) then
    step = 0. 
  else 
    step = 1.
  endif
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  
end function step


real function linear(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  if (.not. validvector) then 
     linear = huge(1.0)
  else if (params(1) .gt. 0.0) then
    linear = params(1) 
  else 
    linear = 0.
  endif
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  
end function linear


real function gauss1(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  quit = .false.
  gauss1 = params(1)*params(1) + params(2)*params(2)
  if (.not. validvector) gauss1=huge(1.0)
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  
end function gauss1


real function gauss2(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  
  fcall = fcall + 1
  !if (fcall .lt. 30000) then 
    quit = .false.
  !else
  !  quit = .true.
  !endif
  gauss2 = (1.-params(1))*(1.-params(1)) + (5.-params(2))*(5.-params(2)) - 8.14897
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  if (.not. validvector) gauss2=huge(1.0)
  
end function gauss2

real function gauss8(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
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
  gauss8=0.d0
  do i = 1,8
    gauss8 = gauss8+params(i)*params(i)
  enddo
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]
  if (.not. validvector) gauss8=huge(1.0)
  
end function gauss8

real function gauss4(params, fcall, quit, validvector)

  real, dimension(size(lowerbounds)+nDerived), intent(inout) :: params
  integer, intent(inout) :: fcall
  logical, intent(out) :: quit
  logical, intent(in) :: validvector
  integer :: i  
  integer:: ierror
  real, dimension(4) :: totals=0., totalr

  fcall = fcall + 1
  !if (fcall .lt. 30000) then 
    quit = .false.
  !else
  !  quit = .true.
  !endif
  gauss4=0.d0
  do i = 1,4
    gauss4 = gauss4+params(i)*params(i)
  enddo
  params(size(lowerbounds)+1:) = [2.*params(1),params(1)+params(2)]

  if (.not. validvector) gauss4=huge(1.0)
  
end function gauss4

!Example prior distributions

!Flat prior distribution for all parameters
real function flatprior(X)

  real, dimension(size(lowerbounds)), intent(in) :: X
  flatprior = 1.0 / dPrior

end function flatprior

end module examples



program dedriver !Tester for general differential evolution

  use examples

  implicit none

  call run_de(gauss4, flatprior, lowerbounds, upperbounds, path, nDerived=nDerived, doBayesian=.true., &
  resume=.false., lambda=0.8, Ztolerance=0.1, bndry=1)

end program dedriver
