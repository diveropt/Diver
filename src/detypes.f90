module detypes

use iso_c_binding, only: c_ptr

implicit none

integer, parameter, public :: dp = kind(1.0d0)          !definition of 'double precision' used throughout
!integer, parameter, public :: dp = kind(1.0)


!type definitions

type deparams                                           !differential evolution parameters
                                                        ! (remember to expand io::save_run_params and io::read_state if you expand this type)
   integer NP                                           !population size
   real(dp), allocatable, dimension(:) ::  F            !mutation scale factors
   integer Fsize                                        !number of entries in F
   real(dp) lambda                                      !mutation scale factor for best-to-rand/current
   logical current                                      !true: use current/best-to-current mutation
   real(dp) Cr                                          !crossover rate
   logical expon                                        !when true, use exponential crossover (else use binomial)
   integer bconstrain                                   !boundary constraints for selection
   logical jDE                                          !when true, use self-adaptive Cr and F for DE
   logical lambdajDE                                    !when true, use self-adaptive Cr, lambda, and F for DE
   logical removeDuplicates                             !when true, weeds out duplicate vectors within a generation
end type deparams

type codeparams                                         !code parameters (remember to expand io::save_run_params and io::read_state if you expand this type)
   type (deparams) DE                                   !differential evolution parameters
   real(dp), allocatable, dimension(:) :: lowerbounds   !lower bounds on parameter space
   real(dp), allocatable, dimension(:) :: upperbounds   !upper bounds on parameter space
   integer :: D, D_derived, D_discrete                  !dimension of parameter space; dimension of derived space, dimension of discrete parameter space
   integer, allocatable, dimension(:) :: discrete       !lists the discrete dimensions of parameter space (size D_discrete)
   logical :: partitionDiscrete                         !split the population evenly amongst discrete parameters and evolve separately
   integer, allocatable, dimension(:) :: repeat_scales  !population scale on which partitioned parameters repeat when partitionDiscrete = true
   integer :: subpopNP                                  !subpopulation in each partition of the population when partitionDiscrete = true
   integer :: numciv, numgen                            !maximum number of civilizations, generations
   real(dp) :: convthresh                               !threshold for convergence (smoothed fractional improvement in the mean population value)
   integer :: convsteps                                 !number of steps to smooth over for testing convergence
   real(dp) :: tol                                      !tolerance in log-evidence
   real(dp) :: maxNodePop                               !maximum population to allow in a cell before partitioning it
   logical :: calcZ                                     !calculate evidence or not
   logical :: outputSamples                             !write samples as output
   integer :: savefreq                                  !frequency with which to save progress
   integer :: mpirank                                   !rank of current process (0 if no MPI)
   integer :: mpiprocs                                  !number of processes running (1 if no MPI)
   integer :: mpipopchunk                               !number of vectors for each process to work on (NP if no MPI)
   integer :: init_population_strategy                  !initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
   logical :: discard_unfit_points                      !recalculate any trial vector whose fitness is above max_acceptable_value
   integer :: max_initialisation_attempts               !maximum number of times to try to find a valid vector for each slot in the initial population. Also applies to later generations if discard_unfit_points = .true.
   real(dp) :: max_acceptable_value                     !maximum fitness to accept for the initial generation if init_population_strategy > 0.
   type(c_ptr) :: context                               !context pointer
   integer :: verbose                                   !level of verbosity: 0=quiet, 1=basic, 2=civ-level info, 3=verbose, negative for mpirank!=0
   integer :: convergence_criterion                     !indicates which convergence criterion has been selected (see convergence.f90 for codes)
   integer :: seed                                      !base seed for random number generation; non-positive or absent means seed from the system clock
   real(dp) :: meanlike                                 !the normalized average fitness of the population for the last generation
   real(dp), allocatable, dimension(:) :: improvements  !fracdiff stored for convsteps most recent steps (to be smoothed over)
end type codeparams

type population
  !add array of strings for names?
  real(dp), allocatable, dimension(:,:) :: vectors                       !dimension(NP, D)
  real(dp), allocatable, dimension(:) :: values, weights, multiplicities !dimension(NP)
  real(dp), allocatable, dimension(:,:) :: vectors_and_derived           !dimension(NP, D+D_derived)
  real(dp), allocatable, dimension(:) :: FjDE, CrjDE, lambdajDE          !dimension(NP)
end type population



!interfaces for the likelihood and prior functions

abstract interface
   !the likelihood function to be minimised -- assumed to be -ln(likelihood)
   real(dp) function MinusLogLikeFunc(params, fcall, quit, validvector, context)
     use iso_c_binding, only: c_ptr
     import dp
     implicit none
     real(dp), dimension(:), intent(inout) :: params
     integer, intent(inout) :: fcall
     logical, intent(out) :: quit
     logical, intent(in) :: validvector
     type(c_ptr), intent(inout) :: context
   end function MinusLogLikeFunc
end interface

abstract interface
   !the prior function
   real(dp) function PriorFunc(X, context)
     use iso_c_binding, only: c_ptr
     import dp
     implicit none
     real(dp), dimension(:), intent(in) :: X
     type(c_ptr), intent(inout) :: context
   end function PriorFunc
end interface


end module detypes
