module detypes


implicit none

logical, parameter :: verbose = .false.                  !print verbose output

integer, parameter, public :: dp = kind(1.0d0)          !definition of 'double precision' used throughout


type deparams                 				!differential evolution parameters 
                                                        ! (remember to expand io::save_state and io::resume if you expand this type)
   integer NP                 				!population size
   real(dp), allocatable, dimension(:) ::  F 		!mutation scale factors
   integer Fsize                                        !number of entries in F
   real(dp) lambda                			!mutation scale factor for best-to-rand/current
   logical current            				!true: use current/best-to-current mutation
   real(dp) Cr                    			!crossover rate
   logical expon               				!when true, use exponential crossover (else use binomial)
   integer bconstrain                                   !boundary constraints for selection
   logical jDE                                          !when true, use self-adaptive Cr and F for DE
   logical removeDuplicates                             !when true, weeds out duplicate vectors within a generation
end type deparams

type codeparams                 			!code parameters (remember to expand io::save_state and io::resume if you expand this type)
   type (deparams) DE					!differential evolution parameters
   integer :: D, D_derived, D_discrete			!dimension of parameter space; dimension of derived space, dimension of discrete parameter space
   integer, allocatable, dimension(:) :: discrete       !lists the discrete dimensions of parameter space (size D_discrete)
   integer :: numciv, numgen				!maximum number of civilizations, generations
   real(dp) :: tol					!tolerance in log-evidence
   real(dp) :: maxNodePop                               !maximum population to allow in a cell before partitioning it
   logical :: calcZ      				!calculate evidence or not
   integer :: savefreq					!frequency with which to save progress
   integer :: mpirank					!rank of current process (0 if no MPI)
   integer :: mpipopchunk				!number of vectors for each process to work on (NP if no MPI)
end type codeparams

type population
  !add array of strings for names?
  real(dp), allocatable, dimension(:,:) :: vectors 			 !dimension(NP, D)
  real(dp), allocatable, dimension(:) :: values, weights, multiplicities !dimension(NP)
  real(dp), allocatable, dimension(:,:) :: vectors_and_derived		 !dimension(NP, D+D_derived)
  real(dp), allocatable, dimension(:) :: FjDE, CrjDE                     !dimension(NP)
end type population

end module detypes
