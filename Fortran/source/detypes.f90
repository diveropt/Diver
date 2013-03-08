module detypes


implicit none

logical, parameter :: verbose = .false.                 !print verbose output


!FIXME: specify kinds in the rest of the program.
integer, parameter, public :: i4b = selected_int_kind(9)

integer, parameter, public :: sp = kind(1.0)
integer, parameter, public :: dp = kind(1.0d0)


type deparams                 				!differential evolution parameters 
                                                        ! (remember to expand io::save_state and io::resume if you expand this type)
   integer NP                 				!population size
   real, allocatable, dimension(:) ::  F 		!mutation scale factors
   integer Fsize                                        !number of entries in F
   real lambda                				!mutation scale factor for best-to-rand/current
   logical current            				!true: use current/best-to-current mutation
   real Cr                    				!crossover rate
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
   real :: tol						!tolerance in log-evidence
   real :: maxNodePop                                   !maximum population to allow in a cell before partitioning it
   logical :: calcZ      				!calculate evidence or not
   integer :: savefreq					!frequency with which to save progress
   integer :: mpirank					!rank of current process (0 if no MPI)
   integer :: mpipopchunk				!number of vectors for each process to work on (NP if no MPI)
end type codeparams

type population
  !add array of strings for names?
  real, allocatable, dimension(:,:) :: vectors 				!dimension(NP, D)
  real, allocatable, dimension(:) :: values, weights, multiplicities 	!dimension(NP)
  real, allocatable, dimension(:,:) :: derived				!dimension(NP, D_derived)
  real, allocatable, dimension(:) :: FjDE, CrjDE                        !dimension(NP)
end type population

end module detypes
