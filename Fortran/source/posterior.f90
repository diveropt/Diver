module posterior

use detypes

implicit none

!Linked binary tree node type
type Node
  type(Node), pointer :: trunk, branchA, branchB
  integer :: population = 0
  double precision, allocatable :: bounds(:,:)
end type Node

private
public getweights

contains 

!Calculates posterior pdf weights of points in a new generation
subroutine getweights(X,prior)

  type(population), intent(inout) :: X   !current generation of target vectors
  real prior
  external prior

  X%weights = 1.0

end subroutine


end module posterior
