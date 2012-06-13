module posterior

use detypes

implicit none

!Linked binary tree node type
type Node
  type(Node), pointer :: trunk => null(), branchA => null(), branchB => null()
  type(Point), pointer :: firstnewpt => null(), lastnewpt => null()
  real :: population = 0
  integer :: newpopulation = 0
  integer :: branchesDifferInDim = 0
  real, allocatable :: upperbounds(:), lowerbounds(:)
end type Node
type(Node), pointer :: root

!Entry in linked list of binary tree nodes
type ListNode
  type(Node), pointer :: thisnode => null()
  type(ListNode), pointer :: next => null(), prev => null()  
end type ListNode
type(ListNode), pointer :: firstListNode => null(), currentListNode => null()

!Entry in linked list of individuals (points in parameter space)
type Point
  type (Point), pointer :: next => null()
  real, pointer :: vector(:) => null()
  real, pointer :: weight => null()
end type

real, allocatable  :: ranges(:)			!Prior box side lengths
integer :: D = 0				!Dimension of parameter space being scanned
integer :: totalCells = 0			!Number of cells in binary space partitioning
real, parameter :: maxNodePop = 3.0 		!Population above which to do binary partitioning
logical :: debug = .false.			!Debugging flag for posterior routines

public initree, getweights
private climbTree, growTree, growBranches, addToEndOfPtList

contains 

  subroutine initree(lowerbounds,upperbounds)
  !Initialises the root of the tree and the starting point of the list of tree nodes
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space 
    totalCells = 1
    D = size(lowerbounds)
    allocate(root, ranges(D))
    allocate(root%upperbounds(D), root%lowerbounds(D))
    root%upperbounds = upperbounds
    root%lowerbounds = lowerbounds
    ranges = upperbounds - lowerbounds
    allocate(firstListNode)
  end subroutine initree


  subroutine getweights(X,prior)
  !Calculates evidence weights of points in a new generation

    type(population), target, intent(inout) :: X  !current generation of target vectors
    real prior					!prior pdf function
    external prior
    type(Point), pointer :: individual	 	!pointer to a holder for an individual point in parameter space
    integer :: NP, i				!size of generation, iteration variable

    nullify(individual)

    NP = size(X%values)	

    currentListNode => firstListNode
    nullify(currentListNode%next)

    if (debug) write(*,*) 'current number of cells: ',totalCells
    if (debug) write(*,*) 'about to climb'

    !FIXME check thread safety of parallelising following loop
    do i = 1, NP

      !Sort each of the new points into the existing tree structure, saving each
      !node where a new point is added in a linked list.
      if (debug) write(*,*) 'about to allocate individual',i,'of',NP
      if (debug) write(*,*) 'individual is already associated: ',associated(individual)

      allocate(individual)
      individual%vector => X%vectors(i,:)
      individual%weight => X%weights(i)
      
      if (debug) write(*,*) 'about to set off initial climber',i,'of',NP
      call climbTree(individual, root)
    enddo

    !Work through the linked list and check if addional branching is required, 
    !deallocating each entry as we go, and adding more to the end where required
    call growTree(firstListNode)

    !$OMP PARALLEL DO
    do i = 1, NP
      X%weights(i) = X%weights(i) * dble(totalCells) * prior(X%vectors(i,:))
    enddo
    !$END OMP PARALLEL DO

    nullify(individual)

  end subroutine getweights


  recursive subroutine climbTree(individual,currentNode)
  !Climbs a single level in the binary space partition tree from currentNode, 
  !or establishes that currentNode is at the top of the tree.

    type(Point), pointer, intent(in) :: individual
    type(Node), pointer, intent(in) :: currentNode

    if (currentNode%branchesDifferInDim .eq. 0) then	
    !Node is tip of the tree and has no branches, so...

      if (debug) write(*,*) 'this is the canopy!'
      !increment the population of this node
      currentNode%population = currentNode%population + 1.0
      if (debug) write(*,*) 'new population: ',currentNode%newpopulation
      call addToEndOfPtList(individual,currentNode)

      !If this is the first time the current node has had a new point added to it, then
      !save the node to the linked list of nodes that need to be later checked for splitting
      if (.not. associated(currentNode%firstnewpt%next)) then 
        !Create the next entry in the linked list of nodes
        allocate(currentListNode%next)
        !Set the next list entry to point back to the current list entry
        currentListNode%next%prev => currentListNode
        !Change current list entry to point to the newly-created one
        currentListNode => currentListNode%next
        !Save the location of the current node in the newly-created list entry 
        currentListNode%thisNode => currentNode

      endif

    else
    !We need to move on to the next branch  

      !Choose branch A or branch B depending on the location of the point in the magic (=splitting) dimension 
      if (individual%vector(currentNode%branchesDifferInDim) .lt. &
       currentNode%branchA%upperbounds(currentNode%branchesDifferInDim)) then
        if (debug) write(*,*) 'moving up branch A...'

        call climbTree(individual,currentNode%branchA)

      else 
        if (debug) write(*,*) 'moving up branch B...'
        call climbTree(individual,currentNode%branchB)

      endif

    endif

  end subroutine climbTree


  subroutine addToEndOfPtList(individual, currentNode)
  !Adds an additional individual point to the end of the linked list of new points associated with a node 

    type(Point), pointer, intent(in) :: individual
    type(Node), pointer, intent(in) :: currentNode

    if (currentNode%newpopulation .eq. 0) then
      !If this is the first new point for this node, put it in the first entry in the list of points to save
      currentNode%firstnewpt => individual
      currentNode%lastnewpt => currentNode%firstnewpt

    else
      !Otherwise, make a new entry at the end of the list of points 
      currentNode%lastnewpt%next => individual
      currentNode%lastnewpt => currentNode%lastnewpt%next

    endif

    !Terminate the list of points properly (required if individual was previosuly part of another list) 
    nullify(currentNode%lastnewpt%next)

    !Increment the counter of new population points in the current node
    currentNode%newpopulation = currentNode%newpopulation + 1

  end subroutine addToEndofPtList


  recursive subroutine growTree(workingListNode)
  !Steps through linked list of nodes potentially needing to be partitioned, and
  !does the partitioning

    type(ListNode), pointer, intent(inout) :: workingListNode
    type(ListNode), pointer :: temp 
    type(Point), pointer :: individual, temppt
    real :: nodeWeight
    integer :: i

    nullify (temp)
    nullify (individual)
    nullify (temppt)

    if (debug) write(*,*) 'attempting to grow tree'
  
    !Check if this is the first (=permanent) node in the list
    if (.not. associated(workingListNode%prev)) then
      !it is, so step on to the second node
      call growTree(workingListNode%next)

    else
      if (debug) write(*,*) 'this is a later node'
      !it isn't.  Check if this node needs partitioning.

      if (workingListNode%thisNode%population .gt. maxNodePop) then
        !Current node needs to grow new branches

        call growBranches(workingListNode%thisNode)

      else

        !No new branches required.  The weightings of the points in the 
        !new population of this node are given by the weighting of the node
        nodeWeight = product(workingListNode%thisNode%upperbounds - workingListNode%thisNode%lowerbounds)
        individual => workingListNode%thisNode%firstnewpt

        do i = 1, workingListNode%thisNode%newpopulation
          individual%weight = nodeWeight
          temppt => individual%next
          if (debug) write(*,*) 'deallocating individual'
          deallocate(individual)
          if (associated(temppt)) individual => temppt

        enddo

      endif

      !Save pointer to next node in the list
      temp => workingListNode%next

      !Done with growing the tree for this node; set its new population to zero and remove it from the linked list
      workingListNode%thisNode%newpopulation = 0
      deallocate(workingListNode)

      !Check if it is the last node in the list
      if (associated(temp)) then
        !If not, step onwards to the next node in the linked list
        call growTree(temp)
      endif

    endif

    nullify (temp)
    nullify (individual)
    nullify (temppt)

  end subroutine growTree


  recursive subroutine growBranches(currentNode)
  !Adds two new branches from a node, sorts the old population of the trunk node into them evenly,
  !and sends the new individuals in the trunk climbing up the new branches

    type(Node), pointer, intent(in) :: currentNode
    type(Point), pointer :: individual, temppt
    integer :: i, splitIndex(1)

    nullify(individual)
    nullify(temppt)

    !Create the two new nodes and set their boundaries to those of the trunk node
    allocate(currentNode%branchA, currentNode%branchB)
    allocate(currentNode%branchA%upperbounds(D), currentNode%branchA%lowerbounds(D))
    allocate(currentNode%branchB%upperbounds(D), currentNode%branchB%lowerbounds(D))
    currentNode%branchA%trunk => currentNode
    currentNode%branchB%trunk => currentNode
    currentNode%branchA%upperbounds = currentNode%upperbounds
    currentNode%branchA%lowerbounds = currentNode%lowerbounds
    currentNode%branchB%upperbounds = currentNode%upperbounds
    currentNode%branchB%lowerbounds = currentNode%lowerbounds

    !Work out which dimension to partition the current node in to make the two new nodes
    splitIndex = minloc((currentNode%upperbounds-currentNode%lowerbounds)/ranges)
    currentNode%branchesDifferInDim = splitIndex(1)
  
    !Set branch A to the lower half of this partition, branch B to the upper half
    currentNode%branchA%upperbounds(splitIndex) = &
     0.5*(currentNode%upperbounds(splitIndex)+currentNode%lowerbounds(splitIndex))
    currentNode%branchB%lowerbounds(splitIndex) = currentNode%branchA%upperbounds(splitIndex)

    !Give each new branch half the old population
    currentNode%branchA%population = 0.5*(currentNode%population - dble(currentNode%newpopulation))
    currentNode%branchB%population = currentNode%branchA%population

    !Send the points in the new population off up the tree
    individual => currentNode%firstnewpt

    if (currentNode%newpopulation .gt. 1) then
      do i = 1, currentNode%newpopulation - 1
 
        !Need to save the next point in the list for next iteration in the loop,  
        !as the current point will be appended to the end of a different list 
        !after it is sent up the tree
        temppt => individual%next
        !Send current new point on up the tree to another node

        call climbTree(individual,currentNode)
        !Set next point to send up the tree

        individual => temppt

     enddo

    endif
    call climbTree(individual,currentNode)
 
    !Set the population of the current node back to zero
    currentNode%population = 0.0
    currentNode%newpopulation = 0

    !Clear the pointers to the linked list of new individuals in the current cell
    nullify(currentNode%firstnewpt, currentNode%lastnewpt)

    !Increment the counter of the total number of cells
    totalCells = totalCells + 1

    nullify(individual)
    nullify(temppt)

  end subroutine growBranches


end module posterior
