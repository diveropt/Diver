module post

use detypes

implicit none

!Entry in linked list of individuals (points in parameter space)
type Point
  type (Point), pointer :: next => null()
  real(dp), pointer :: vector(:) => null()
  real(dp), pointer :: weight => null()
end type

!Linked binary tree node type
type Node
  type (Node), pointer :: trunk => null(), branchA => null(), branchB => null()
  type (Point), pointer :: firstnewpt => null(), lastnewpt => null()
  real(dp) :: population = 0.0
  integer :: newpopulation = 0
  integer :: branchesDifferInDim = 0
  real(dp), allocatable :: upperbounds(:), lowerbounds(:)
end type Node
type (Node), pointer :: root

!Entry in linked list of binary tree nodes
type ListNode
  type (Node), pointer :: thisnode => null()
  type (ListNode), pointer :: next => null(), prev => null()
end type ListNode
type (ListNode), pointer :: firstListNode => null(), currentListNode => null()

!Entry in linked list of duplicate individuals (points in parameter space)
type Duplicate
  type (Duplicate), pointer :: next => null()
  type (Point), pointer :: point => null()
  type (Node), pointer :: stoppingNode => null()
end type
type (Duplicate), pointer :: firstDuplicate => null(), lastDuplicate => null()

real(dp), allocatable  :: ranges(:)             !Prior box side lengths
integer :: D = 0                                !Dimension of parameter space being scanned
integer :: totalCells = 0                       !Number of cells in binary space partitioning
logical :: debug = .false.                      !Debugging flag for posterior routines
real(dp) :: maxNodePop                          !Population at which to divide cells
real(dp), parameter :: duplicate_tol = 10.      !Scaling factor for epsilon in duplicate test

public iniTree, clearTree, growTree
private burnTree, climbTree, tendTree, growBranches, addToEndOfPtList, maxNodePop

contains

  subroutine iniTree(lowerbounds,upperbounds,maxpop)
  !Initialises the root of the tree and the starting point of the list of tree nodes
    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds      !boundaries of parameter space
    real(dp), intent(in) :: maxpop                                      !population at which to divide cells
    totalCells = 1
    D = size(lowerbounds)
    allocate(root, ranges(D))
    allocate(root%upperbounds(D), root%lowerbounds(D))
    root%upperbounds = upperbounds
    root%lowerbounds = lowerbounds
    ranges = upperbounds - lowerbounds
    allocate(firstListNode)
    maxNodePop = maxpop
  end subroutine iniTree


  subroutine clearTree
  !Start the tree over
    totalCells = 1
    call burnTree(root)
  end subroutine clearTree


  recursive subroutine burnTree(currentNode)
  !Burns off the tree everywhere upwards of currentNode

    type(Node), pointer, intent(in) :: currentNode

    if (currentNode%branchesDifferInDim .ne. 0) then
      !Node is not the tip of the tree, so has branches.  Get 'em.
      call burnTree(currentNode%branchA)
      call burnTree(currentNode%branchB)
      deallocate(currentNode%branchA%upperbounds,currentNode%branchA%lowerbounds)
      deallocate(currentNode%branchB%upperbounds,currentNode%branchB%lowerbounds)
      deallocate(currentNode%branchA,currentNode%branchB)
    endif

  end subroutine burnTree


  real(dp) function getWeight(vector, prior, context)
  !Finds the posterior weight of a single individual, using the existing tree structure

    use iso_c_binding, only: c_ptr

    type(c_ptr), intent(inout) :: context         !context pointer
    real(dp), intent (in), target :: vector(:)    !target vector
    procedure(PriorFunc) :: prior                 !prior pdf function
    type (Point), pointer :: individual
    real(dp), target :: weight

    allocate(individual)
    individual%vector => vector
    individual%weight => weight

    call climbTree(individual, root, justLooking=.true.)
    getWeight = individual%weight * real(totalCells, kind=dp) * prior(vector, context)

    deallocate(individual)

  end function getWeight


  subroutine growTree(X,prior,context)
  !Grows the tree with points in a new generation, and calculates their posterior weights

    use iso_c_binding, only: c_ptr

    type(population), target :: X                 !current generation of target vectors
    type(c_ptr), intent(inout) :: context         !context pointer
    procedure(PriorFunc) :: prior                 !prior pdf function
    type(Point), pointer :: individual            !pointer to a holder for an individual point in parameter space
    integer :: NP, i                              !size of generation, iteration variable

    nullify(individual)

    NP = size(X%values)

    currentListNode => firstListNode
    nullify(currentListNode%next)

    if (debug) write(*,*) 'current number of cells: ',totalCells
    if (debug) write(*,*) 'about to climb'

    !TODO parallelise the following loop (maybe not so easy)
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

    !Work through the linked list and check if additional branching is required,
    !deallocating each entry as we go, and adding more to the end where required
    call tendTree(firstListNode)

    !Repeat if any duplicates have been identified
    do while (associated(firstDuplicate))
      if (debug) write(*,*) 'Dealing with duplicates...'
      !Reset the list of nodes to be checked for branching
      currentListNode => firstListNode
      nullify(currentListNode%next)
      !Run any stragglers (duplicates) up the tree
      call straggleUpTree(firstDuplicate)
      !Work through any new entries in the linked list of nodes needing to be
      !checked for branching.
      call tendTree(firstListNode)
    end do

    !$OMP PARALLEL DO
    do i = 1, NP
      X%weights(i) = X%weights(i) * real(totalCells, kind=dp) * prior(X%vectors(i,:), context)
    enddo
    !$END OMP PARALLEL DO

    nullify(individual)

  end subroutine growTree


  recursive subroutine tendTree(workingListNode)
  !Steps through linked list of nodes potentially needing to be partitioned, and
  !does the partitioning

    type(ListNode), pointer :: workingListNode
    type(ListNode), pointer :: temp
    type(Point), pointer :: individual, temppt
    real(dp) :: nodeWeight
    integer :: i

    nullify (temp)
    nullify (individual)
    nullify (temppt)

    if (debug) write(*,*) 'attempting to grow tree'

    !Check if this is the first (=permanent) node in the list
    if (.not. associated(workingListNode%prev)) then
      !it is, so step on to the second node

      if (debug) write(*,*) 'this is the first node in the list'
      call tendTree(workingListNode%next)

    else
      if (debug) write(*,*) 'this is a later node'
      !it isn't.  Check if this node needs partitioning.

      if (workingListNode%thisNode%population .gt. maxNodePop) then
        !Current node needs to grow new branches

        call growBranches(workingListNode%thisNode)

      else

        if (debug) write(*,*) 'done growing here.'

        !No new branches required.  The weightings of the points in the
        !new population of this node are given by the weighting of the node
        nodeWeight = product(workingListNode%thisNode%upperbounds - workingListNode%thisNode%lowerbounds)
        individual => workingListNode%thisNode%firstnewpt

        if (debug) write(*,*) 'clearing individuals from this node'
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

      !Done with growing the tree for this node for now; set its new population to zero
      workingListNode%thisNode%newpopulation = 0

      !Remove this node from the linked list
      deallocate(workingListNode)

      !Check if this is the last node in the list
      if (associated(temp)) then
        !If not, step onwards to the next node in the linked list
        call tendTree(temp)
      endif

    endif

    nullify (temp)
    nullify (individual)
    nullify (temppt)

  end subroutine tendTree


  recursive subroutine straggleUpTree(straggler)
  !Sends stragglers climbing up the tree at a late stage, when
  !they cannot interfere with duplicates in the same
  !generation who had already climbed.

    type(Duplicate), pointer :: straggler

    !Send current straggler on up the tree to another node
    call climbTree(straggler%point, straggler%stoppingNode)

    !If this is not the last straggler, send the next one up
    if (associated(straggler%next)) call straggleUpTree(straggler%next)

    !Remove this duplicate from the linked list
    deallocate(straggler)

  end subroutine straggleUpTree


  recursive subroutine climbTree(individual,currentNode,justLooking)
  !Climbs a single level in the binary space partition tree from currentNode,
  !or establishes that currentNode is at the top of the tree.

    type(Point), pointer :: individual
    type(Node), pointer :: currentNode
    !Indication that the individual is here for a good time, not a long time
    logical, optional :: justLooking
    logical :: onlyLooking

    !Default is to assume individual is here to stay
    if (present(justLooking)) then
      onlyLooking = justLooking
    else
      onlyLooking = .false.
    endif

    if (currentNode%branchesDifferInDim .eq. 0) then
    !Node is tip of the tree and has no branches, so...

      if (debug) write(*,*) 'this is the canopy!'

      if (onlyLooking) then

        !Find the weight of this node and return it
        if (debug) write(*,*) 'just looking...'
        individual%weight = product(currentNode%upperbounds - currentNode%lowerbounds)

      else

        !increment the population of this node
        currentNode%population = currentNode%population + 1.0
        call addToEndOfPtList(individual,currentNode)
        if (debug) write(*,*) 'new population: ',currentNode%newpopulation

        !If this is the first time the current node has had a new point added to it, then
        !save the node to the linked list of nodes that need to be later checked for splitting
        if (.not. associated(currentNode%firstnewpt%next)) then
          if (debug) write(*,*) 'Adding node to list for later checking'
          !Create the next entry in the linked list of nodes
          allocate(currentListNode%next)
          !Set the next list entry to point back to the current list entry
          currentListNode%next%prev => currentListNode
          !Change current list entry to point to the newly-created one
          currentListNode => currentListNode%next
          !Save the location of the current node in the newly-created list entry
          currentListNode%thisNode => currentNode
        endif

      endif

    else
    !We need to move on to the next branch

      !Choose branch A or branch B depending on the location of the point in the magic (=splitting) dimension
      if (individual%vector(currentNode%branchesDifferInDim) .lt. &
       currentNode%branchA%upperbounds(currentNode%branchesDifferInDim)) then
        if (debug) write(*,*) 'moving up branch A...'
        call climbTree(individual,currentNode%branchA,justLooking=onlyLooking)

      else
        if (debug) write(*,*) 'moving up branch B...'
        call climbTree(individual,currentNode%branchB,justLooking=onlyLooking)

      endif

    endif

  end subroutine climbTree


  subroutine addToEndOfPtList(individual, currentNode)
  !Adds an additional individual point to the end of the linked list of new points associated with a node

    type(Point), pointer :: individual
    type(Node), pointer :: currentNode

    if (currentNode%newpopulation .eq. 0) then
      !If this is the first new point for this node, put it in the first entry in the list of points to save
      currentNode%firstnewpt => individual
      currentNode%lastnewpt => currentNode%firstnewpt

    else
      !Otherwise, make a new entry at the end of the list of points
      currentNode%lastnewpt%next => individual
      currentNode%lastnewpt => currentNode%lastnewpt%next

    endif

    !Terminate the list of points properly (required if individual was previously part of another list)
    nullify(currentNode%lastnewpt%next)

    !Increment the counter of new population points in the current node
    currentNode%newpopulation = currentNode%newpopulation + 1

  end subroutine addToEndofPtList


  recursive subroutine growBranches(currentNode)
  !Adds two new branches from a node, sorts the old population of the trunk node into them evenly,
  !and sends the new individuals in the trunk climbing up the new branches

    type(Node), pointer :: currentNode
    type(Point), pointer :: individual, temppt
    integer :: i, j, splitIndex(1)
    real(dp) :: ptRegistry(D, currentNode%newpopulation)
    logical :: dupeFlags(currentNode%newpopulation)

    nullify(individual)
    nullify(temppt)
    ptRegistry = 0.
    dupeFlags = .false.

    if (debug) write(*,*) 'Growing branches...'

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
    splitIndex = maxloc((currentNode%upperbounds-currentNode%lowerbounds)/ranges)
    currentNode%branchesDifferInDim = splitIndex(1)

    !Set branch A to the lower half of this partition, branch B to the upper half
    currentNode%branchA%upperbounds(splitIndex) = &
     0.5_dp*(currentNode%upperbounds(splitIndex)+currentNode%lowerbounds(splitIndex))
    currentNode%branchB%lowerbounds(splitIndex) = currentNode%branchA%upperbounds(splitIndex)

    !Give each new branch half the old population
    currentNode%branchA%population = 0.5_dp*(currentNode%population - real(currentNode%newpopulation, kind=dp))
    currentNode%branchB%population = currentNode%branchA%population

    !Send the points in the new population off up the tree
    individual => currentNode%firstnewpt

    do i = 1, currentNode%newpopulation

      !Need to save the next point in the list for next iteration in the loop,
      !as the current point will be appended to the end of a different list
      !after it is sent up the tree
      if (i .ne. currentNode%newpopulation) temppt => individual%next

      !Work out if the current individual is a duplicate of any already sent in this loop
      forall(j=1:i) dupeFlags(j) = all(abs(ptRegistry(:,j) - individual%vector) .le. &
       duplicate_tol*epsilon(individual%vector)*abs(individual%vector))

      !Don't allow duplicates to be sent climbing at this stage
      if (.not. any(dupeFlags)) then

        !Save details of current point
        ptRegistry(:,i) = individual%vector

        !Send current new point on up the tree to another node
        call climbTree(individual,currentNode)

      else

        if (debug) write(*,*) 'Duplicate found.'

        !Remove the individual from the list of points to be sent climbing
        nullify(individual%next)

        !Save the duplicate point for later
        if (.not. associated(firstDuplicate)) then
          !If this is the first duplicate point in the list, put it in the first entry in the list of duplicates
          allocate(firstDuplicate)
          lastDuplicate => firstDuplicate
        else
          !Otherwise, make a new entry at the end of the list of duplicate points
          allocate(lastDuplicate%next)
          lastDuplicate => lastDuplicate%next
        endif
        lastDuplicate%point => individual
        lastDuplicate%stoppingNode => currentNode

      endif

      !Set next point to try sending up the tree
      if (i .ne. currentNode%newpopulation) individual => temppt

    enddo

    !Set the population of the current node back to zero
    currentNode%population = 0.0_dp
    currentNode%newpopulation = 0

    !Clear the pointers to the linked list of new individuals in the current cell
    nullify(currentNode%firstnewpt, currentNode%lastnewpt)

    !Increment the counter of the total number of cells
    totalCells = totalCells + 1

    nullify(individual)
    nullify(temppt)

  end subroutine growBranches


end module post
