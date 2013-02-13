module evidence

use detypes
use posterior

implicit none

contains

  !Get posterior weights and update evidence on the fly
  subroutine updateEvidence(X, Z, prior, oldsamples)
  
    type(population), intent(inout) :: X		!current generation
    real, intent(inout) :: Z				!evidence
    real, external :: prior 				!prior funtion
    integer, intent(inout) :: oldsamples		!previous (running) number of samples
    integer :: totsamples				!total number of samples
    
    !Find weights for posterior pdf / evidence calculation
    call growTree(X,prior)
    
    !Find total number of samples
    totsamples = oldsamples + size(X%weights)

    !Calculate multiplicity for outputting in chains
    X%multiplicities = X%weights*exp(-X%values)/dble(totsamples)

    !Update evidence
    Z = Z*dble(oldsamples)/dble(totsamples) + sum(X%multiplicities)

    !Update number of samples for next time
    oldsamples = totsamples

  end subroutine updateEvidence

  
  !Recalculate evidence and all posterior weights at the end of a civilisation
  subroutine polishEvidence(Z, oldZ, prior)

    real :: Z, oldZ
    real, external :: prior 				!prior funtion

    !open the chain file
     !open

    !loop over the points in the chain
    
      !read in each point 
       !vector =
      !use the tree to get a new weight for the point
       !weight = getWeight(vector,prior)
      !save the new weight of the point to disk
       !write
      !add the contribution of the point with the new weight   
       !Z = Z + weight*exp(-value)

    !end loop
    !Z = Z/dble(npts)

  end subroutine polishEvidence


end module evidence
