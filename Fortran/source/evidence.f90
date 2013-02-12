module evidence

use detypes
use posterior

implicit none

contains

  !Get posterior weights and update evidence on the fly
  subroutine updateEvidence(X, Z, prior, oldsamples, newsamples)
  
    type(population), intent(inout) :: X		!current generation
    real, intent(inout) :: Z				!evidence
    real, external :: prior 				!prior funtion
    integer, intent(inout) :: oldsamples		!previous (running) number of samples
    integer, intent(in) :: newsamples 			!additional number of samples this time
    integer :: totsamples				!total number of samples
    
    !Find weights for posterior pdf / evidence calculation
    call getweights(X,prior)
    
    !Find total number of samples
    totsamples = oldsamples + newsamples

    !Calculate multiplicity for outputting in chains
    X%multiplicities = X%weights*exp(-X%values)/dble(totsamples)

    !Update evidence
    Z = Z*dble(oldsamples)/dble(totsamples) + sum(X%multiplicities)

    !Update number of samples for next time
    oldsamples = totsamples

  end subroutine updateEvidence

  
  !Recalculate evidence and all posterior weights at the end of a civilisation
  subroutine reassessEvidence(Z, oldZ)

    real :: Z, oldZ

  end subroutine reassessEvidence


end module evidence
