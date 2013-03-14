module evidence

use detypes
use post

implicit none

private
public updateEvidence, polishEvidence

integer, parameter :: samlun = 1, devolun=2

contains

  !Get posterior weights and update evidence on the fly
  subroutine updateEvidence(X, Z, Zmsq, Zerr, prior, oldsamples)
  
    type(population), intent(inout) :: X		!current generation
    real, intent(inout) :: Z, Zmsq, Zerr		!evidence, mean square of weights, error on evidence
    real, external :: prior 				!prior funtion
    real :: sampleratio, totsamples                     !ratio of old samples to total samples, total samples
    integer, intent(inout) :: oldsamples		!previous (running) number of samples
    integer :: inttotsamples				!total number of samples (integer)
    
    !Find weights for posterior pdf / evidence calculation
    call growTree(X,prior)
    
    !Find total number of samples and ratio to the old number
    inttotsamples = oldsamples + size(X%weights)
    totsamples = dble(inttotsamples)
    sampleratio = dble(oldsamples)/totsamples

    !Calculate multiplicity for outputting in chains
    X%multiplicities = X%weights*exp(-X%values)/totsamples

    !Update evidence
    Z = Z*sampleratio + sum(X%multiplicities)

    !Update the mean square of the weights and the standard deviation of the evidence
    Zmsq = Zmsq*sampleratio + sum(X%multiplicities*X%multiplicities*totsamples)
    Zerr = sqrt((Zmsq - Z*Z)/totsamples)
    
    !Update number of samples for next time
    oldsamples = inttotsamples

  end subroutine updateEvidence

  
  !Recalculate evidence and all posterior weights at the end of a run
  subroutine polishEvidence(Z, Zmsq, Zerr, prior, Nsamples, path, run_params, update)

    type(codeparams), intent(in) :: run_params
    real, intent(inout) :: Z, Zmsq, Zerr
    real, external :: prior 				
    real :: lnlike, multiplicity, vector(run_params%D), derived(run_params%D_derived)
    integer, intent(in) :: Nsamples
    integer :: filestatus, reclen, civ, gen, i
    character(len=*), intent(in) :: path
    character(len=31) :: formatstring
    character(len=1) :: LF
    logical, intent(in) :: update

    !organise the read/write format
    write(formatstring,'(A18,I4,A9)') '(2E20.9,2x,2I6,2x,', run_params%D+run_params%D_derived, 'E20.9,A1)'  
    reclen = 57 + 20*(run_params%D+run_params%D_derived)

    !open the chain file
    open(unit=samlun, file=trim(path)//'.sam', &
     iostat=filestatus, status='OLD', access='DIRECT', recl=reclen, form='FORMATTED')
    if (filestatus .ne. 0) stop ' Error opening .sam file. Quitting...' 
    
    !loop over the points in the sam file
    Z = 0.
    Zmsq = 0.
    do i = 1, Nsamples
      !read in each point 
      read(samlun,formatstring,rec=i) multiplicity, lnlike, civ, gen, vector, derived, LF
      !Could implement a skip out if this is the first generation (burn in), but this gen should not be in the sam file anyway
      !if (gen .eq. 1) cycle
      !use the tree to get a new weight for the point
      multiplicity = getWeight(vector,prior)*exp(-lnlike)/dble(Nsamples) 
      !save the new multiplicity of the point to disk
      if (update) write(samlun,formatstring,rec=i) multiplicity, lnlike, civ, gen, vector, derived, LF
      !add the contribution of the point with the new multiplicity to the evidence   
      Z = Z + multiplicity
      !add the contribution of the point with the new multiplicity to the error
      Zmsq = Zmsq + multiplicity*multiplicity
    enddo
    Zmsq = Zmsq*dble(Nsamples)
    Zerr = sqrt((Zmsq - Z*Z)/dble(Nsamples))   

    close(samlun)

  end subroutine polishEvidence


end module evidence
