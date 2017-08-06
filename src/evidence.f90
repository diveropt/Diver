module evidence

use detypes
use deutils
use post

implicit none

private
public updateEvidence, polishEvidence

integer :: rawlun, samlun, devolun

contains

  !Get posterior weights and update evidence on the fly
  subroutine updateEvidence(X, Z, Zmsq, Zerr, prior, context, oldsamples)
  
    use iso_c_binding, only: c_ptr

    type(population), intent(inout) :: X                !current generation
    real(dp), intent(inout) :: Z, Zmsq, Zerr            !evidence, mean square of weights, error on evidence
    procedure(PriorFunc), optional :: prior             !prior function
    real(dp) :: sampleratio, totsamples                 !ratio of old samples to total samples, total samples
    integer, intent(inout) :: oldsamples                !previous (running) number of samples
    integer :: inttotsamples                            !total number of samples (integer)
    type(c_ptr), intent(inout) :: context               !context pointer
    
    !Find weights for posterior pdf / evidence calculation
    call growTree(X,prior,context)
    
    !Find total number of samples and ratio to the old number
    inttotsamples = oldsamples + size(X%weights)
    totsamples = real(inttotsamples, kind=dp)
    sampleratio = real(oldsamples, kind=dp)/totsamples

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
  subroutine polishEvidence(Z, Zmsq, Zerr, prior, context, Nsamples, path, run_params, update)

    use iso_c_binding, only: c_ptr

    type(codeparams), intent(in) :: run_params
    real(dp), intent(inout) :: Z, Zmsq, Zerr
    procedure(PriorFunc) :: prior
    real(dp) :: lnlike, multiplicity, vector(run_params%D), vectors_and_derived(run_params%D+run_params%D_derived)
    integer, intent(in) :: Nsamples
    type(c_ptr), intent(inout) :: context
    integer :: filestatus, reclen_raw, reclen_sam, civ, gen, i
    character(len=*), intent(in) :: path
    character(len=31) :: formatstring_raw 
    character(len=31) :: formatstring_sam
    character(len=1)  :: LF
    logical, intent(in) :: update
    logical :: dosam

    !determine whether to bother with .sam file
    dosam = ((run_params%D_derived .ne. 0) .or. (size(run_params%discrete) .ne. 0))

    !organise the read/write formats
    write(formatstring_raw,'(A18,I4,A9)') '(2E20.9,2x,2I6,2x,', run_params%D, 'E20.9,A1)'  
    reclen_raw = 57 + 20*run_params%D
    if (dosam) then
       write(formatstring_sam,'(A18,I4,A9)') '(2E20.9,2x,2I6,2x,', run_params%D+run_params%D_derived, 'E20.9,A1)'
       reclen_sam = reclen_raw + 20*run_params%D_derived
    endif

    !open the chain files
    open(newunit=rawlun, file=trim(path)//'.raw', &
     iostat=filestatus, status='OLD', access='DIRECT', recl=reclen_raw, form='FORMATTED')
    if (filestatus .ne. 0) call quit_all_processes(' Error opening .raw file. Quitting...')
    if (dosam) then
      open(newunit=samlun, file=trim(path)//'.sam', &
       iostat=filestatus, status='OLD', access='DIRECT', recl=reclen_sam, form='FORMATTED')
      if (filestatus .ne. 0) call quit_all_processes(' Error opening .sam file. Quitting...')
    endif
    
    !loop over the points in the raw and sam files
    Z = 0.
    Zmsq = 0.
    do i = 1, Nsamples
      !read in each point 
      if (dosam) read(samlun,formatstring_sam,rec=i) multiplicity, lnlike, civ, gen, vectors_and_derived, LF
      read(rawlun,formatstring_raw,rec=i) multiplicity, lnlike, civ, gen, vector, LF
      !Could implement a skip out if this is the first generation (burn in), but this gen should not be in the raw/sam file anyway
      !if (gen .eq. 1) cycle
      !use the tree to get a new weight for the point
      multiplicity = getWeight(vector,prior,context)*exp(-lnlike)/real(Nsamples, kind=dp) 
      !save the new multiplicity of the point to disk
      if (update) then
        write(rawlun,formatstring_raw,rec=i) multiplicity, lnlike, civ, gen, vector, LF
        if (dosam) write(samlun,formatstring_sam,rec=i) multiplicity, lnlike, civ, gen, vectors_and_derived, LF
      endif
      !add the contribution of the point with the new multiplicity to the evidence   
      Z = Z + multiplicity
      !add the contribution of the point with the new multiplicity to the error
      Zmsq = Zmsq + multiplicity*multiplicity
    enddo
    Zmsq = Zmsq*real(Nsamples, kind=dp)
    Zerr = sqrt((Zmsq - Z*Z)/real(Nsamples, kind=dp))   

    close(rawlun)
    if (dosam) close(samlun)

  end subroutine polishEvidence


end module evidence
