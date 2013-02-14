module io

use detypes

implicit none

private
public save_all, io_begin, resume

integer, parameter :: samlun = 1, devolun=2

contains


subroutine io_begin(path, civ, gen, Z, Zerr, Nsamples, run_params, restart)

  character(len=*), intent(in) :: path
  integer, intent(inout) :: civ, gen, Nsamples
  real, intent(inout) :: Z, Zerr
  type(codeparams), intent(inout) :: run_params
  logical, intent(in), optional :: restart
  integer :: filestatus  
 
  if (present(restart) .and. restart) then
    call resume(path, civ, gen, Z, Zerr, Nsamples, run_params)
  else
    !Create .sam and .devo files
    write(*,*) 'Creating DEvoPack output files at '//trim(path)//'.*'
    open(unit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='REPLACE')
    open(unit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='REPLACE')
    if (filestatus .ne. 0) stop ' Error creating output files. Quitting...'
    close(samlun)
    close(devolun)
  endif

end subroutine io_begin


subroutine save_all(X, path, civ, gen, Z, Zerr, Nsamples, run_params, final)

  type(population), intent(in) :: X
  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples
  real, intent(in) :: Z, Zerr
  type(codeparams), intent(in) :: run_params
  logical, intent(in), optional :: final

  if (.not. present(final) .or. (present(final) .and. .not. final)) call save_samples(X, path, civ, gen, run_params)  
  call save_state(path, civ, gen, Z, Zerr, Nsamples, run_params)

end subroutine save_all


subroutine save_samples(X, path, civ, gen, run_params)

  type(population), intent(in) :: X
  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen
  type(codeparams), intent(in) :: run_params
  integer :: filestatus, i
  character(len=28) :: formatstring

  open(unit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='OLD', POSITION='APPEND')
  if (filestatus .ne. 0) stop ' Error opening sam file.  Quitting...'
  write(formatstring,'(A18,I4,A6)') '(2E16.5,2x,2I6,2x,', run_params%D+run_params%D_derived, 'E16.5)'
  do i = 1, size(X%weights)
    write(samlun,formatstring) X%multiplicities(i), X%values(i), civ, gen, X%vectors(i,:), X%derived(i,:)
  enddo
  close(samlun)

end subroutine save_samples


subroutine save_state(path, civ, gen, Z, Zerr, Nsamples, run_params)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples
  real, intent(in) :: Z, Zerr
  type(codeparams), intent(in) :: run_params
  integer :: filestatus
  character(len=12) :: formatstring
  
  !Save restart info
  open(unit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='OLD')
  if (filestatus .ne. 0) stop ' Error opening devo file.  Quitting...'

  write(devolun,'(2I6)') 	civ, gen					!current civilisation, generation
  write(devolun,'(2E16.5)') 	Z, Zerr						!current evidence and uncertainty
  write(devolun,'(I6)') 	Nsamples					!total number of independent samples so far

  write(devolun,'(I6)') 	run_params%DE%NP               			!population size
  write(devolun,'(L1)') 	run_params%DE%jDE            			!true: use jDE
  write(devolun,'(I4)')         run_params%DE%Fsize                             !number of mutation scale factors

  if (run_params%DE%Fsize .ne. 0) then
    write(formatstring,'(A1,I4,A6)') '(',run_params%DE%Fsize,'E16.5)'
    write(devolun,formatstring) run_params%DE%F			 		!mutation scale factors
  endif 

  write(devolun,'(E16.5)') 	run_params%DE%lambda        			!mutation scale factor for best-to-rand/current
  write(devolun,'(L1)') 	run_params%DE%current            		!true: use current/best-to-current mutation
  write(devolun,'(E16.5)') 	run_params%DE%Cr            			!crossover rate
  write(devolun,'(L1)')  	run_params%DE%expon               		!when true, use exponential crossover (else use binomial)
  write(devolun,'(I6)')  	run_params%DE%bconstrain               		!boundary constraint to use
  write(devolun,'(2I6)') 	run_params%D, run_params%D_derived		!dim of parameter space (known from the bounds given); dim of derived space
  write(devolun,'(2I6)') 	run_params%numciv, run_params%numgen		!maximum number of civilizations, generations
  write(devolun,'(E16.5)') 	run_params%tol					!tolerance in log-evidence
  write(devolun,'(E16.5)') 	run_params%maxNodePop				!maximum population to allow in a cell before partitioning it
  write(devolun,'(L1)') 	run_params%calcZ				!calculate evidence or not
  write(devolun,'(I6)') 	run_params%savefreq				!frequency with which to save progress

  close(devolun)

end subroutine save_state


subroutine read_state(path, civ, gen, Z, Zerr, Nsamples, run_params)

  real, intent(out) :: Z, Zerr
  integer, intent(out) :: civ, gen, Nsamples
  integer :: filestatus
  character(len=*), intent(in) :: path
  character(len=12) :: formatstring
  type(codeparams), intent(out) :: run_params
  
  !Read in run info
  open(unit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='READ', status='OLD')
  if (filestatus .ne. 0) stop ' Error opening devo file.  Quitting...'

  read(devolun,'(2I6)') 	civ, gen					!current civilisation, generation
  read(devolun,'(2E16.5)') 	Z, Zerr						!current evidence and uncertainty
  read(devolun,'(I6)') 		Nsamples					!total number of independent samples so far

  read(devolun,'(I6)') 		run_params%DE%NP               			!population size
  read(devolun,'(L1)') 		run_params%DE%jDE            			!true: use jDE
  read(devolun,'(I4)')          run_params%DE%Fsize                             !number of mutation scale factors

  if (run_params%DE%Fsize .ne. 0) then
    read(devolun,'(I4)') run_params%DE%Fsize 
    write(formatstring,'(A1,I4,A6)') '(',run_params%DE%Fsize,'E16.5)'
    read(devolun,formatstring) run_params%DE%F			 		!mutation scale factors
  endif 

  read(devolun,'(E16.5)') 	run_params%DE%lambda        			!mutation scale factor for best-to-rand/current
  read(devolun,'(L1)') 		run_params%DE%current            		!true: use current/best-to-current mutation
  read(devolun,'(E16.5)') 	run_params%DE%Cr            			!crossover rate
  read(devolun,'(L1)')  	run_params%DE%expon               		!when true, use exponential crossover (else use binomial)
  read(devolun,'(I6)')  	run_params%DE%bconstrain               		!boundary constraint to use
  read(devolun,'(2I6)') 	run_params%D, run_params%D_derived		!dim of parameter space (known from the bounds given); dim of derived space
  read(devolun,'(2I6)') 	run_params%numciv, run_params%numgen		!maximum number of civilizations, generations
  read(devolun,'(E16.5)') 	run_params%tol					!tolerance in log-evidence
  read(devolun,'(E16.5)') 	run_params%maxNodePop				!maximum population to allow in a cell before partitioning it
  read(devolun,'(L1)') 		run_params%calcZ				!calculate evidence or not
  read(devolun,'(I6)') 		run_params%savefreq				!frequency with which to save progress

  close(devolun)

end subroutine read_state


!Resumes from a previous run
subroutine resume(path, civ, gen, Z, Zerr, Nsamples, run_params)

  character(len=*), intent(in) :: path
  integer, intent(inout) :: civ, gen, Nsamples
  real, intent(inout) :: Z, Zerr
  type(codeparams), intent(inout) :: run_params
  type(codeparams) :: run_params_restored

  write(*,*) 'Restoring from previous run...'
  !Read the run state
  call read_state(path, civ, gen, Z, Zerr, Nsamples, run_params_restored)
  !FIXME Do some error-checking on overrides/disagreements between run_params
  run_params = run_params_restored

end subroutine resume


end module io
