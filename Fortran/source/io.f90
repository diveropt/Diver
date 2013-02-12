module io

use detypes

implicit none

private
public save_all, io_begin

integer, parameter :: samlun = 1, devolun=2

contains


subroutine io_begin(path, civ, gen, Z, Zold, Nsamples, convcount, run_params, restart)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params
  logical, intent(in), optional :: restart
  integer :: filestatus  
 
  if (present(restart) .and. restart) then
    call resume(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)
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


subroutine save_all(X, path, civ, gen, Z, Zold, Nsamples, convcount, run_params, final)

  type(population), intent(in) :: X
  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params
  logical, intent(in), optional :: final

  if (.not. present(final) .or. (present(final) .and. .not. final)) call save_samples(X, path, civ, gen, run_params)  
  call save_state(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

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


subroutine save_state(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params
  integer :: filestatus
  character(len=12) :: formatstring
  
  !Save restart info
  open(unit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='OLD')
  if (filestatus .ne. 0) stop ' Error opening devo file.  Quitting...'

  write(devolun,'(2I6)') 	civ, gen					!current civilisation, generation
  write(devolun,'(2E16.5)') 	Z, Zold						!current evidence, evidence from previous gen
  write(devolun,'(I6)') 	Nsamples					!total number of independent samples so far
  write(devolun,'(I6)') 	convcount					!number of times delta ln Z < tol in a row so far

  write(devolun,'(I6)') 	run_params%DE%NP               			!population size
  write(formatstring,'(A1,I4,A6)') '(',size(run_params%DE%F),'E16.5)'
  if (.not. run_params%DE%jDE) then
    write(devolun,formatstring) run_params%DE%F	                                !mutation scale factors
  end if
  write(devolun,'(E16.5)') 	run_params%DE%lambda        			!mutation scale factor for best-to-rand/current
  write(devolun,'(L1)') 	run_params%DE%current            		!true: use current/best-to-current mutation
  write(devolun,'(E16.5)') 	run_params%DE%Cr            			!crossover rate
  write(devolun,'(L1)')  	run_params%DE%expon               		!when true, use exponential crossover (else use binomial)
  write(devolun,'(2I6)') 	run_params%D, run_params%D_derived		!dimension of parameter space (known from the bounds given); dimension of derived space
  write(devolun,'(2I6)') 	run_params%numciv, run_params%numgen		!maximum number of civilizations, generations
  write(devolun,'(E16.5)') 	run_params%tol					!tolerance in log-evidence
  write(devolun,'(I6)') 	run_params%convcountreq				!number of times delta ln Z < tol in a row for convergence
  write(devolun,'(L1)') 	run_params%calcZ				!calculate evidence or not
  write(devolun,'(I6)') 	run_params%savefreq				!frequency with which to save progress

  close(devolun)

end subroutine save_state


!Resumes from a previous run
subroutine resume(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params

  write(*,*) 'Restoring from previous run...'
  !FIXME mirror save_state, do some error-checking on overrides/disagreements with run_params

end subroutine resume


end module io
