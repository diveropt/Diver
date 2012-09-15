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
    if (filestatus .ne. 0) stop 'Error creating output files. Quitting...'
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

  if (.not. present(final) .or. (present(final) .and. .not. final)) call save_samples(X, civ, gen)  
  call save_state(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

end subroutine save_all


subroutine save_samples(X, civ, gen)

  type(population), intent(in) :: X
  integer, intent(in) :: civ, gen

!.sam
!mult like civ gen params derived

end subroutine save_samples


subroutine save_state(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params


  !FIXME save restart info
!path.devo:
!civ, gen
!Z, Zold
!Nsamples
!convcount
!run_params

end subroutine save_state


!Resumes from a previous run
subroutine resume(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, convcount
  real, intent(in) :: Z, Zold
  type(codeparams), intent(in) :: run_params

  write(*,*) 'Restoring from previous run...'

end subroutine resume


end module io
