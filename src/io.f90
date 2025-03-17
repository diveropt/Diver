module io

use detypes
use deutils
use converge

implicit none

private
public io_begin, save_all, save_run_params, resume

integer :: rawlun, samlun, devolun, rparamlun
real(dp), parameter :: Ftolscale = 100., Bndtolscale = 100.

contains


subroutine io_begin(gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path, restart)

  integer, intent(inout) :: gen, Nsamples, Nsamples_saved, fcall
  type(codeparams), intent(inout) :: run_params
  integer :: filestatus
  type(population), intent(inout) :: X, BF
  character(len=*), intent(in), optional :: path
  logical, intent(in), optional :: restart

  logical           :: restart_

  restart_ = .false.
  if (present(restart)) restart_ = restart

  if (restart_) then
    if (.not. present(path)) then
      call quit_de('Error: Resuming a Diver run requires the path argument to be set to the location of the previous run files.')
    endif
    call resume(path, gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF)
  else if (run_params%mpirank .eq. 0 .and. .not. run_params%disableIO) then
    if (.not. present(path)) then
      call quit_de('Error: The path argument must be set unless disableIO = true and not attempting to resume an old run.')
    endif
    !Create output .raw file
    if (run_params%outputRaw) then
      if (run_params%verbose .ge. 1) write(*,*) 'Creating Diver .raw file at '//trim(path)//'.raw'
      open(newunit=rawlun, file=trim(path)//'.raw', iostat=filestatus, action='WRITE', status='REPLACE')
      close(rawlun)
    endif
    !Create output .sam file only if there are discrete parameters or derived quantities to save.
    if (run_params%outputSam .and. ( (run_params%D_derived .ne. 0) .or. (size(run_params%discrete) .ne. 0) )) then
      if (run_params%verbose .ge. 1) write(*,*) 'Creating Diver .sam file at '//trim(path)//'.sam'
      open(newunit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='REPLACE')
      if (filestatus .ne. 0) call quit_all_processes(' Error creating .sam file. Quitting...')
      close(samlun)
    endif
  endif

end subroutine io_begin


subroutine save_all(X, BF, gen, Nsamples, Nsamples_saved, fcall, run_params, path, final)

  type(population), intent(in) :: X, BF
  integer, intent(inout) :: Nsamples_saved
  integer, intent(in) :: gen, Nsamples, fcall
  type(codeparams), intent(in) :: run_params
  character(len=*), intent(in), optional :: path
  logical, intent(in), optional :: final

  logical         :: final_

  final_ = .false.
  if (present(final)) final_ = final

  if (.not. final_) then
    Nsamples_saved = Nsamples_saved + run_params%DE%NP
    call save_samples(X, gen, run_params, path=path)
  endif
  call save_state(gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path=path)

end subroutine save_all


subroutine save_samples(X, gen, run_params, path)

  type(population), intent(in) :: X
  integer, intent(in) :: gen
  type(codeparams), intent(in) :: run_params
  integer :: filestatus, i
  character(len=28) :: formatstring_raw
  character(len=28) :: formatstring_sam
  character(len=*), intent(in), optional :: path

  if (run_params%disableIO) return

  if (run_params%outputRaw) then
    open(newunit=rawlun, file=trim(path)//'.raw', iostat=filestatus, action='WRITE', status='OLD', POSITION='APPEND')
    if (filestatus .ne. 0) call quit_all_processes(' Error opening raw file.  Quitting...')
    write(formatstring_raw,'(A18,I4,A6)') '(E20.9,2x,I6,2x,', run_params%D, 'E20.9)'
    do i = 1, size(X%values)
      write(rawlun,formatstring_raw) X%values(i), gen, X%vectors(i,:)
    enddo
    close(rawlun)
  endif

  if (run_params%outputSam) then
    if ((run_params%D_derived .ne. 0) .or. (size(run_params%discrete) .ne. 0)) then
      open(newunit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='OLD', POSITION='APPEND')
      if (filestatus .ne. 0) call quit_all_processes(' Error opening sam file.  Quitting...')
      write(formatstring_sam,'(A18,I4,A6)') '(E20.9,2x,I6,2x,', run_params%D+run_params%D_derived, 'E20.9)'
      do i = 1, size(X%values)
        write(samlun,formatstring_sam) X%values(i), gen, X%vectors_and_derived(i,:)
      enddo
      close(samlun)
    endif
  endif

end subroutine save_samples


subroutine save_run_params(run_params, path)

  type(codeparams), intent(in) :: run_params
  character(len=*), intent(in), optional :: path
  integer :: filestatus
  logical :: exists
  character(len=31) :: formatstring

  if (run_params%disableIO) return

  inquire(file=trim(path)//'.rparam',exist=exists)
  if (exists) then
     open(newunit=rparamlun, file=trim(path)//'.rparam', iostat=filestatus, action='WRITE', status='OLD')
  else
     open(newunit=rparamlun, file=trim(path)//'.rparam', iostat=filestatus, action='WRITE', status='REPLACE')
  endif
  if (filestatus .ne. 0) call quit_all_processes(' Error opening rparam file.  Quitting...')

  write(rparamlun,'(I6)')     run_params%DE%NP                          !population size
  write(rparamlun,'(L1)')     run_params%DE%jDE                         !true: use jDE
  write(rparamlun,'(L1)')     run_params%DE%lambdajDE                   !true: use jDE with self-adaptive lambda parameter
  write(rparamlun,'(I4)')     run_params%DE%Fsize                       !number of mutation scale factors

  if (run_params%DE%Fsize .ne. 0 .and. .not. run_params%DE%jDE) then
    write(formatstring,'(A1,I4,A6)') '(',run_params%DE%Fsize,'E20.9)'
    write(rparamlun,formatstring) run_params%DE%F                       !mutation scale factors
  endif

  write(rparamlun,'(E20.9)')  run_params%DE%lambda                      !mutation scale factor for best-to-rand/current
  write(rparamlun,'(L1)')     run_params%DE%current                     !true: use current/best-to-current mutation
  write(rparamlun,'(E20.9)')  run_params%DE%Cr                          !crossover rate
  write(rparamlun,'(L1)')     run_params%DE%expon                       !when true, use exponential crossover (else use binomial)
  write(rparamlun,'(I6)')     run_params%DE%bconstrain                  !boundary constraint to use
  write(rparamlun,'(2I6)')    run_params%D, run_params%D_derived        !dim of parameter space (known from the bounds given); dim of derived space
  write(formatstring,'(A1,I4,A6)') '(',run_params%D,'E20.9)'
  write(rparamlun,formatstring) run_params%lowerbounds                  !lower bounds of prior box
  write(rparamlun,formatstring) run_params%upperbounds                  !upper bounds of prior box
  write(rparamlun,'(I6)')     run_params%D_discrete                     !dimenension of discrete parameter space
  if (run_params%D_discrete .ne. 0) then
     write(formatstring,'(A1,I4,A3)') '(',run_params%D_discrete,'I6)'
     write(rparamlun,formatstring) run_params%discrete                  !discrete dimensions
     write(rparamlun,'(L1)')  run_params%partitionDiscrete              !split the population amongst discrete parameters and evolve separately
     if (run_params%partitionDiscrete) then
        write(rparamlun,formatstring) run_params%repeat_scales          !scales on which partitioned parameters repeat
        write(rparamlun,'(I6)') run_params%subpopNP                     !subpopulation NP for partitioned parameters
     endif
  endif
  write(rparamlun,'(I6)')     run_params%numgen                         !maximum number of generations
  write(rparamlun,'(E20.9)')  run_params%convthresh                     !threshold for gen-level convergence
  write(rparamlun,'(I6)')     run_params%convsteps                      !number of steps to smooth over when checking convergence
  write(rparamlun,'(L1)')     run_params%disableIO                      !disable all IO or not
  write(rparamlun,'(L1)')     run_params%outputRaw                      !output raw parameter samples to a .raw file or not
  write(rparamlun,'(L1)')     run_params%outputSam                      !output rounded and derived parameter samples to a .sam file or not
  write(rparamlun,'(I6)')     run_params%savefreq                       !frequency with which to save progress
  write(rparamlun,'(L1)')     run_params%DE%removeDuplicates            !true: remove duplicate vectors in a generation
  write(rparamlun,'(I6)')     run_params%verbose                        !amount of output to print to the screen
  write(rparamlun,'(I6)')     run_params%convergence_criterion          !indicates which convergence criterion has been selected (see convergence.f90 for codes)

  close(rparamlun)

end subroutine save_run_params


subroutine save_state(gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path)

  integer, intent(in) :: gen, Nsamples, Nsamples_saved, fcall
  type(codeparams), intent(in) :: run_params
  integer :: filestatus
  logical :: exists
  character(len=31) :: formatstring
  type(population), intent(in) :: X, BF
  character(len=*), intent(in), optional :: path

  if (run_params%disableIO) return

  !Save restart info
  inquire(file=trim(path)//'.devo',exist=exists)
  if (exists) then
     open(newunit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='OLD')
  else
     open(newunit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='REPLACE')
  endif
  if (filestatus .ne. 0) call quit_all_processes(' Error opening devo file.  Quitting...')

  write(devolun,'(I10)')      gen                                       !current generation
  write(devolun,'(3I10)')     Nsamples, Nsamples_saved, fcall           !total number of independent samples so far, num saved, num function calls

  write(devolun,'(E20.9)')    BF%values(1)                              !current best-fit
  write(formatstring,'(A1,I4,A6)') '(',run_params%D,'E20.9)'
  write(devolun,formatstring) BF%vectors(1,:)                           !current best-fit vector
  write(formatstring,'(A1,I4,A6)') '(',run_params%D+run_params%D_derived,'E20.9)'
  write(devolun,formatstring) BF%vectors_and_derived(1,:)               !reprocessed vector and derived parameters at current best fit

  write(formatstring,'(A1,I8,A6)') '(',run_params%DE%NP,'E20.9)'
  write(devolun,formatstring) X%values                                  !current population fitnesses
  write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*run_params%D,'E20.9)'
  write(devolun,formatstring) X%vectors                                 !current population
  write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived),'E20.9)'
  write(devolun,formatstring) X%vectors_and_derived                     !current reprocessed vector and derived values

  if (run_params%DE%jDE) then                                           !for self-adaptive F, Cr, optional lambda
    write(formatstring,'(A1,I16,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived)*run_params%DE%Fsize,'E20.9)'
    write(devolun,formatstring) X%FjDE                                  !current population F values
    write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived),'E20.9)'
    write(devolun,formatstring) X%CrjDE                                 !current population Cr values
    if (run_params%DE%lambdajDE) then
       write(devolun, formatstring) X%lambdajDE                         !current population lambda values
    end if
  end if

  if (run_params%convergence_criterion == meanimprovement) then
     write(devolun,'(E20.9)')    run_params%meanlike                    !the average fitness of the population for the last generation
     write(formatstring,'(A1,I4,A6)') '(',run_params%convsteps,'E20.9)'
     write(devolun,formatstring) run_params%improvements                !fractional diff in the mean, for convsteps most recent steps
  endif

  close(devolun)

end subroutine save_state


subroutine read_state(path, gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

  integer, intent(out) :: gen, Nsamples, Nsamples_saved, fcall
  integer :: filestatus, inNP
  logical :: exists
  character(len=*), intent(in) :: path
  character(len=31) :: formatstring
  type(codeparams), intent(inout) :: run_params
  type(population), intent(inout) :: X, BF

  !Read in run parameters
  inquire(file=trim(path)//'.rparam',exist=exists)
  if (.not. exists) call quit_all_processes(trim(path)//'.rparam does not exist. Cannot resume Diver.')
  open(newunit=rparamlun, file=trim(path)//'.rparam', iostat=filestatus, action='READ', status='OLD')
  if (filestatus .ne. 0) call quit_all_processes(' Error opening rparam file.  Quitting...')

  read(rparamlun,'(I6)')     inNP                                       !population size
  if (run_params%DE%NP .ne. inNP) then
     write(*,*) 'Error: NP differs in current and previous run. '
     write(*,*) 'Current:  ',run_params%DE%NP
     write(*,*) 'Previous: ',inNP
     call quit_all_processes('Please modify NP and try again.')
  endif
  run_params%DE%NP = inNP
  read(rparamlun,'(L1)')     run_params%DE%jDE                          !true: use jDE
  read(rparamlun,'(L1)')     run_params%DE%lambdajDE                    !true: use jDE with self-adaptive lambda
  read(rparamlun,'(I4)')     run_params%DE%Fsize                        !number of mutation scale factors

  if (run_params%DE%Fsize .ne. 0 .and. .not. run_params%DE%jDE) then
    allocate(run_params%DE%F(run_params%DE%Fsize))
    write(formatstring,'(A1,I4,A6)') '(',run_params%DE%Fsize,'E20.9)'
    read(rparamlun,formatstring) run_params%DE%F                        !mutation scale factors
  endif

  read(rparamlun,'(E20.9)')  run_params%DE%lambda                       !mutation scale factor for best-to-rand/current
  read(rparamlun,'(L1)')     run_params%DE%current                      !true: use current/best-to-current mutation
  read(rparamlun,'(E20.9)')  run_params%DE%Cr                           !crossover rate
  read(rparamlun,'(L1)')     run_params%DE%expon                        !when true, use exponential crossover (else use binomial)
  read(rparamlun,'(I6)')     run_params%DE%bconstrain                   !boundary constraint to use
  read(rparamlun,'(2I6)')    run_params%D, run_params%D_derived         !dim of parameter space (known from the bounds given); dim of derived space
  write(formatstring,'(A1,I4,A6)') '(',run_params%D,'E20.9)'
  allocate(run_params%lowerbounds(run_params%D), run_params%upperbounds(run_params%D))
  read(rparamlun,formatstring) run_params%lowerbounds                   !lower bounds of prior box
  read(rparamlun,formatstring) run_params%upperbounds                   !upper bounds of prior box
  read(rparamlun,'(I6)')     run_params%D_discrete                      !dimension of discrete parameter space
  if (run_params%D_discrete .gt. 0) then
     allocate(run_params%discrete(run_params%D_discrete))
     write(formatstring,'(A1,I4,A6)') '(',run_params%D_discrete,'I6)'
     read(rparamlun,formatstring) run_params%discrete                   !discrete dimensions in parameter sapce
     read(rparamlun,'(L1)')  run_params%partitionDiscrete               !split the population amongst discrete parameters and evolve separately
     if (run_params%partitionDiscrete) then
        read(rparamlun,formatstring) run_params%repeat_scales           !scales on which partitioned parameters repeat
        read(rparamlun,'(I6)') run_params%subpopNP                      !subpopulation NP for partitioned parameters
    endif
  else
     allocate(run_params%discrete(0))
  endif
  read(rparamlun,'(I6)')     run_params%numgen                          !maximum number of generations
  read(rparamlun,'(E20.9)')  run_params%convthresh                      !threshold for gen-level convergence
  read(rparamlun,'(I6)')     run_params%convsteps                       !number of steps to smooth over when checking convergence
  read(rparamlun,'(L1)')     run_params%disableIO                       !disable all IO or not
  read(rparamlun,'(L1)')     run_params%outputRaw                       !output raw parameter samples to a .raw file or not
  read(rparamlun,'(L1)')     run_params%outputSam                       !output rounded and derived parameter samples to a .sam file or not
  read(rparamlun,'(I6)')     run_params%savefreq                        !frequency with which to save progress
  read(rparamlun,'(L1)')     run_params%DE%removeDuplicates             !true: remove duplicate vectors in a generation
  read(rparamlun,'(I6)')     run_params%verbose                         !amount of output to print to the screen
  read(rparamlun,'(I6)')     run_params%convergence_criterion           !indicates which convergence criterion has been selected (see convergence.f90 for codes)

  close(rparamlun)

  !Read in run status info
  inquire(file=trim(path)//'.devo',exist=exists)
  if (.not. exists) call quit_all_processes(trim(path)//'.devo does not exist. Cannot resume Diver.')
  open(newunit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='READ', status='OLD')
  if (filestatus .ne. 0) call quit_all_processes(' Error opening devo file.  Quitting...')

  read(devolun,'(I10)')      gen                                        !current generation
  read(devolun,'(3I10)')     Nsamples, Nsamples_saved, fcall            !total number of independent samples so far, num saved, num function calls

  read(devolun,'(E20.9)')    BF%values(1)                               !current best-fit
  write(formatstring,'(A1,I4,A6)') '(',run_params%D,'E20.9)'
  read(devolun,formatstring) BF%vectors(1,:)                            !current best-fit vector
  write(formatstring,'(A1,I4,A6)') '(',run_params%D+run_params%D_derived,'E20.9)'
  read(devolun,formatstring) BF%vectors_and_derived(1,:)                !reprocessed vector and derived parameters at current best fit

  write(formatstring,'(A1,I8,A6)') '(',run_params%DE%NP,'E20.9)'
  read(devolun,formatstring) X%values                                   !current population fitnesses
  write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*run_params%D,'E20.9)'
  read(devolun,formatstring) X%vectors                                  !current population
  write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived),'E20.9)'
  read(devolun,formatstring) X%vectors_and_derived                      !current reprocessed vector and derived values

  if (run_params%DE%jDE) then                                           !for self-adaptive F, Cr, optional lambda
    write(formatstring,'(A1,I16,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived)*run_params%DE%Fsize,'E20.9)'
    read(devolun,formatstring) X%FjDE                                   !current population F values
    write(formatstring,'(A1,I12,A6)') '(',run_params%DE%NP*(run_params%D+run_params%D_derived),'E20.9)'
    read(devolun,formatstring) X%CrjDE                                  !current population Cr values
    if (run_params%DE%lambdajDE) then
       read(devolun, formatstring) X%lambdajDE                          !current population lambda values
    end if
  end if

  if (run_params%convergence_criterion == meanimprovement) then
     read(devolun,'(E20.9)') run_params%meanlike                        !the average fitness of the population for the last generation
     write(formatstring,'(A1,I4,A6)') '(',run_params%convsteps,'E20.9)'
     allocate(run_params%improvements(run_params%convsteps))
     read(devolun,formatstring) run_params%improvements                 !fractional diff in the mean, for convsteps most recent steps
  endif

  close(devolun)

end subroutine read_state


!Resumes from a previous run
subroutine resume(path, gen, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

  character(len=*), intent(in) :: path
  integer, intent(inout) :: gen, Nsamples, Nsamples_saved, fcall
  integer :: reclen, filestatus, i, j, passoverlen
  character(len=31) :: formatstring
  character(len=1) :: LF
  type(codeparams), intent(inout) :: run_params
  type(codeparams) :: run_params_restored
  type(population), intent(inout) :: X, BF
  type(population) :: Y

  if (run_params%verbose .ge. 1) write(*,*) 'Restoring from previous run...'

  !Read the run state
  run_params_restored%DE%NP = run_params%DE%NP
  call read_state(path, gen, Nsamples, Nsamples_saved, fcall, run_params_restored, X, BF)
  if (run_params_restored%convergence_criterion == meanimprovement) then
    run_params%meanlike = run_params_restored%meanlike
    passoverlen = min(run_params%convsteps,run_params_restored%convsteps)
    run_params%improvements(1:passoverlen) = run_params_restored%improvements(1:passoverlen)
    if (passoverlen .lt. run_params%convsteps) run_params%improvements(passoverlen+1:) = 1.0_dp
  endif

  !Do some error-checking on overrides/disagreements between run_params
  if (run_params%D .ne. run_params_restored%D) &
   call quit_de('Restored and new runs have different dimensionality.')
  if (run_params%D_derived .ne. run_params_restored%D_derived) &
   call quit_de('Restored and new runs have different number of derived params.')
  if (run_params%D_discrete .ne. run_params_restored%D_discrete) &
   call quit_de('Restored and new runs have different number of discrete parameters.')
  if ( any(run_params%discrete .ne. run_params_restored%discrete)) &
   call quit_de('Restored and new runs have different discrete parameters.')

  !Check for full generations
  if (mod(Nsamples_saved,run_params%DE%NP) .ne. 0) then
    call quit_de('Error: resumed run does not contain only full generations - file likely corrupted.')
  endif

  !Make sure we haven't already passed the number of gens
  if (gen .ge. run_params%numgen) call quit_de('Max number of generations already reached.')

  if (run_params%verbose .ge. 1) write(*,*) 'Restored successfully.'

end subroutine resume


end module io
