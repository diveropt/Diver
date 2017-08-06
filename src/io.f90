module io

use detypes
use deutils
use evidence
use converge

implicit none

private
public io_begin, save_all, save_run_params, resume

integer :: rawlun, samlun, devolun, rparamlun
real(dp), parameter :: Ztolscale = 100., Ftolscale = 100., Bndtolscale = 100.

contains


subroutine io_begin(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior, restart)

  character(len=*), intent(in) :: path
  integer, intent(inout) :: civ, gen, Nsamples, Nsamples_saved, fcall
  real(dp), intent(inout) :: Z, Zmsq, Zerr, Zold
  type(codeparams), intent(inout) :: run_params
  logical, intent(in), optional :: restart
  integer :: filestatus
  type(population), intent(inout) :: X, BF
  procedure(PriorFunc), optional :: prior

  logical           :: restart_

  restart_ = .false.
  if (present(restart)) restart_ = restart

  if (restart_) then
    if (present(prior)) then
      call resume(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior=prior)
    else
      call resume(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF)
    endif
  else if (run_params%mpirank .eq. 0) then
    if (run_params%outputSamples) then
      !Create .raw and .sam files.  .rparam and .devo files are created only at first save, so their existence indicates whether
      !resuming is allowed or not.
      if (run_params%verbose .ge. 1) write(*,*) 'Creating Diver output files at '//trim(path)//'.*'
      open(newunit=rawlun, file=trim(path)//'.raw', iostat=filestatus, action='WRITE', status='REPLACE')
      if ( (run_params%D_derived .ne. 0) .or. (size(run_params%discrete) .ne. 0) ) then
         open(newunit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='REPLACE')
      end if
      if (filestatus .ne. 0) call quit_all_processes(' Error creating output files. Quitting...')
      close(rawlun)
      if (run_params%D_derived .ne. 0) close(samlun)
    endif
  endif

end subroutine io_begin


subroutine save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, final)

  type(population), intent(in) :: X, BF
  character(len=*), intent(in) :: path
  integer, intent(inout) :: Nsamples_saved
  integer, intent(in) :: civ, gen, Nsamples, fcall
  real(dp), intent(in) :: Z, Zmsq, Zerr, Zold
  type(codeparams), intent(in) :: run_params
  logical, intent(in), optional :: final

  logical         :: final_

  final_ = .false.
  if (present(final)) final_ = final

  if (.not. final_) then
    Nsamples_saved = Nsamples_saved + run_params%DE%NP
    call save_samples(X, path, civ, gen, run_params)
  endif
  call save_state(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

end subroutine save_all


subroutine save_samples(X, path, civ, gen, run_params)

  type(population), intent(in) :: X
  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen
  type(codeparams), intent(in) :: run_params
  integer :: filestatus, i
  character(len=28) :: formatstring_raw
  character(len=28) :: formatstring_sam

  if (.not. run_params%outputSamples) return

  open(newunit=rawlun, file=trim(path)//'.raw', iostat=filestatus, action='WRITE', status='OLD', POSITION='APPEND')
  if (filestatus .ne. 0) call quit_all_processes(' Error opening raw file.  Quitting...')
  write(formatstring_raw,'(A18,I4,A6)') '(2E20.9,2x,2I6,2x,', run_params%D, 'E20.9)'
  do i = 1, size(X%weights)
    write(rawlun,formatstring_raw) X%multiplicities(i), X%values(i), civ, gen, X%vectors(i,:)
  enddo
  close(rawlun)

  if ( (run_params%D_derived .ne. 0) .or. (size(run_params%discrete) .ne. 0) ) then
    open(newunit=samlun, file=trim(path)//'.sam', iostat=filestatus, action='WRITE', status='OLD', POSITION='APPEND')
    if (filestatus .ne. 0) call quit_all_processes(' Error opening sam file.  Quitting...')
    write(formatstring_sam,'(A18,I4,A6)') '(2E20.9,2x,2I6,2x,', run_params%D+run_params%D_derived, 'E20.9)'
    do i = 1, size(X%weights)
      write(samlun,formatstring_sam) X%multiplicities(i), X%values(i), civ, gen, X%vectors_and_derived(i,:)
    enddo
    close(samlun)
  endif

end subroutine save_samples


subroutine save_run_params(path, run_params)

  character(len=*), intent(in) :: path
  type(codeparams), intent(in) :: run_params
  integer :: filestatus
  logical :: exists
  character(len=31) :: formatstring

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
  write(rparamlun,'(2I6)')    run_params%numciv, run_params%numgen      !maximum number of civilizations, generations
  write(rparamlun,'(E20.9)')  run_params%convthresh                     !threshold for gen-level convergence
  write(rparamlun,'(I6)')     run_params%convsteps                      !number of steps to smooth over when checking convergence
  write(rparamlun,'(E20.9)')  run_params%tol                            !tolerance in log-evidence
  write(rparamlun,'(E20.9)')  run_params%maxNodePop                     !maximum population to allow in a cell before partitioning it
  write(rparamlun,'(L1)')     run_params%calcZ                          !calculate evidence or not
  write(rparamlun,'(L1)')     run_params%outputSamples                  !output parameter samples or not
  write(rparamlun,'(I6)')     run_params%savefreq                       !frequency with which to save progress
  write(rparamlun,'(L1)')     run_params%DE%removeDuplicates            !true: remove duplicate vectors in a generation
  write(rparamlun,'(I6)')     run_params%verbose                        !amount of output to print to the screen
  write(rparamlun,'(I6)')     run_params%convergence_criterion          !indicates which convergence criterion has been selected (see convergence.f90 for codes)

  close(rparamlun)

end subroutine save_run_params


subroutine save_state(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

  character(len=*), intent(in) :: path
  integer, intent(in) :: civ, gen, Nsamples, Nsamples_saved, fcall
  real(dp), intent(in) :: Z, Zmsq, Zerr, Zold
  type(codeparams), intent(in) :: run_params
  integer :: filestatus
  logical :: exists
  character(len=31) :: formatstring
  type(population), intent(in) :: X, BF

  !Save restart info
  inquire(file=trim(path)//'.devo',exist=exists)
  if (exists) then
     open(newunit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='OLD')
  else
     open(newunit=devolun, file=trim(path)//'.devo', iostat=filestatus, action='WRITE', status='REPLACE')
  endif
  if (filestatus .ne. 0) call quit_all_processes(' Error opening devo file.  Quitting...')

  write(devolun,'(2I10)')     civ, gen                                  !current civilisation, generation
  write(devolun,'(4E20.9)')   Z, Zmsq, Zerr, Zold                       !current evidence, mean square, stat. uncertainty, approx Z if Z=corrected
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


subroutine read_state(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

  real(dp), intent(out) :: Z, Zmsq, Zerr, Zold
  integer, intent(out) :: civ, gen, Nsamples, Nsamples_saved, fcall
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
  read(rparamlun,'(2I6)')    run_params%numciv, run_params%numgen       !maximum number of civilizations, generations
  read(rparamlun,'(E20.9)')  run_params%convthresh                      !threshold for gen-level convergence
  read(rparamlun,'(I6)')     run_params%convsteps                       !number of steps to smooth over when checking convergence
  read(rparamlun,'(E20.9)')  run_params%tol                             !tolerance in log-evidence
  read(rparamlun,'(E20.9)')  run_params%maxNodePop                      !maximum population to allow in a cell before partitioning it
  read(rparamlun,'(L1)')     run_params%calcZ                           !calculate evidence or not
  read(rparamlun,'(L1)')     run_params%outputSamples                   !output parameter samples or not
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

  read(devolun,'(2I10)')     civ, gen                                   !current civilisation, generation
  read(devolun,'(4E20.9)')   Z, Zmsq, Zerr, Zold                        !current evidence, mean square, stat. uncertainty, approx Z if Z=corrected
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
subroutine resume(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior)

  character(len=*), intent(in) :: path
  integer, intent(inout) :: civ, gen, Nsamples, Nsamples_saved, fcall
  integer :: reclen, filestatus, i, j, passoverlen
  real(dp), intent(inout) :: Z, Zmsq, Zerr, Zold
  procedure(PriorFunc), optional :: prior
  real(dp) :: Z_new, Zmsq_new, Zerr_new, Z_3, Zmsq_3, Zerr_3
  character(len=31) :: formatstring
  character(len=1) :: LF
  logical :: require_Z_match = .true.
  type(codeparams), intent(inout) :: run_params
  type(codeparams) :: run_params_restored
  type(population), intent(inout) :: X, BF
  type(population) :: Y

  if (run_params%verbose .ge. 1) write(*,*) 'Restoring from previous run...'

  !Read the run state
  run_params_restored%DE%NP = run_params%DE%NP
  call read_state(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params_restored, X, BF)
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

  if (run_params%calcZ) then
    if (.not. run_params_restored%outputSamples) call quit_de('Error: cannot resume in Bayesian mode if samples were not output.')
    if (.not. run_params_restored%calcZ) call quit_de('Error: cannot resume in Bayesian mode from non-Bayesian run.')
    if (.not. present(prior)) call quit_de('Error: evidence calculation requested without specifying a prior.')
    if (any(abs(run_params%upperbounds-run_params_restored%upperbounds)/run_params%upperbounds .ge. &
                Bndtolscale*epsilon(run_params%upperbounds))                                   .or. &
        any(abs(run_params%lowerbounds-run_params_restored%lowerbounds)/run_params%lowerbounds .ge. &
                Bndtolscale*epsilon(run_params%lowerbounds)) )                                 then
       call quit_de('Error: cannot resume in Bayesian mode with a modified prior box.')
    end if
    if ( ((run_params%convthresh .ne. run_params_restored%convthresh) .or. &
         (run_params%convsteps .ne. run_params_restored%convsteps))  .and. (run_params%verbose .ge. 1) ) then
       write(*,*) 'WARNING: changing the generation-level convergence parameters mid-run may make evidence inaccurate.'
    end if
    if (run_params%MaxNodePop .ne. run_params_restored%MaxNodePop) then
       call quit_de('Error: you cannot change MaxNodePopulation mid-run!')
    end if
    if (.not. (run_params%DE%jDE .or. run_params_restored%DE%jDE)) then
      if ((run_params%DE%Fsize .ne. run_params_restored%DE%Fsize) .and. (run_params%verbose .ge. 1)) then
        write(*,*) 'WARNING: changing the number of F parameters mid-run may make evidence inaccurate.'
      elseif (run_params%DE%Fsize .ne. 0) then
        if ( any(abs(run_params%DE%F-run_params_restored%DE%F)/run_params%DE%F .ge. Ftolscale*epsilon(run_params%DE%F))) then
          if (run_params%verbose .ge. 1) write(*,*) 'WARNING: changing F values mid-run may make evidence inaccurate.'
        endif
      endif
    endif
    if ( any ( (/ run_params%DE%lambda     .ne.   run_params_restored%DE%lambda,     &
                  run_params%DE%current    .neqv. run_params_restored%DE%current,    &
                  run_params%DE%Cr         .ne.   run_params_restored%DE%Cr,         &
                  run_params%DE%expon      .neqv. run_params_restored%DE%expon,      &
                  run_params%DE%bconstrain .ne.   run_params_restored%DE%bconstrain, &
                  run_params%DE%jDE        .neqv. run_params_restored%DE%jDE,        &
                  run_params%DE%lambdajDE  .neqv. run_params_restored%DE%lambdajDE   /) ) ) then
      if (run_params%verbose .ge. 1) write(*,*) 'WARNING: changing DE algorithm mid-run may make evidence inaccurate!'
    endif
  endif

  if (mod(Nsamples_saved,run_params%DE%NP) .ne. 0) then
    call quit_de('Error: resumed run does not contain only full generations - file likely corrupted.')
  endif
  if (Nsamples .ne. Nsamples_saved) then
     if (run_params%verbose .ge. 1 .and. run_params%calcZ) then
        write(*,*) 'WARNING: running evidence from restored chain will differ to saved value, '
        write(*,*) 'as not all points used for the previous error calculation were saved.'
     end if
    require_Z_match = .false.
  endif

  allocate(Y%vectors(run_params%DE%NP, run_params%D))
  allocate(Y%vectors_and_derived(run_params%DE%NP, run_params%D+run_params%D_derived))
  allocate(Y%values(run_params%DE%NP), Y%weights(run_params%DE%NP), Y%multiplicities(run_params%DE%NP))


  !Rebuild the binary spanning tree by reading the points in by generation and sending them climbing

  !Organise the read/write format
  write(formatstring,'(A18,I4,A9)') '(2E20.9,2x,2I6,2x,', run_params%D, 'E20.9,A1)'
  reclen = 57 + 20*run_params%D

  !open the chain file
  open(newunit=rawlun, file=trim(path)//'.raw', &
   iostat=filestatus, status='OLD', access='DIRECT', action='READ', recl=reclen, form='FORMATTED')
  if (filestatus .ne. 0) call quit_all_processes(' Error opening .raw file. Quitting...')

  Z_new = 0.0_dp
  Zmsq_new = 0.0_dp
  Zerr_new = 0.0_dp
  Nsamples = 0

  !loop over the generations in the raw file to recreate the BSP tree
  do i = 1, Nsamples_saved/run_params%DE%NP
    !read in a generation
    do j = 1, run_params%DE%NP
      !read in a point
      read(rawlun,formatstring,rec=(i-1)*run_params%DE%NP+j) Y%multiplicities(j), Y%values(j), civ, gen, &
       Y%vectors(j,:), LF
    enddo
    !Update the evidence calculation
    if (run_params%calcZ) call updateEvidence(Y, Z_new, Zmsq_new, Zerr_new, prior, run_params%context, Nsamples)
  enddo

  close(rawlun)

  !Make sure we haven't already passed the number civs or gens
  if (civ .gt. run_params%numciv) call quit_de('Max number of civilisations already reached.')
  if (civ .eq. run_params%numciv .and. gen .ge. run_params%numgen) call quit_de('Max number of generations already reached.')

  !Check agreement of the evidence things with what was read in from devo file
  if (run_params%calcZ .and. require_Z_match) then
    if (any(abs((/(Z_new-Z)/Z, (Zmsq_new-Zmsq)/Zmsq, (Zerr_new - Zerr)/Zerr/)) .gt. Ztolscale*epsilon(Z))) then
      call polishEvidence(Z_3, Zmsq_3, Zerr_3, prior, run_params%context, Nsamples_saved, path, run_params, .false.)
      if (any(abs((/(Z_3-Z)/Z, (Zmsq_3-Zmsq)/Zmsq, (Zerr_3 - Zerr)/Zerr/)) .gt. Ztolscale*epsilon(Z))) then
        write(*,*) ' Evidence variables look fishy...'
        write(*,'(A24,3F16.5)') '  From devo file: ', log(Z), log(Zmsq), log(Zerr)
        write(*,'(A24,3F16.5)') '  From samples: ', log(Z_new), log(Zmsq_new), log(Zerr_new)
        write(*,'(A24,3F16.5)') '  From polished samples: ',log(Z_3), log(Zmsq_3), log(Zerr_3)
        call quit_de(' Error: evidence variables in devo file do not exactly match sample file:')
        !Z = Z_new; Zmsq = Zmsq_new; Zerr = Zerr_new
      else
        if (run_params_restored%tol .le. run_params%tol .and. run_params_restored%numciv .ge. run_params%numciv) then
          call quit_de(' This run was already completed.  Quitting...')
        else
          write(*,*) ' This run was converged already, but I will try to do a tighter job...'
          Z = Z_new; Zmsq = Zmsq_new; Zerr = Zerr_new
        endif
      endif
    endif
  endif

  if (run_params%verbose .ge. 1) write(*,*) 'Restored successfully.'

end subroutine resume


end module io
