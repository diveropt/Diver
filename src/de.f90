module de

use detypes
use deutils
use init
use io
use converge
use selection
use mutation
use crossover
use post
use evidence

implicit none

#ifdef MPI
  include 'mpif.h'
#endif

public

contains


  !Main differential evolution routine.
  subroutine diver(func, lowerbounds, upperbounds, path, nDerived, discrete, partitionDiscrete,            &
                   maxciv, maxgen, NP, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE,               &
                   convthresh, convsteps, removeDuplicates, doBayesian, prior, maxNodePop, Ztolerance,     &
                   savecount, resume, outputSamples, init_population_strategy, discard_unfit_points,       &
                   max_initialisation_attempts, max_acceptable_value, seed, context, verbose)

    use iso_c_binding, only: c_ptr

    procedure(MinusLogLikeFunc)      :: func                    !the likelihood function to be minimised -- assumed to be -ln(likelihood)
    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds !boundaries of parameter space
    character(len=*), intent(in)     :: path                    !path to save samples, resume files, etc
    integer, intent(in), optional    :: nDerived                !input number of derived quantities to output
    integer, dimension(:), intent(in), optional :: discrete     !a vector listing all discrete dimensions of parameter space
    logical, intent(in), optional    :: partitionDiscrete       !split the population evenly amongst discrete parameters and evolve separately
    integer, intent(in), optional    :: maxciv                  !maximum number of civilisations
    integer, intent(in), optional    :: maxgen                  !maximum number of generations per civilisation
    integer, intent(in), optional    :: NP                      !population size (individuals per generation)
    real(dp), dimension(:), intent(in), optional :: F           !scale factor(s).  Note that this must be entered as an array.
    real(dp), intent(in), optional   :: Cr                      !crossover factor
    real(dp), intent(in), optional   :: lambda                  !mixing factor between best and rand/current
    logical, intent(in), optional    :: current                 !use current vector for mutation
    logical, intent(in), optional    :: expon                   !use exponential crossover
    integer, intent(in), optional    :: bndry                   !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional    :: jDE                     !use self-adaptive choices for rand/1/bin parameters as described in Brest et al 2006
    logical, intent(in), optional    :: lambdajDE               !use self-adaptive choices for rand-to-best/1/bin parameters; based on Brest et al 2006
    real(dp), intent(in), optional   :: convthresh              !threshold for generation-level convergence
    integer, intent(in), optional    :: convsteps               !number of steps to smooth over when checking convergence
    logical, intent(in), optional    :: removeDuplicates        !weed out duplicate vectors within a single generation
    logical, intent(in), optional    :: doBayesian              !calculate log evidence and posterior weightings
    procedure(PriorFunc), optional   :: prior                   !the prior function
    real(dp), intent(in), optional   :: maxNodePop              !population at which node is partitioned in binary space partitioning for posterior
    real(dp), intent(in), optional   :: Ztolerance              !input tolerance in log-evidence
    integer, intent(in), optional    :: savecount               !save progress every savecount generations
    logical, intent(in), optional    :: resume                  !restart from a previous run
    integer, intent(in), optional    :: init_population_strategy!initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
    logical, intent(in), optional    :: discard_unfit_points    !recalculate any trial vector whose fitness is above max_acceptable_value
    integer, intent(in), optional    :: max_initialisation_attempts !maximum number of times to try to find a valid vector for each slot in the initial population.
    real(dp), intent(in), optional   :: max_acceptable_value    !maximum fitness to accept for the initial generation if init_population_strategy > 0. Also applies to later generations if discard_unfit_points = .true.
    logical, intent(in), optional    :: outputSamples           !write samples as output
    integer, intent(in), optional    :: seed                    !base seed for random number generation; non-positive or absent means seed from the system clock
    type(c_ptr), intent(inout), optional :: context             !context pointer, used for passing info from the caller to likelihood/prior

    integer, intent(in), optional    :: verbose                 !output verbosity: 0=only error messages, 1=basic info, 2=civ-level info, 3+=population info

    type(codeparams) :: run_params                              !carries the code parameters

    type(population), target :: X, BF                           !population of target vectors, best-fit vector
    type(population) :: Xnew, Xsub                              !population for the next generation, partitioned subpopulation
    real(dp), dimension(size(lowerbounds)) :: V, U              !donor, trial vectors
    real(dp) :: trialF, trialCr, triallambda                    !adaptive F and Cr for jDE, lambda for lambdajDE

    integer :: fcall=0, accept=0                                !fcall counts function calls, accept counts acceptance rate
    integer :: totfcall = 0, totaccept = 0                      !for function calls & acceptance rates for all processes
    integer :: civ, gen, m                                      !civ, gen, n for iterating civilisation, generation, population loops
    integer :: n, nsub                                          !current member of population being evolved (same as m unless using MPI), subpop version
    integer :: civstart=1, genstart=1                           !starting values of civ, gen

    real(dp) :: Z=0., Zmsq=0., Zerr=0., Zold=0.                 !evidence
    integer :: Nsamples = 0                                     !number of statistically independent samples from posterior
    integer :: Nsamples_saved = 0                               !number of samples saved to .sam file so far
    logical :: quit = .false.                                   !flag passed from user function to indicate need to stop
    logical :: acceptable_trial_vector                          !dictates the calculation/recalculation of trial vectors

    integer :: ierror                                           !MPI error code
    logical :: mpi_already_init                                 !MPI initialization
    real(dp) :: t1, t2                                          !for timing

    call cpu_time(t1)

#ifdef MPI
    call MPI_Initialized(mpi_already_init, ierror)              !check if MPI has been initialized by the calling routine
    if (.not. mpi_already_init) call MPI_Init(ierror)
#endif

    !This will need to be changed to an input parameter if alternative covergence criteria are actually implemented
    run_params%convergence_criterion = meanimprovement

    !Assign specified or default values to run_params and print out information to the screen
    call param_assign(run_params, lowerbounds, upperbounds, nDerived=nDerived, discrete=discrete,    &
                      partitionDiscrete=partitionDiscrete, maxciv=maxciv, maxgen=maxgen, NP=NP,      &
                      F=F, Cr=Cr, lambda=lambda, current=current, expon=expon, bndry=bndry, jDE=jDE, &
                      lambdajDE=lambdajDE, convthresh=convthresh, convsteps=convsteps,               &
                      removeDuplicates=removeDuplicates, doBayesian=doBayesian,                      &
                      maxNodePop=maxNodePop, Ztolerance=Ztolerance, savecount=savecount,             &
                      outputSamples=outputSamples, init_population_strategy=init_population_strategy,&
                      discard_unfit_points=discard_unfit_points,                                     &
                      max_initialisation_attempts=max_initialisation_attempts,                       &
                      max_acceptable_value=max_acceptable_value, seed=seed, context=context,         &
                      verbose=verbose)

    if (run_params%calcZ .and. .not. present(prior)) then
       call quit_de('Error: evidence calculation requested without specifying a prior.')
    end if

    !seed the random number generator(s) from the system clock
    call init_all_random_seeds(run_params%DE%NP/run_params%mpipopchunk, run_params%mpirank, run_params%seed)

    !Allocate vector population: X is the full population, Xnew is the size of the population each process handles
    !Xsub is a subset of the population that has the same values of discrete parameters when partitionDiscrete is used.
    allocate(X%vectors(run_params%DE%NP, run_params%D))
    allocate(X%vectors_and_derived(run_params%DE%NP, run_params%D+run_params%D_derived))
    allocate(X%values(run_params%DE%NP), X%weights(run_params%DE%NP), X%multiplicities(run_params%DE%NP))
    allocate(Xsub%vectors(run_params%subpopNP, run_params%D), Xsub%values(run_params%subpopNP))
    allocate(Xnew%vectors(run_params%mpipopchunk, run_params%D))
    allocate(Xnew%vectors_and_derived(run_params%mpipopchunk, run_params%D+run_params%D_derived))
    allocate(Xnew%values(run_params%mpipopchunk))

    !for self-adaptive DE (jDE), allocate space for the populations of parameters
    if (run_params%DE%jDE) then
       allocate(X%FjDE(run_params%DE%NP))
       allocate(X%CrjDE(run_params%DE%NP))
       allocate(Xsub%FjDE(run_params%subpopNP))
       allocate(Xnew%FjDE(run_params%mpipopchunk))
       allocate(Xnew%CrjDE(run_params%mpipopchunk))

       !allocate self-adaptive lambda
       if (run_params%DE%lambdajDE) then
          allocate(X%lambdajDE(run_params%DE%NP))
          allocate(Xsub%lambdajDE(run_params%subpopNP))
          allocate(Xnew%lambdajDE(run_params%mpipopchunk))
       end if
    endif

    !Allocate best-fit containers
    allocate(BF%vectors(1, run_params%D), BF%vectors_and_derived(1, run_params%D+run_params%D_derived), BF%values(1))

    !If required, initialise the linked tree used for estimating the evidence and posterior
    if (run_params%calcZ) call iniTree(lowerbounds,upperbounds,run_params%maxNodePop)

    !Initialise internal variables
    BF%values(1) = huge(BF%values(1))

    !Resume from saved run or initialise save files for a new one.
    call io_begin(path, civstart, genstart, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, totfcall, &
                  run_params, X, BF, prior=prior, restart=resume)

    !Tidy a few things up if resuming.
    if (present(resume)) then
       if (resume) then
         genstart = genstart + 1
         !Partition totfcall into fcall in each process.  This should always be an even split, but just in case...
         if(run_params%mpirank .ne. 0) then
           fcall = floor(float(totfcall/run_params%mpiprocs))
         else
           fcall = totfcall - (run_params%mpiprocs-1)*floor(float(totfcall/run_params%mpiprocs))
         endif
       endif
    endif

    !Run a number of sequential DE optimisations, exiting either after a set number of
    !runs through or after the evidence has been calculated to a desired accuracy
    civloop: do civ = civstart, run_params%numciv

       if (run_params%calcZ) then
          !Break out if posterior/evidence is converged
          if (evidenceDone(Z,Zerr,run_params%tol)) then
            if (civ .eq. civstart) gen = genstart-1
            exit
          endif
       endif

       !Split if the quit flag has been raised
       if (quit) exit

       if (run_params%verbose .ge. 2) then
          write (*,*) '-----------------------------'
          write (*,*) 'Civilisation: ', civ
       end if

       !Internal (normal) DE loop: calculates population for each generation
       genloop: do gen = genstart, run_params%numgen

          if (run_params%verbose .ge. 2) then
             write (*,*) '  -----------------------------'
             write (*,*) '  Generation: ', gen
          end if

          if (gen .eq. 1) then

             !Initialise the convergence criterion
             call init_convergence(run_params)

             !initialise the first generation
             call initialize(X, Xnew, run_params, func, fcall, quit, accept)

             !sync quit flags
             quit = sync(quit)

             !update accept and fcall
             call update_acceptance(accept, fcall, totaccept, totfcall, run_params%verbose .ge. 2, run_params%DE%NP)

             !update best fits
             call newBFs(X,BF)

             !don't use initial generation for estimating evidence, as it biases the BSP.  Therefore no call to updateEvidence.

             !save things
             if (run_params%mpirank .eq. 0) then
                if (civ .eq. 1) call save_run_params(path, run_params)
                if (run_params%savefreq .eq. 1) then
                   call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, huge(Z), Nsamples, Nsamples_saved, totfcall, run_params)
                endif
             endif

          else

             accept = 0

             poploop: do m=1, run_params%mpipopchunk                       !evolves individual members of the population

                n = run_params%mpipopchunk*run_params%mpirank + m          !current member of the population (same as m if no MPI)

                !keep making trial vectors until one has a fitness less than max_acceptable_value, if discard_unfit_points = .true.
                do
                   if (run_params%partitionDiscrete) then
                      call getSubpopulation(X, Xsub, n, nsub, run_params)     !restrict donor pool to this member's subpopulation
                      call mutate(Xsub, V, nsub, run_params, trialF, triallambda) !create new donor vector V
                   else
                      call mutate(X, V, n, run_params, trialF, triallambda)   !create new donor vector V
                   endif

                   call gencrossover(X, V, U, n, run_params, trialCr)         !trial vectors

                   !choose next generation
                   call selector(X, Xnew, U, trialF, triallambda, trialCr, m, n, run_params, func, fcall, quit, accept, &
                                 acceptable_trial_vector)

                   if (acceptable_trial_vector) exit

                end do

                if (abs(run_params%verbose) .ge. 3) then
                   if (run_params%DE%lambdajDE) then
                      write (*,*) n, Xnew%vectors_and_derived(m, :), '->', Xnew%values(m), '|', &
                                  Xnew%lambdajDE(m), Xnew%FjDE(m), Xnew%CrjDE(m)
                   else if (run_params%DE%jDE) then
                      write (*,*) n, Xnew%vectors_and_derived(m, :), '->', Xnew%values(m), '|', Xnew%FjDE(m), Xnew%CrjDE(m)
                   else
                      write (*,*) n, Xnew%vectors_and_derived(m, :), '->', Xnew%values(m)
                   end if
                end if

             end do poploop

             !sync quit flags
             quit = sync(quit)

             !replace old generation with newly calculated one
             call replace_generation(X, Xnew, run_params, func, fcall, quit, accept, init=.false.)

             !debugging code: choose random new population members uniformly from the allowed parameter ranges
             !call initialize(X, Xnew, run_params, func, fcall, quit)

             !update accept and fcall
             call update_acceptance(accept, fcall, totaccept, totfcall, run_params%verbose .ge. 2, run_params%DE%NP)

             !Update best fits
             call newBFs(X,BF)

             !Update the evidence calculation
             if (run_params%calcZ) call updateEvidence(X, Z, Zmsq, Zerr, prior, run_params%context, Nsamples)

             !Do periodic save
             if ((mod(gen,run_params%savefreq) .eq. 0) .and. (run_params%mpirank .eq. 0)) then
                call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, huge(Z), Nsamples, Nsamples_saved, totfcall, run_params)
             endif

          endif

          if (quit .and. run_params%mpirank .eq. 0) then
             if (run_params%verbose .gt. 0) write(*,*) 'Quit requested by objective function - Diver will save and exit now.'
             call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, huge(Z), Nsamples, Nsamples_saved, totfcall, run_params, &
                          final=(mod(gen,run_params%savefreq) .eq. 0) )
          endif

          !Check generation-level convergence: if satisfied, or quit flag set, exit genloop
          if (converged(X, run_params) .or. quit) exit

       end do genloop

       genstart = 1

       !Update best fits
       call newBFs(X,BF)

       if (run_params%verbose .ge. 3) then
          write (*,*)
          write (*,*) '  ============================='
       end if
       if (run_params%verbose .ge. 2) then
          write (*,*) '  Number of generations in this civilisation: ', min(gen,run_params%numgen)
          write (*,*) '  Best final vector in this civilisation: ', roundvector(BF%vectors(1,:), run_params)
          write (*,*) '  Value at best final vector in this civilisation: ', BF%values(1)
          write (*,*) '  Cumulative function calls: ', totfcall
       end if

    enddo civloop

    !Correct civ in cases where the loop has gone through at least once
    if (civ .ne. civstart) civ = civ - 1

    if (run_params%verbose .ge. 1) then
       write (*,*) '============================='
       write (*,'(A25,I4)') ' Number of civilisations: ', min(civ,run_params%numciv)
       write (*,*) 'Best final vector: ', roundvector(BF%vectors(1,:), run_params)
       write (*,*) 'Value at best final vector: ', BF%values(1)
       if (run_params%calcZ) then
          write (*,'(A23,E13.6,A5,E13.6,A7)') ' approx. ln(Evidence): ', log(Z), ' +/- ', log(Z/(Z-Zerr)), ' (stat)'
       end if
       write (*,*) 'Total Function calls: ', totfcall
    end if

    !Polish the evidence
    if (run_params%calcZ .and. run_params%mpirank .eq. 0 .and. Nsamples_saved .gt. 0) then
      Zold = Z
      call polishEvidence(Z, Zmsq, Zerr, prior, run_params%context, Nsamples_saved, path, run_params, .true.)
      if (run_params%verbose .ge. 1) then
         write (*,'(A25,E13.6)') ' corrected ln(Evidence): ', log(Z)
         write (*,'(A25,E13.6,A6)') '                     +/- ', abs(log(Z/Zold)), ' (sys)'
         write (*,'(A25,E13.6,A7)') '                     +/- ', log(Z/(Z-Zerr)), ' (stat)'
      end if
    endif

    !Do final save operation
    if (run_params%mpirank .eq. 0 ) then
       call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, totfcall, run_params, &
                     final = ( (mod(gen,run_params%savefreq) .eq. 0) .or. (civ .eq. civstart) ) )
    end if

    !Clean up and shut down.
    if (allocated(run_params%improvements))  deallocate(run_params%improvements)
    if (allocated(X%FjDE))                   deallocate(X%FjDE)
    if (allocated(X%CrjDE))                  deallocate(X%CrjDE)
    if (allocated(X%lambdajDE))              deallocate(X%lambdajDE)
    if (allocated(Xsub%FjDE))                deallocate(Xsub%FjDE)
    if (allocated(Xsub%lambdajDE))           deallocate(Xsub%lambdajDE)
    if (allocated(Xnew%FjDE))                deallocate(Xnew%FjDE)
    if (allocated(Xnew%CrjDE))               deallocate(Xnew%CrjDE)
    if (allocated(Xnew%lambdajDE))           deallocate(Xnew%lambdajDE)
    if (allocated(run_params%DE%F))          deallocate(run_params%DE%F)
    if (allocated(run_params%discrete))      deallocate(run_params%discrete)
    if (allocated(run_params%repeat_scales)) deallocate(run_params%repeat_scales)
                                             deallocate(X%vectors, X%values, X%weights)
                                             deallocate(X%vectors_and_derived, X%multiplicities)
                                             deallocate(Xsub%vectors, Xsub%values)
                                             deallocate(Xnew%vectors, Xnew%values)
                                             deallocate(BF%vectors, BF%values, BF%vectors_and_derived)
    if (run_params%calcZ) call clearTree

    call cpu_time(t2)

#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
#endif
    if (abs(run_params%verbose) .ge. 1) then
       write (*,'(A26,I4,A2,F10.2)') ' Total seconds for process ', run_params%mpirank, ': ', t2-t1
    end if
#ifdef MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
#endif

#ifdef MPI
    if (.not. mpi_already_init) call MPI_Finalize(ierror)
#endif

  end subroutine diver


end module de
