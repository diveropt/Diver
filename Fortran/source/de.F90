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

#ifdef USEMPI
  include 'mpif.h'
#endif

private
public run_de

contains 


  !Main differential evolution routine.  
  subroutine run_de(func, lowerbounds, upperbounds, path, nDerived, discrete, partitionDiscrete, &
                    maxciv, maxgen, NP, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE,           &
                    removeDuplicates, doBayesian, prior, maxNodePop, Ztolerance, savecount, resume)

    real(dp), dimension(:), intent(in) :: lowerbounds, upperbounds !boundaries of parameter space
    character(len=*), intent(in)   :: path			!path to save samples, resume files, etc  
    integer, intent(in), optional  :: nDerived	 		!input number of derived quantities to output
    integer, dimension(:), intent(in), optional :: discrete     !a vector listing all discrete dimensions of parameter space
    logical, intent(in), optional  :: partitionDiscrete         !split the population evenly amongst discrete parameters and evolve separately
    integer, intent(in), optional  :: maxciv 			!maximum number of civilisations
    integer, intent(in), optional  :: maxgen 			!maximum number of generations per civilisation
    integer, intent(in), optional  :: NP 			!population size (individuals per generation)
    real(dp), dimension(:), intent(in), optional :: F		!scale factor(s).  Note that this must be entered as an array.
    real(dp), intent(in), optional :: Cr 			!crossover factor
    real(dp), intent(in), optional :: lambda 			!mixing factor between best and rand/current
    logical, intent(in), optional  :: current 			!use current vector for mutation
    logical, intent(in), optional  :: expon 			!use exponential crossover
    integer, intent(in), optional  :: bndry                     !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional  :: jDE                       !use self-adaptive choices for rand/1/bin parameters as described in Brest et al 2006
    logical, intent(in), optional  :: lambdajDE                 !use self-adaptive choices for rand-to-best/1/bin parameters; based on Brest et al 2006
    logical, intent(in), optional  :: removeDuplicates          !weed out duplicate vectors within a single generation
    logical, intent(in), optional  :: doBayesian                !calculate log evidence and posterior weightings
    real(dp), intent(in), optional :: maxNodePop                !population at which node is partitioned in binary space partitioning for posterior
    real(dp), intent(in), optional :: Ztolerance		!input tolerance in log-evidence
    integer, intent(in), optional  :: savecount			!save progress every savecount generations
    logical, intent(in), optional  :: resume			!restart from a previous run
     
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

    real(dp), dimension(size(lowerbounds)) :: bestvector        !for calculating final best fit
    real(dp) :: bestvalue
    real(dp), allocatable :: bestvecderived(:)
    integer :: bestloc(1)

    real(dp) :: Z=0., Zmsq=0., Zerr = 0.                        !evidence
    integer :: Nsamples = 0                                     !number of statistically independent samples from posterior
    integer :: Nsamples_saved = 0                               !number of samples saved to .sam file so far
    logical :: quit						!flag passed from user function to indicate need to stop 

    integer :: ierror		                                !MPI error code
    real(dp) :: t1, t2                                          !for timing

    interface
    !the likelihood function to be minimised -- assumed to be -ln(likelihood)
       real(dp) function func(params, fcall, quit, validvector)
          use detypes
          implicit none
          real(dp), dimension(:), intent(inout) :: params
          integer, intent(inout) :: fcall 
          logical, intent(out) :: quit
          logical, intent(in) :: validvector
       end function func
    end interface
	
    optional :: prior
    interface
    !the prior function
       real(dp) function prior(X)
          use detypes
          implicit none
          real(dp), dimension(:), intent(in) :: X
       end function prior
    end interface


    call cpu_time(t1)

#ifdef USEMPI
    call MPI_Init(ierror)
#endif

    !Assign specified or default values to run_params and print out information to the screen
    call param_assign(run_params, lowerbounds, upperbounds, nDerived=nDerived, discrete=discrete, &
                      partitionDiscrete=partitionDiscrete, maxciv=maxciv, maxgen=maxgen, NP=NP, &
                      F=F, Cr=Cr, lambda=lambda, current=current, expon=expon, bndry=bndry, jDE=jDE, & 
                      lambdajDE=lambdajDE, removeDuplicates=removeDuplicates, doBayesian=doBayesian, &
                      maxNodePop=maxNodePop, Ztolerance=Ztolerance, savecount=savecount)

    !seed the random number generator(s) from the system clock
    call init_all_random_seeds(run_params%DE%NP/run_params%mpipopchunk, run_params%mpirank)

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
    allocate(BF%vectors(1, run_params%D), BF%vectors_and_derived(1, run_params%D+run_params%D_derived))
    allocate(BF%values(1), bestvecderived(run_params%D+run_params%D_derived))
    
    !If required, initialise the linked tree used for estimating the evidence and posterior
    if (run_params%calcZ) call iniTree(lowerbounds,upperbounds,run_params%maxNodePop)

    !Initialise internal variables
    BF%values(1) = huge(BF%values(1))

    !Resume from saved run or initialise save files for a new one.
    !FIXME The switching I've implmented here is ridiculous.  Surely there is a better way, some option forwarding or similar??
    if (present(resume)) then 
       if (present(prior)) then
         call io_begin(path, civstart, genstart, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior=prior, &
          restart=resume)
       else
         call io_begin(path, civstart, genstart, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, X, BF, &
          restart=resume)
       endif
       if (resume) genstart = genstart + 1
    else 
       if (present(prior)) then
         call io_begin(path, civstart, genstart, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior=prior)
       else
         call io_begin(path, civstart, genstart, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, X, BF)
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

       if (run_params%mpirank .eq. 0) then
          if (verbose) write (*,*) '-----------------------------'
          if (verbose) write (*,*) 'Civilisation: ', civ
       end if
       
       !Internal (normal) DE loop: calculates population for each generation
       genloop: do gen = genstart, run_params%numgen 
          if (run_params%mpirank .eq. 0) then
             if (verbose) write (*,*) '  -----------------------------'
             if (verbose) write (*,*) '  Generation: ', gen
          end if
     
          if (gen .eq. 1) then 

            !Initialise the first generation
            call initialize(X, Xnew, run_params, fcall, func, quit)
            !Don't use initial generation for estimating evidence, as it biases the BSP
            if ((civ .eq. 1) .and. (run_params%mpirank .eq. 0)) call save_run_params(path, run_params)
            
          else
             
             accept = 0

             poploop: do m=1, run_params%mpipopchunk                       !evolves individual members of the population

                n = run_params%mpipopchunk*run_params%mpirank + m          !current member of the population (same as m if no MPI)

                if (run_params%partitionDiscrete) then
                   call getSubpopulation(X, Xsub, n, nsub, run_params)     !restrict donor pool to this member's subpopulation
                   call mutate(Xsub, V, nsub, run_params, trialF, triallambda) !create new donor vector V
                else
                   call mutate(X, V, n, run_params, trialF, triallambda)   !create new donor vector V
                endif
               
                call gencrossover(X, V, U, n, run_params, trialCr)         !trial vectors 
    
                !choose next generation          
                call selector(X, Xnew, U, trialF, triallambda, trialCr, m, n, run_params, fcall, func, quit, accept)

                if (verbose) then
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

             !replace old generation with newly calculated one
             call replace_generation(X, Xnew, run_params, func, fcall, quit, accept, init=.false.)

             !debugging code: choose random new population members uniformly from the allowed parameter ranges
             !call initialize(X, Xnew, run_params, fcall, func, quit)

#ifdef USEMPI
             call MPI_Allreduce(accept, totaccept, 1, MPI_integer, MPI_sum, MPI_COMM_WORLD, ierror)
#else
             totaccept = accept
#endif

             if (verbose .and. run_params%mpirank .eq. 0) write (*,*) '  Acceptance rate: ', totaccept/real(run_params%DE%NP)

             !Update the evidence calculation
             if (run_params%calcZ) call updateEvidence(X, Z, Zmsq, Zerr, prior, run_params%priorVolume, Nsamples)

             !Do periodic save
             if ((mod(gen,run_params%savefreq) .eq. 0) .and. (run_params%mpirank .eq. 0)) then
                call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params)
             endif

          endif

          if (quit) then
             write(*,*) 'Quit requested by objective function - saving and exiting.'
             if (run_params%mpirank .eq. 0) then 
                call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, &
                             final=(mod(gen,run_params%savefreq) .eq. 0) )
             end if
             call quit_de()
          endif

          if (converged(X, gen, run_params)) exit    !Check generation-level convergence: if satisfied, exit genloop
                                                     !PS, comment: it looks like the convergence of the evidence *requires*
                                                     !that the generation-level convergence check is done, as continuing to
                                                     !evolve a population after it has converged just results in many more 
                                                     !copies of the same point ending up in the database, which seems to
                                                     !start to introduce a bias in the evidence.

       end do genloop
       
       genstart = 1

       !best fit for this civilization
       bestloc = minloc(X%values)
       bestvector = X%vectors(bestloc(1),:)
       bestvecderived = X%vectors_and_derived(bestloc(1),:)
       bestvalue = minval(X%values)
          
       !Update current best fit
       if (bestvalue .le. BF%values(1)) then
          BF%values(1) = bestvalue
          BF%vectors(1,:) = bestvector
          BF%vectors_and_derived(1,:) = bestvecderived
       endif

       !get the total number of function calls
#ifdef USEMPI
       call MPI_Allreduce(fcall, totfcall, 1, MPI_integer, MPI_sum, MPI_COMM_WORLD, ierror)
#else
       totfcall = fcall
#endif
 
       if (run_params%mpirank .eq. 0) then
          !bestvecderived(:run_params%D) = avgvector
          if (verbose) write (*,*)
          if (verbose) write (*,*) '  ============================='
          if (verbose) write (*,*) '  Number of generations in this civilisation: ', min(gen,run_params%numgen)
          if (verbose) write (*,*) '  Best final vector in this civilisation: ', roundvector(bestvector, run_params)
          if (verbose) write (*,*) '  Value at best final vector in this civilisation: ', bestvalue
          if (verbose) write (*,*) '  Cumulative function calls: ', totfcall
       end if

    enddo civloop

    !Correct civ in cases where the loop has gone through at least once
    if (civ .ne. civstart) civ = civ - 1

    if (run_params%mpirank .eq. 0) then
       write (*,*) '============================='
       write (*,'(A25,I4)') 'Number of civilisations: ', min(civ,run_params%numciv)
       write (*,*) 'Best final vector: ', roundvector(BF%vectors(1,:), run_params)
       write (*,*) 'Value at best final vector: ', BF%values(1)
       if (run_params%calcZ) write (*,*)   'ln(Evidence): ', log(Z), ' +/- ', log(Z/(Z-Zerr))
       write (*,*) 'Total Function calls: ', totfcall
    end if

    !Polish the evidence
    if (run_params%calcZ .and. run_params%mpirank .eq. 0 .and. Nsamples_saved .gt. 0) then
      call polishEvidence(Z, Zmsq, Zerr, prior, Nsamples_saved, path, run_params, .true.)     
      write (*,*)   'corrected ln(Evidence): ', log(Z), ' +/- ', log(Z/(Z-Zerr))
    endif

    !Do final save operation
    if (run_params%mpirank .eq. 0 ) then
       call save_all(X, BF, path, civ, gen, Z, Zmsq, Zerr, Nsamples, Nsamples_saved, fcall, run_params, &
                     final = ( (mod(gen,run_params%savefreq) .eq. 0) .or. (civ .eq. civstart) ) )
    end if

    deallocate(X%vectors, X%values, X%weights, X%vectors_and_derived, X%multiplicities)
    deallocate(Xsub%vectors, Xsub%values)
    deallocate(Xnew%vectors, Xnew%values)
    if (allocated(X%FjDE)) deallocate(X%FjDE)
    if (allocated(X%CrjDE)) deallocate(X%CrjDE)
    if (allocated(X%lambdajDE)) deallocate(X%lambdajDE)
    if (allocated(Xsub%FjDE)) deallocate(Xsub%FjDE)
    if (allocated(Xsub%lambdajDE)) deallocate(Xsub%lambdajDE)
    if (allocated(Xnew%FjDE)) deallocate(Xnew%FjDE)
    if (allocated(Xnew%CrjDE)) deallocate(Xnew%CrjDE)
    if (allocated(Xnew%lambdajDE)) deallocate(Xnew%lambdajDE)
    if (allocated(run_params%DE%F)) deallocate(run_params%DE%F)
    if (allocated(run_params%discrete)) deallocate(run_params%discrete)
    if (allocated(run_params%repeat_scales)) deallocate(run_params%repeat_scales)
    deallocate(BF%vectors, BF%values, BF%vectors_and_derived)
    if (run_params%calcZ) call clearTree

#ifdef USEMPI
    call MPI_Finalize(ierror)
#endif

    call cpu_time(t2)

    write (*,'(A23,I4,A2,F7.2)') 'Total time for process ', run_params%mpirank, ': ', t2-t1

  end subroutine run_de


end module de
