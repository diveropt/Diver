module de

use init
use io
use converge
use selection
use mutation
use crossover
use evidence

implicit none

private
public run_de

contains 


  !Main differential evolution routine.  
  subroutine run_de(func, prior, lowerbounds, upperbounds, path, nDerived, maxciv, maxgen, NP, F, Cr, lambda, current, expon, &
                    bndry, jDE, doBayesian, maxNodePop, Ztolerance, tolcount, savecount, resume)
    real, external :: func, prior 				!function to be minimized (assumed -ln[likelihood]), prior function
    real, dimension(:), intent(in) :: lowerbounds, upperbounds	!boundaries of parameter space
    character(len=*), intent(in)   :: path			!path to save samples, resume files, etc  
    integer, intent(in), optional  :: nDerived	 		!input number of derived quantities to output
    integer, intent(in), optional  :: maxciv 			!maximum number of civilisations
    integer, intent(in), optional  :: maxgen 			!maximum number of generations per civilisation
    integer, intent(in), optional  :: NP 			!population size (individuals per generation)
    real, dimension(:), intent(in), optional :: F		!scale factor(s).  Note that this must be entered as an array.
    real, intent(in), optional     :: Cr 			!crossover factor
    real, intent(in), optional     :: lambda 			!mixing factor between best and rand/current
    logical, intent(in), optional  :: current 			!use current vector for mutation
    logical, intent(in), optional  :: expon 			!use exponential crossover
    integer, intent(in), optional  :: bndry                     !boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection
    logical, intent(in), optional  :: jDE                       !use self-adaptive choices for rand/1/bin parameters as described in Brest et al 2006
    logical, intent(in), optional  :: doBayesian                !calculate log evidence and posterior weightings
    real, intent(in), optional     :: maxNodePop                !population at which node is partitioned in binary space partitioning for posterior
    real, intent(in), optional     :: Ztolerance		!input tolerance in log-evidence
    integer, intent(in), optional  :: tolcount	 		!input number of times delta ln Z < tol in a row for convergence
    integer, intent(in), optional  :: savecount			!save progress every savecount generations
    logical, intent(in), optional  :: resume			!restart from a previous run
     
    type(codeparams) :: run_params                              !carries the code parameters 

    type(population), target :: X, BF                           !population of target vectors, best-fit vector  
    real, dimension(size(lowerbounds)) :: V, U                  !donor, trial vectors
    real :: trialF, trialCr                                     !adaptive F and Cr for jDE

    integer :: fcall, accept                                    !fcall counts function calls, accept counts acceptance rate
    integer :: civ, gen, n                                      !civ, gen, n for iterating civilisation, generation, population loops

    real, dimension(size(lowerbounds)) :: avgvector, bestvector !for calculating final average and best fit
    real :: bestvalue
    real, allocatable :: bestderived(:)
    integer :: bestloc(1)

    real :: Zold, Z = 0.                                        !evidence
    integer :: Nsamples = 0                                     !number of statistically independent samples from posterior
    integer :: convcount = 0                                    !number of times delta ln Z < tol in a row so far
    
    write (*,*) '============================='
    write (*,*) ' ******** Begin DE *********'

    !seed the random number generator from the system clock
    call random_seed()
    
    !Assign specified or default values to run_params, bconstrain
    call param_assign(run_params, lowerbounds, upperbounds, nDerived=nDerived, maxciv=maxciv, maxgen=maxgen, NP=NP, F=F, &
                       Cr=Cr, lambda=lambda, current=current, expon=expon, bndry=bndry, jDE=jDE, doBayesian=doBayesian, &
                       maxNodePop=maxNodePop, Ztolerance=Ztolerance, tolcount=tolcount, savecount=savecount)

    !Resume from saved run or initialise save files for a new one
    if (present(resume)) then
       call io_begin(path, civ, gen, Z, Zold, Nsamples, convcount, run_params, restart=resume)
    else
       call io_begin(path, civ, gen, Z, Zold, Nsamples, convcount, run_params)
    endif

    !Allocate vector population
    allocate(X%vectors(run_params%DE%NP, run_params%D))
    allocate(X%derived(run_params%DE%NP, run_params%D_derived))
    allocate(X%values(run_params%DE%NP), X%weights(run_params%DE%NP), X%multiplicities(run_params%DE%NP)) 

    !for self-adaptive DE (jDE), allocate space for the populations of parameters
    if (run_params%DE%jDE) then
       allocate(X%FjDE(run_params%DE%NP))
       allocate(X%CrjDE(run_params%DE%NP))
    endif
    
    !Allocate best-fit containers
    allocate(BF%vectors(1, run_params%D), BF%derived(1, run_params%D_derived), BF%values(1), bestderived(run_params%D_derived))
    
    !If required, initialise the linked tree used for estimating the evidence and posterior
    if (run_params%calcZ) call initree(lowerbounds,upperbounds,run_params%maxNodePop)

    !Initialise internal variables
    fcall = 0
    BF%values(1) = huge(BF%values(1))

    !Run a number of sequential DE optimisations, exiting either after a set number of
    !runs through or after the evidence has been calculated to a desired accuracy
    civloop: do civ = 1, run_params%numciv

       if (verbose) write (*,*) '-----------------------------'
       if (verbose) write (*,*) 'Civilisation: ', civ

       !Initialise the first generation
       call initialize(X, run_params, lowerbounds, upperbounds, fcall, func)
       if (run_params%calcZ) call updateEvidence(X, Z, prior, Nsamples, run_params%DE%NP)        
       
       !Internal (normal) DE loop: calculates population for each generation
       genloop: do gen = 2, run_params%numgen 

          if (verbose) write (*,*) '  -----------------------------'
          if (verbose) write (*,*) '  Generation: ', gen
   
          accept = 0

          !$OMP PARALLEL DO
          poploop: do n=1, run_params%DE%NP                             !evolves one member of the population

             call mutate(X, V, n, run_params, trialF)                   !create new donor vector V
             call gencrossover(X, V, U, n, run_params, trialCr)  	!trial vectors

             !choose next generation of target vectors
             call selector(X, U, trialF, trialCr, n, lowerbounds, upperbounds, run_params, fcall, func, accept)
            
             if(run_params%DE%jDE) then 
                if (verbose) write (*,*) n, X%vectors(n, :), '->', X%values(n), '|', X%FjDE(n), X%CrjDE(n)
             else
                if (verbose) write (*,*) n, X%vectors(n, :), '->', X%values(n)
             end if

          end do poploop
          !$END OMP PARALLEL DO
 
          if (verbose) write (*,*) '  Acceptance rate: ', accept/real(run_params%DE%NP)

          if (run_params%calcZ) then
             call updateEvidence(X, Z, prior, Nsamples, run_params%DE%NP)
          endif

          !Do periodic save
          if (mod(gen,run_params%savefreq) .eq. 0) call save_all(X, path, civ, gen, Z, Zold, Nsamples, convcount, run_params)

          if (converged(X, gen)) exit                !Check generation-level convergence: if satisfied, exit loop
                                                     !PS, comment: it looks like the convergence of the evidence *requires*
                                                     !that the generation-level convergence check is done, as continuing to
                                                     !evolve a population after it has converged just results in many more 
                                                     !copies of the same point ending up in the database, which seems to
                                                     !start to introduce a bias in the evidence.
       end do genloop
       
       avgvector = sum(X%vectors, dim=1)/real(run_params%DE%NP)
       bestloc = minloc(X%values)
       bestvector = X%vectors(bestloc(1),:)
       bestderived = X%derived(bestloc(1),:)
       bestvalue = minval(X%values)
          
       !Update current best fit
       if (bestvalue .le. BF%values(1)) then
          BF%values(1) = bestvalue
          BF%vectors(1,:) = bestvector
          BF%derived(1,:) = bestderived
       endif
          
       if (verbose) write (*,*)
       if (verbose) write (*,*) '  ============================='
       if (verbose) write (*,*) '  Number of generations in this civilisation: ', min(gen,run_params%numgen)
       if (verbose) write (*,*) '  Average final vector in this civilisation: ', avgvector
       if (verbose) write (*,*) '  Value at average final vector in this civilisation: ', func(avgvector, bestderived, fcall) 
       if (verbose) write (*,*) '  Best final vector in this civilisation: ', bestvector
       if (verbose) write (*,*) '  Value at best final vector in this civilisation: ', bestvalue
       if (verbose) write (*,*) '  Cumulative function calls: ', fcall
      
       if (run_params%calcZ) then
         !Reassess the previous civilisations' contribution to the evidence
         call reassessEvidence(Z,Zold)
         !Break out if posterior/evidence is converged
         if (run_params%calcZ .and. evidenceDone(Z,Zold,run_params%tol,convcount,run_params%convcountreq)) exit
       endif

    enddo civloop

    write (*,*) '============================='
    write (*,*) 'Number of civilisations: ', min(civ,run_params%numciv)
    write (*,*) 'Best final vector: ', BF%vectors(1,:)
    write (*,*) 'Value at best final vector: ', BF%values(1)
    if (run_params%calcZ) write (*,*)   'ln(Evidence): ', log(Z)
    write (*,*) 'Total Function calls: ', fcall

    !Do final save operation
    call save_all(X, path, civ, gen, Z, Zold, Nsamples, convcount, run_params, final=.true.)

    deallocate(X%vectors, X%values, X%weights, X%derived, X%multiplicities) 
    if (allocated(X%FjDE)) deallocate(X%FjDE)
    if (allocated(X%CrjDE)) deallocate(X%CrjDE)
    if (allocated(run_params%DE%F)) deallocate(run_params%DE%F)
    deallocate(BF%vectors, BF%values, BF%derived)

  end subroutine run_de

end module de
