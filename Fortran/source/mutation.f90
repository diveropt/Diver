module mutation 

use detypes
use deutils

implicit none

private
public getSubpopulation, mutate, init_FjDE, init_lambdajDE

real(dp), parameter :: tauF=0.1_dp, Fl=0.1_dp, Fu=0.9_dp !jDE control parameters as in Brest et al 2006
real(dp), parameter :: taulambda=0.1_dp !lambdajDE control parameter as above (but since 0<=lambda<=1, no need for Fu, Fl analogs)

contains


  !Determine the subpopulation of points that have the same partitioned discrete parameters
  subroutine getSubpopulation(X, Xsub, n, nsub, run_params)
    type(population), intent(in) :: X
    type(population), intent(inout) :: Xsub                      
    type(codeparams), intent(in) :: run_params
    integer, intent(in) :: n
    integer, intent(inout) :: nsub
    integer :: discrete_vals(run_params%D_discrete), subpopulationIndex, i, j
    logical :: same_subpop

    subpopulationIndex = 0
    forall(j=1:run_params%D_discrete) discrete_vals(j) = nint(X%vectors(n,run_params%discrete(j)))

    do i = 1, run_params%DE%NP
       same_subpop = .true.
       do j = 1, run_params%D_discrete
          same_subpop = same_subpop .and. ( nint(X%vectors(i,run_params%discrete(j))) .eq. discrete_vals(j) )
       enddo
       if (same_subpop) then
          subpopulationIndex = subpopulationIndex + 1
          if (i == n) nsub = subpopulationIndex
          Xsub%vectors(subpopulationIndex,:) = X%vectors(i,:)
          Xsub%values(subpopulationIndex) = X%values(i)

          if (run_params%DE%jDE) then
             Xsub%FjDE(subpopulationIndex) = X%FjDE(i)
             if (run_params%DE%lambdajDE) then
                Xsub%lambdajDE(subpopulationIndex) = X%lambdajDE(i)
             end if
          end if

       endif
    enddo    

    if (subpopulationIndex .ne. run_params%subpopNP) then
       if (run_params%verbose .ge. 1) then
          write(*,*) 'Apparent subpopulation size: ',subpopulationIndex
          write(*,*) 'Expected subpopulation size: ',run_params%subpopNP
       end if
       call quit_de('ERROR: subpopulation size does not match run_params%subpopNP!')
    endif

  end subroutine getSubpopulation


  subroutine mutate(X, V, n, run_params, trialF, triallambda) 
    type(population), intent(in) :: X                      !valid set of target vectors
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D), intent(out) :: V    !donor vector
    integer, intent(in) :: n                               !index of current vector in X
    real(dp), intent(out) :: trialF, triallambda

    if (run_params%DE%lambdajDE) then
       trialF = newF(X, n)
       triallambda = newlambda(X, n)
       V = lambdajDEmutation(X, n, run_params, trialF, triallambda)
    else if (run_params%DE%jDE) then
       trialF = newF(X, n)
       triallambda = 0.0_dp
       V = jDEmutation(X, n, run_params, trialF)
    else
       trialF = 0.0_dp
       triallambda = 0.0_dp
       V = genmutation(X, n, run_params)
    endif

  end subroutine mutate


  !general mutation strategy: V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
  function genmutation(X, n, run_params) 
    type(population), intent(in) :: X                !current generation of target vectors
    integer, intent(in) :: n                         !index of current vector
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D) :: genmutation !donor vector
    integer, dimension(2*run_params%DE%Fsize) :: r   !index of random vectors X_J, X_K
    integer ri                                       !index of (random or current) vector X_I
    integer,  dimension(1) :: rbest                  !index of best vector (array to make minloc() happy)
    integer q                                        !use to iterate over the scale factors F(q)
    real(dp), dimension(run_params%D) :: sumF        !the summed difference vector over the F's

    !assign rbest
    if (run_params%DE%lambda .gt. 0.0_dp) then
       rbest = minloc(X%values) 
    else
       rbest = (/n/)                          !don't want to restrict options for ri, r(q)
    end if

    !assign ri
    if(run_params%DE%current .or. (run_params%DE%lambda .eq. 1.0_dp)) then 
       ri = n                                !use current target vector for mutation (for lambda=1, don't need unique ri)
    else                                     !for rand/ or rand-to-best/
       do                                    !choose random ri not equal to n or rbest (if lambda>0)
          call random_int(ri, 1, run_params%subpopNP)
          if ( any( (/n, rbest(1)/) .eq. ri) ) then 
             cycle                           !ri not unique, keep trying
          else
             exit                            !finished picking ri
          end if
       end do 
    end if

    !assign unique r(q)'s from population
    do q=1, 2*run_params%DE%Fsize
       do
          call random_int(r(q), 1, run_params%subpopNP)
          if ( any( (/r(1:q-1), n, ri, rbest(1)/) .eq. r(q)) ) then
             cycle                           !continue picking new r(q)'s until unique
          else   
             exit                            !r(q) unique, begin picking r(q+1) 
          end if
       end do 
    end do

    !find the difference vector associated with the F's:
    sumF(:) = 0.0_dp
    do q=1, run_params%DE%Fsize
       sumF(:) = run_params%DE%F(q)*(X%vectors(r(2*q-1),:) - X%vectors(r(2*q),:)) + sumF(:)
    end do

    !V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
    genmutation(:) =  run_params%DE%lambda*X%vectors(rbest(1), :) + (1-run_params%DE%lambda)*X%vectors(ri,:) + sumF(:)

  end function genmutation
  

!rand/1 or rand-to-best/1 mutation using self-adaptive F parameter
  function jDEmutation(X, n, run_params, trialF)
    type(population), intent(in) :: X
    integer, intent(in) :: n 
    type(codeparams), intent(in) :: run_params
    real(dp), intent(in) :: trialF
    real(dp), dimension(run_params%D) :: jDEmutation
    integer :: r1, r2, r3
    integer, dimension(1) :: rbest


    !assign rbest if lambda is nonzero
    if (run_params%DE%lambda .gt. 0.0_dp) then
       rbest = minloc(X%values) 
    else
       rbest = (/n/)                          !don't want to restrict options for ri, r(q)
    end if

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    do
       call random_int(r1, 1, run_params%subpopNP)  !pick 1st vector from population; must not equal n
       !if (r1 .ne. n) exit
       if ( all(r1 .ne. (/n, rbest/)) ) exit
    end do
    do                                              !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, run_params%subpopNP) 
       if ( all(r2 .ne. (/n, r1, rbest/)) ) exit
    end do
    do                                              !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, run_params%subpopNP) 
       if ( all(r3 .ne. (/n, r1, r2, rbest/)) ) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    !jDEmutation = X%vectors(r1, :) + trialF*(X%vectors(r2,:) - X%vectors(r3,:))
    jDEmutation = run_params%DE%lambda*X%vectors(rbest(1), :) + (1-run_params%DE%lambda)*X%vectors(r1,:) & 
                  + trialF*(X%vectors(r2,:) - X%vectors(r3,:))

  end function jDEmutation


!rand-to-best/1 mutation using self-adaptive F parameter
  function lambdajDEmutation(X, n, run_params, trialF, triallambda)
    type(population), intent(in) :: X
    integer, intent(in) :: n 
    type(codeparams), intent(in) :: run_params
    real(dp), intent(in) :: trialF, triallambda
    real(dp), dimension(run_params%D) :: lambdajDEmutation
    integer :: r1, r2, r3
    integer, dimension(1) :: rbest

    !best vector
    rbest = minloc(X%values) 

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    do
       call random_int(r1, 1, run_params%subpopNP)  !pick 1st vector from population; must not equal n
       if ( all(r1 .ne. (/n, rbest/)) ) exit
    end do
    do                                              !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, run_params%subpopNP) 
       if ( all(r2 .ne. (/n, r1, rbest/)) ) exit
    end do
    do                                              !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, run_params%subpopNP) 
       if ( all(r3 .ne. (/n, r1, r2, rbest/)) ) exit
    end do
    !V = lambda*Xrbest + (1 - lambda)*Xr1 + F*(Xr2 - Xr3)
    lambdajDEmutation = triallambda*X%vectors(rbest(1), :) + (1 - triallambda)*X%vectors(r1,:) & 
                  + trialF*(X%vectors(r2,:) - X%vectors(r3,:))

  end function lambdajDEmutation


  function newF(X, n)          !choose a trial value for the F parameter
    type(population), intent(in) :: X
    integer, intent(in) :: n
    real(dp) :: newF
    real(dp) :: rand1, rand2

    call random_number(rand1)
    if (rand1 .lt. tauF) then     
       call random_number(rand2)
       newF = Fl + rand2*Fu    !selects a new value between 0.1 and 1.0
    else
       newF = X%FjDE(n)        !inherit F from previous generation
    endif

  end function newF


  function newlambda(X, n)         !choose a trial value for the lambda parameter
    type(population), intent(in) :: X
    integer, intent(in) :: n
    real(dp) :: newlambda
    real(dp) :: rand1, rand2

    call random_number(rand1)
    if (rand1 .lt. taulambda) then     
       call random_number(rand2)
       newlambda = rand2           !selects a new value between 0.0 and 1.0
    else
       newlambda = X%lambdajDE(n)  !inherit lambda from previous generation
    endif

  end function newlambda


  function init_FjDE(size)
    integer, intent(in) :: size
    real(dp), dimension(size) :: init_FjDE
    real(dp), dimension(size) :: rand
    
    call random_number(rand)
    init_FjDE =  Fl + rand*Fu

  end function init_FjDE


  function init_lambdajDE(size)
    integer, intent(in) :: size
    real(dp), dimension(size) :: init_lambdajDE
    real(dp), dimension(size) :: rand
    
    call random_number(rand)
    init_lambdajDE =  rand
    
  end function init_lambdajDE
  

  subroutine random_int(harvest, min, max) !choose a random integer between min and max, inclusive
    integer, intent(out) :: harvest
    integer, intent(in) :: min, max
    real(dp) :: range
    real(dp) :: r
   
    range = max - min
    call random_number(r)
    r = r*range + min
    harvest = nint(r)
  end subroutine random_int


end module mutation
