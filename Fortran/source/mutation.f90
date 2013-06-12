module mutation 

use detypes

implicit none

private
public mutate, init_FjDE

real(dp), parameter :: tau=0.1, Fl=0.1, Fu=0.9 !jDE control parameters as in Brest et al 2006

contains

  subroutine mutate(X, V, n, run_params, trialF) 
    type(population), intent(in) :: X                      !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D), intent(out) :: V        !donor vector
    integer, intent(in) :: n                               !index of current vector
    real(dp), intent(out) :: trialF

    if(run_params%DE%jDE) then
       trialF = newF(X, n)
       V = jDEmutation(X, n, run_params, trialF)
    else
       trialF = 0.0_dp
       V = genmutation(X, n, run_params)
    endif

  end subroutine mutate


  !general mutation strategy: V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
  function genmutation(X, n, run_params) 
    type(population), intent(in) :: X                !current generation of target vectors
    integer, intent(in) :: n                         !index of current vector
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D) :: genmutation     !donor vector
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
          call random_int(ri, 1, run_params%DE%NP)
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
          call random_int(r(q), 1, run_params%DE%NP)
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
  

!rand/1 mutation using self-adaptive F parameter
  function jDEmutation(X, n, run_params, trialF)
    type(population), intent(in) :: X
    integer, intent(in) :: n 
    type(codeparams), intent(in) :: run_params
    real(dp), intent(in) :: trialF
    real(dp), dimension(run_params%D) :: jDEmutation
    integer :: r1, r2, r3

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    do
       call random_int(r1, 1, run_params%DE%NP)     !pick 1st vector from population; must not equal n
       if (r1 .ne. n) exit
    end do
    do                                              !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, run_params%DE%NP) 
       if ( all(r2 .ne. [n, r1]) ) exit
    end do
    do                                              !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, run_params%DE%NP) 
       if ( all(r3 .ne. [n, r1, r2]) ) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    jDEmutation = X%vectors(r1, :) + trialF*(X%vectors(r2,:) - X%vectors(r3,:))

  end function jDEmutation


  function newF(X, n)              !choose a trial value for the F parameter
    type(population), intent(in) :: X
    integer, intent(in) :: n
    real(dp) :: newF
    real(dp) :: rand1, rand2

    call random_number(rand1)
    if (rand1 .lt. tau) then     
       call random_number(rand2)
       newF = Fl + rand2*Fu    !selects a new value between 0.1 and 1.0
    else
       newF = X%FjDE(n)        !inherit F from previous generation
    endif

  end function newF

  function init_FjDE(run_params)
    type(codeparams), intent(in) ::run_params
    real(dp), dimension(run_params%mpipopchunk) :: init_FjDE
    real(dp), dimension(run_params%mpipopchunk) :: rand
    
    call random_number(rand)
    init_FjDE =  Fl + rand*Fu

  end function init_FjDE


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
