module mutation 

use detypes

implicit none

private
public mutate

contains

  subroutine mutate(X, V, n, run_params, trialF) 
    type(population), intent(in) :: X                      !current generation of target vectors
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D), intent(out) :: V        !donor vector
    integer, intent(in) :: n                               !index of current vector
    real, intent(out) :: trialF

    if(run_params%DE%jDE) then
       trialF = newF(X, n)
       V = jDEmutation(X, n, run_params, trialF)
    else
       trialF = 0.
       V = genmutation(X, n, run_params)
    endif

  end subroutine mutate


  !general mutation strategy: V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
  function genmutation(X, n, run_params) 
    type(population), intent(in) :: X                !current generation of target vectors
    integer, intent(in) :: n                         !index of current vector
    type(codeparams), intent(in) :: run_params
    real, dimension(run_params%D) :: genmutation  !donor vector
    integer :: totF                                  !number of F scale factors
    integer, dimension(2*size(run_params%DE%F)) :: r !index of random vectors X_J, X_K
    integer ri                                       !index of (random or current) vector X_I
    integer q
    real, dimension(run_params%D) :: sumF            !the summed difference vector over the F's

    totF=size(run_params%DE%F)

    !assign unique r(q)'s from population
    do q=1, 2*totF
       do
          call random_int(r(q), 1, run_params%DE%NP)
          if ( any( (/r(1:q-1), n/) .eq. r(q)) ) then
             cycle                           !continue picking new r(q)'s until unique
          else   
             exit                            !r(q) unique, begin picking r(q+1) 
          end if
       end do 
    end do

    !assign ri
    if(run_params%DE%current) then 
       ri = n                                !use current target vector for mutation
    else
       do                                    !choose random unique ri
          call random_int(ri, 1, run_params%DE%NP)
          if ( any( (/r(:), n/) .eq. ri) ) then 
             cycle 
          else
             exit 
          end if
       end do 
    end if

    !find the difference vector associated with the F's:
    do q=1, totF
       sumF(:) = run_params%DE%F(q)*(X%vectors(r(q),:) - X%vectors(r(2*q),:))
    end do
    
    !V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
    genmutation(:) =  run_params%DE%lambda*bestvector(X, run_params) + (1-run_params%DE%lambda)*X%vectors(ri,:) + sumF(:)

  end function genmutation
  

!rand/1 mutation using self-adaptive F parameter
  function jDEmutation(X, n, run_params, trialF)
    type(population), intent(in) :: X
    integer, intent(in) :: n 
    type(codeparams), intent(in) :: run_params
    real, intent(in) :: trialF
    real, dimension(run_params%D) :: jDEmutation
    integer :: r1, r2, r3

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    call random_int(r1, 1, run_params%DE%NP)        !pick 1st vector from population
    do                                              !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, run_params%DE%NP) 
       if (r2 .ne. r1) exit
    end do
    do                                              !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, run_params%DE%NP) 
       if ((r3 .ne. r1) .and. (r3 .ne. r2)) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    jDEmutation = X%vectors(r1, :) + trialF*(X%vectors(r2,:) - X%vectors(r3,:))

  end function jDEmutation


  function newF(X, n)              !choose a trial value for the F parameter
    type(population), intent(in) :: X
    integer, intent(in) :: n
    real :: newF
    real :: rand1, rand2
    real, parameter :: tau=0.1, Fl=0.1, Fu=0.9 !control parameters as in Brest et al 2006

    call random_number(rand1)
    if (rand1 .lt. tau) then     
       call random_number(rand2)
       newF = Fl + rand2*Fu    !selects a new value between 0.1 and 1.0
    else
       newF = X%FjDE(n)        !inherit F from previous generation
    endif

  end function newF


  function bestvector(X, run_params)            !returns the best vector in the population
    type(population), intent(in) :: X
    type(codeparams), intent(in) :: run_params
    integer, dimension(1) :: locbest        !the index of the best vector
    real, dimension(run_params%D) :: bestvector !the best vector in the current population

    locbest = minloc(X%values)
    bestvector(:) = X%vectors(locbest(1), :)

  end function bestvector



  subroutine random_int(harvest, min, max) !choose a random integer between min and max, inclusive
    integer, intent(out) :: harvest
    integer, intent(in) :: min, max
    real :: range
    real :: r
   
    range = max - min
    call random_number(r)
    r = r*range + min
    harvest = nint(r)
  end subroutine random_int


end module mutation
