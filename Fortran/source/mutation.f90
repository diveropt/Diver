module mutation 

use detypes

implicit none

private
public genmutation

contains

  !general mutation strategy: V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
  function genmutation(X, n, params) 
    type(population), intent(in) :: X         !current generation of target vectors
    integer, intent(in) :: n                  !index of current vector
    type(deparams), intent(in) :: params
    real, dimension(params%D) :: genmutation  !donor vector
    integer :: totF                           !number of F scale factors
    integer, dimension(2*size(params%F)) :: r !index of random vectors X_J, X_K
    integer ri                                !index of (random or current) vector X_I
    integer q
    real, dimension(params%D) :: sumF         !the summed difference vector over the F's

    totF=size(params%F)

    !assign unique r(q)'s from population
    do q=1, 2*totF
       do
          call random_int(r(q), 1, params%NP)
          if( any(r(1:q-1) .eq. r(q)) ) then 
             cycle 
          else
             exit 
          end if
       end do 
    end do

    !assign ri
    if(params%current) then 
       ri = n                             !use current target vector for mutation
    else
       do                                 !choose random unique ri
          call random_int(ri, 1, params%NP)
          if( any(r(:) .eq. ri) ) then 
             cycle 
          else
             exit 
          end if
       end do 
    end if

    !find the difference vector associated with the F's:
    do q=1, totF
       sumF(:) = params%F(q)*(X%vectors(r(q),:) - X%vectors(r(2*q),:))
    end do
    
    !V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
    genmutation(:) =  params%lambda*bestvector(X, params) + (1-params%lambda)*X%vectors(ri,:) + sumF(:)

  end function genmutation
  

  function bestvector(X, params)
    type(population), intent(in) :: X
    type(deparams), intent(in) :: params
    integer, dimension(1) :: locbest        !the index of the best vector
    real, dimension(params%D) :: bestvector !the best vector in the current population

    locbest = minloc(X%values)
    bestvector(:) = X%vectors(locbest(1), :)

  end function bestvector



  function rand1mutation(X, n, params)
    type(population), intent(in) :: X   !current generation of target vectors
    integer, intent(in) :: n            !index of current vector !FIXME n is not used in this function
    type(deparams), intent(in) :: params
    real, dimension(params%D) :: rand1mutation !donor vector
    integer :: r1, r2, r3

    !set each D-dimensional donor vector in V by picking 3 separate random vectors from X
    call random_int(r1, 1, params%NP)    !pick 1st vector from population
    do                            !pick 2nd vector; ensure vectors are distinct
       call random_int(r2, 1, params%NP) 
       if (r2 .ne. r1) exit
    end do
    do                            !pick 3rd vector; ensure vectors are distinct
       call random_int(r3, 1, params%NP) 
       if ((r3 .ne. r1) .and. (r3 .ne. r2)) exit
    end do
    !V = Xr1 + F*(Xr2 - Xr3)
    rand1mutation(:) = X%vectors(r1,:) + params%F*(X%vectors(r2,:) - X%vectors(r3,:))
  end function rand1mutation


!FIXME add other mutation strategies


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
