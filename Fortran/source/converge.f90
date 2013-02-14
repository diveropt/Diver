module converge

use detypes

implicit none

private
public converged, evidenceDone

contains


  logical function evidenceDone(Z,Zerr,tol)

    real, intent(in) :: Z, tol
    real, intent(inout) :: Zerr

    evidenceDone = (log(Z/(Z-Zerr)) .le. tol)

  end function evidenceDone


  logical function converged(X, gen)                  
    
    type(population), intent(in) :: X
    integer, intent(in) :: gen

    logical, parameter :: meanimprovement=.true. !whether to use the mean improvement convergence criteria

    if (meanimprovement) then
       converged = meanimproveconv(X, gen)
    else                                         !FIXME implement other convergence criteria options
       converged = .false.                       !no convergence criteria used
    end if

  end function converged



  !FIXME implement other convergence criteria options



  logical function meanimproveconv(X, gen)       !mean improvement in the average/best fitness of the population

    type(population), intent(in) :: X
    integer, intent(in) :: gen

    real, parameter :: convthresh=0.1             !if the mean improvement is less than this, population has converged
    integer, parameter :: convstep=5              !mean improvement checked using this many steps
    logical, parameter :: avgfitness=.true.       !use the average fitness of the population to check improvement. Otherwise, use the best fitness

    real, save :: oldval, curval                  !the normalized best/average fitness of the population for most recent and previous generation
    real, dimension(convstep), save :: improve    !improvement between oldval, curval over most recent steps, normalized by number of steps stored


    if (verbose) write (*,*) '  Checking convergence...'

    if (gen .eq. 1) then                                          !initialize
       improve = huge(1.0)
       oldval = huge(1.0)
       if (avgfitness) then
          curval = sum(X%values)/(real(size(X%values)*convstep))  !average of new population values, normalized by number of steps stored
       else 
          curval = minval(X%values)/real(convstep)                !best population value, normalized by number of steps stored
       end if
    else                                                                                     
       oldval = curval
       if (avgfitness) then
          curval = sum(X%values)/(real(size(X%values)*convstep))
       else 
          curval = minval(X%values)/real(convstep)              
       end if
    end if

    improve = eoshift(improve, shift=-1, boundary=oldval-curval)  !store new improvement, discard oldest improvement
    
    if (verbose) write (*,*) '  Mean improvement =', sum(improve)
    
    if (sum(improve) .lt. convthresh) then 
       if (verbose) write (*,*) '  Converged.'
       meanimproveconv = .true.
    else
       if (verbose) write (*,*) '  Not converged.'
       meanimproveconv = .false.
    end if

  end function meanimproveconv


end module converge
