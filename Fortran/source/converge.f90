module converge

use detypes

implicit none

private
public converged, evidenceDone

contains


  logical function evidenceDone(Z,Zerr,tol)

    real(dp), intent(in) :: Z, tol
    real(dp), intent(inout) :: Zerr

    evidenceDone = (log(Z/(Z-Zerr)) .le. tol)

  end function evidenceDone


  logical function converged(X, gen, run_params)                  
    
    type(population), intent(in) :: X
    integer, intent(in) :: gen
    type(codeparams), intent(in) :: run_params

    logical, parameter :: meanimprovement=.true. !whether to use the mean improvement convergence criteria

    if (meanimprovement) then                    !FIXME: make choice of convergence criteria part of run_params
       converged = meanimproveconv(X, gen, run_params)
    else                                         !FIXME implement other convergence criteria options
       converged = .false.                       !no convergence criteria used
    end if

  end function converged


  logical function meanimproveconv(X, gen, run_params) !mean improvement in the average/best fitness of the population

    type(population), intent(in) :: X
    integer, intent(in) :: gen
    type(codeparams), intent(in) :: run_params

    real(dp), parameter :: convthresh=0.1             !if the mean improvement is less than this, population has converged
    integer, parameter :: convstep=5              !mean improvement checked using this many steps
    logical, parameter :: avgfitness=.true.       !use the average fitness of the population to check improvement. Otherwise, use the best fitness

    real(dp), save :: oldval, curval                  !the normalized best/average fitness of the population for most recent and previous generation
    real(dp), dimension(convstep), save :: improve    !improvement between oldval, curval over most recent steps, normalized by number of steps stored


    if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Checking convergence...'

    if (gen .eq. 1) then                                          !initialize
       improve = huge(1.0_dp)
       oldval = huge(1.0_dp)
       if (avgfitness) then
          curval = sum(X%values)/(real(size(X%values)*convstep, kind=dp))  !average of new population values, normalized by number of steps stored
       else 
          curval = minval(X%values)/real(convstep, kind=dp)                !best population value, normalized by number of steps stored
       end if
    else                                                                                     
       oldval = curval
       if (avgfitness) then
          curval = sum(X%values)/(real(size(X%values)*convstep, kind=dp))
       else 
          curval = minval(X%values)/real(convstep, kind=dp)              
       end if
    end if

    improve = eoshift(improve, shift=-1, boundary=oldval-curval)  !store new improvement, discard oldest improvement
    
    if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Mean improvement =', sum(improve)
    
    if (sum(improve) .lt. convthresh) then 
       if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Converged.'
       meanimproveconv = .true.
    else
       if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Not converged.'
       meanimproveconv = .false.
    end if

  end function meanimproveconv


end module converge
