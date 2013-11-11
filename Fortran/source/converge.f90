module converge

use detypes

implicit none

logical, parameter :: meanimprovement = .true. !whether to use the mean improvement convergence criteria (currently only option)
logical, parameter :: checkpopres = .false.    !check that the population resolution is ok at all steps. 
                                               !Useful when dealing with duplicates caused *very* small population variance


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

    if (meanimprovement) then                    !FIXME: make choice of convergence criteria part of run_params
       converged = meanimproveconv(X, gen, run_params)
    else                                         !FIXME implement other convergence criteria options
       converged = .false.                       !no convergence criteria used
    end if

    if (checkpopres) call check_population_resolution(X, run_params, converged)

  end function converged

!FIXME: change meanimproveconv to be *fractional* improvement, and pass a tolerance for it from the main calling sequence
!FIXME: add weighted meanimproveconv, pass a threshold from main calling sequence

  logical function meanimproveconv(X, gen, run_params) !mean improvement in the average/best fitness of the population

    type(population), intent(in) :: X
    integer, intent(in) :: gen
    type(codeparams), intent(in) :: run_params

    real(dp), parameter :: convthresh=0.1             !if the mean improvement is less than this, population has converged
    integer, parameter :: convstep=5              !mean improvement checked using this many steps
    logical, parameter :: avgfitness=.true.       !use the average fitness of the population to check improvement. Otherwise, use the best fitness

    real(dp), save :: oldval, curval                  !the normalized best/average fitness of the population for most recent and previous generation
    real(dp) :: valdiff
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

    !make sure we're not subtracting infinity from infinity
    if (curval .gt. 1.e10) then   !FIXME: is 1e10 a good threshold value?
       valdiff = huge(1.0_dp)
    else
       valdiff = oldval-curval
    end if

    improve = eoshift(improve, shift=-1, boundary=valdiff)  !store new improvement, discard oldest improvement
    
    if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Mean improvement =', sum(improve)
    
    if (sum(improve) .lt. convthresh) then 
       if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Converged.'
       meanimproveconv = .true.
    else
       if (verbose .and. (run_params%mpirank .eq. 0)) write (*,*) '  Not converged.'
       meanimproveconv = .false.
    end if

  end function meanimproveconv


  !check if the variance of the population has gotten too small--if close to 
  !floating-point resolution, could have problems with duplicate population members
  subroutine check_population_resolution(X, run_params, converged)
    type(population), intent(in) :: X
    type(codeparams), intent(in) :: run_params
    logical, intent(inout) :: converged
    real(dp), dimension(1, run_params%D) :: avgvector
    real(dp), dimension(run_params%DE%NP, run_params%D) :: diffvectors
    real(dp) :: resolution  !a measure of the floating-point resolution for current population
    integer :: res_pt_count !max number of points allowed inside resolution
    integer i

    avgvector = reshape(sum(X%vectors, dim=1), (/1, run_params%D/))
    avgvector = avgvector/real(run_params%DE%NP, kind=dp)

    if (verbose .and. (run_params%mpirank .eq. 0)) then
       write (*,*) '  Checking population resolution...'
       write (*,*) '  Average vector:', avgvector
    end if

    do i=1,run_params%DE%NP
       diffvectors(i,:) = abs(X%vectors(i,:) - avgvector(1,:))
    end do

    !compare each dimension separately
    do i=1, run_params%D
       if (verbose .and. (run_params%mpirank .eq. 0)) write(*,*) '  Dimension:', i
       resolution = 10.0_dp*spacing(avgvector(1,i))   !gives an idea of when vectors can accidentally take on the same values

       if( any(diffvectors(:,i) .lt. resolution)) then
          if (verbose .and. (run_params%mpirank .eq. 0)) write(*,*) '    WARNING: at least one vector within allowed resolution'

          res_pt_count = run_params%DE%NP/4 !no reason for this value, but it keeps down duplicates
          if (count(diffvectors(:,i) .lt. resolution) .ge. res_pt_count) then  
             converged = .true.
             if (run_params%mpirank .eq. 0) write (*,*) 'WARNING: Points along dimension', i, &
                                                'cannot be resolved further. Ending civilization.'
          end if
       else
          if (verbose .and. (run_params%mpirank .eq. 0)) write(*,*) '    Population resolution okay.'
       end if
    end do

  end subroutine check_population_resolution


end module converge
