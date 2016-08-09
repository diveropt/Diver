module converge

use detypes

implicit none

integer, parameter :: meanimprovement = 0      !indices for convergence criteria (currently only option is mean improvement)
logical, parameter :: checkpopres = .false.    !check that the population resolution is ok at all steps. 
                                               !Useful when dealing with duplicates caused *very* small population variance

private
public init_convergence, converged, evidenceDone, meanimprovement

contains

  subroutine init_convergence(run_params)

    type(codeparams), intent(inout) :: run_params

    if (run_params%convergence_criterion == meanimprovement) then
       run_params%meanlike = huge(1.0_dp)
       run_params%improvements = 1.0_dp
    endif

  end subroutine init_convergence


  logical function evidenceDone(Z,Zerr,tol)

    real(dp), intent(in) :: Z, tol
    real(dp), intent(inout) :: Zerr

    evidenceDone = (log(Z/(Z-Zerr)) .le. tol)

  end function evidenceDone


  logical function converged(X, run_params)                     
    type(population), intent(in) :: X
    type(codeparams), intent(inout) :: run_params

    if (run_params%verbose .ge. 3) write(*,*) '  Checking convergence...'

    select case (run_params%convergence_criterion)
       case (meanimprovement) 
          converged = check_SFIM(X, run_params)
       case default                                 !TODO implement other convergence criteria options
          converged = .false.                       !no convergence criteria used
    end select

    if (converged) then
       if (run_params%verbose .ge. 3) write (*,*) '  Converged.'
    else
       if (run_params%verbose .ge. 3) write (*,*) '  Not converged.'
    end if

    if (checkpopres) call check_population_resolution(X, run_params, converged)
  end function converged


  !tracks the smoothed fractional improvement of the mean value of the population
  !at each generation, and ends the civilization when this goes below a certain threshold
  !Note that this *does not work* for test functions whose minimum is 0
  logical function check_SFIM(X, run_params) result(isConverged)

    type(population), intent(in) :: X
    type(codeparams), intent(inout) :: run_params

    real(dp) :: curval         !the average fitness of the population for the current generation
    real(dp) :: fracdiff       !the fractional difference between this and the previous generation's mean value
    real(dp) :: sfim           !the smoothed fractional improvement in the mean
    real(dp), parameter :: inf_threshold = 0.001*huge(1.0_dp) !to make sure there are no problems with inifinity

    !set the current value to the average of the current population
    curval = sum(X%values)/(real(size(X%values), kind=dp))
    !curval = minval(X%values)            !best population value

    !make sure we don't have problems with infinity
    if (curval .gt. inf_threshold) then     
       run_params%meanlike = inf_threshold
       fracdiff = 1.0_dp
    else
       fracdiff = 1.0_dp - curval/run_params%meanlike  !the fractional improvement between this generation and last generation
       run_params%meanlike = curval
    end if

    run_params%improvements = eoshift(run_params%improvements, shift=-1, boundary=fracdiff)   !store new improvement, discard oldest improvement
    sfim = sum(run_params%improvements)/real(run_params%convsteps, kind=dp)                   !average over the generations stored
    
    if (run_params%verbose .ge. 3) write (*,*) '  Smoothed fractional improvement of the mean =', sfim
    
    !compare to threshold value
    isConverged = (sfim .lt. run_params%convthresh)

  end function check_SFIM



  !check if the variance of the population has gotten too small--if close to 
  !floating-point resolution, could have problems with duplicate population members
  !This is not really useful except for smooth functions
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

    if (run_params%verbose .ge. 3) then
       write (*,*) '  Checking population resolution...'
       write (*,*) '  Average vector:', avgvector
    end if

    do i=1,run_params%DE%NP
       diffvectors(i,:) = abs(X%vectors(i,:) - avgvector(1,:))
    end do

    !compare each dimension separately
    do i=1, run_params%D
       if (run_params%verbose .ge. 3) write(*,*) '  Dimension:', i
       resolution = 10.0_dp*spacing(avgvector(1,i))   !gives an idea of when vectors can accidentally take on the same values

       if( any(diffvectors(:,i) .lt. resolution)) then
          if (run_params%verbose .ge. 3) write(*,*) '    WARNING: at least one vector within allowed resolution'

          res_pt_count = run_params%DE%NP/4 !no reason for this value, but it keeps down duplicates
          if (count(diffvectors(:,i) .lt. resolution) .ge. res_pt_count) then  
             converged = .true.
             if (run_params%verbose .ge. 3) write (*,*) 'WARNING: Points along dimension', i, &
                                                'cannot be resolved further. Ending civilization.'
          end if
       else
          if (run_params%verbose .ge. 3) write(*,*) '    Population resolution okay.'
       end if
    end do

  end subroutine check_population_resolution


end module converge
