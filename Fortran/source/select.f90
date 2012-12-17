module select

use detypes

implicit none

private
public selection

contains 

  subroutine selection(X, U, trialF, trialCr, n, lowerbounds, upperbounds, &
                       run_params, fcall, func, accept)

    type(population), intent(inout) :: X
    integer, intent(inout) :: fcall, accept
    real, dimension(:), intent(in) :: U
    real, intent(in) :: trialF, trialCr
    integer, intent(in) :: n 
    real, dimension(:), intent(in) :: lowerbounds, upperbounds 
    type(codeparams), intent(in) :: run_params   
    real, external :: func

    real :: trialvalue
    real, dimension(size(U)) :: trialvector  
    real, dimension(size(X%derived(1,:))) :: trialderived

    trialderived = 0.

    if (any(U(:) .gt. upperbounds) .or. any(U(:) .lt. lowerbounds)) then 
       !trial vector exceeds parameter space bounds: apply boundary constraints
       select case (run_params%DE%bconstrain)
          case (1)                           !'brick wall'
             trialvalue = huge(1.0)
             trialvector(:) = X%vectors(n,:)
          case (2)                           !randomly re-initialize
             call random_number(trialvector(:))
             trialvector(:) = trialvector(:)*(upperbounds - lowerbounds) + lowerbounds
             trialvalue = func(trialvector, trialderived, fcall)
          case (3)                           !reflection
             trialvector = U
             where (U .gt. upperbounds) trialvector = upperbounds - (U - upperbounds)
             where (U .lt. lowerbounds) trialvector = lowerbounds + (lowerbounds - U)
             trialvalue = func(trialvector, trialderived, fcall)
          case default                       !boundary constraints not enforced
             trialvector = U                
             trialvalue = func(U(:), trialderived, fcall)  
          end select
    else                                     !trial vector is within parameter space bounds, so use it
       trialvector = U                    
       trialvalue = func(U(:), trialderived, fcall)  
    end if

    !when the trial vector is at least as good as the current member  
    !of the population, use the trial vector for the next generation
    if (trialvalue .le. X%values(n)) then
       X%vectors(n,:) = trialvector 
       X%derived(n,:) = trialderived
       X%values(n) = trialvalue
       if (run_params%DE%jDE) then            !in jDE, also keep F and Cr
          X%FjDE(n) = trialF
          X%CrjDE(n) = trialCr
       end if
       accept = accept + 1
    end if

  end subroutine selection

end module select
