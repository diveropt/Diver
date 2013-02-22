module selection

use detypes

implicit none

private
public selector, replace_generation

contains 

  subroutine selector(X, Xtemp, U, trialF, trialCr, n, lowerbounds, upperbounds, &
                       run_params, fcall, func, accept)

    type(population), intent(in) :: X
    type(population), intent(inout) :: Xtemp
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
       Xtemp%vectors(n,:) = trialvector 
       Xtemp%derived(n,:) = trialderived
       Xtemp%values(n) = trialvalue
       if (run_params%DE%jDE) then            !in jDE, also keep F and Cr
          Xtemp%FjDE(n) = trialF
          Xtemp%CrjDE(n) = trialCr
       end if
       accept = accept + 1
    else
       Xtemp%vectors(n,:) = X%vectors(n,:) 
       Xtemp%derived(n,:) = X%derived(n,:)
       Xtemp%values(n) = X%values(n)
       if (run_params%DE%jDE) then
          Xtemp%FjDE(n) = X%FjDE(n)
          Xtemp%CrjDE(n) = X%CrjDE(n)
       end if
    end if

  end subroutine selector


  !replaces old generation (X) with the new generation (Xtemp) calculated during population loop
  subroutine replace_generation(X, Xtemp, run_params)
    type(population), intent(inout) :: X              !old population, will be replaced
    type(population), intent(inout) :: Xtemp          !recently calculated population
    type(codeparams), intent(in) :: run_params
    integer :: k, kmatch                              !indices for vector compared, possible matching vector
    
    !weed out any duplicate vectors to maintain population diversity. One duplicate will be kept and the other 
    !will revert to its value in the previous generation
    if (run_params%DE%removeDuplicates) then
       checkpop: do k=1, run_params%DE%NP-1                                        !look for matches in higher-indexed Xtemp%values  
          if ( any(Xtemp%values(k) .eq. Xtemp%values(k+1:run_params%DE%NP)) ) then !there is at least one possible match

             findmatch: do kmatch=k+1, run_params%DE%NP                            !loop over subpopulation to find the matching vector(s)
                if ( all(Xtemp%vectors(k,:) .eq. Xtemp%vectors(kmatch,:)) ) then   !we've found a duplicate vector
                   if (verbose) write (*,*) '  Duplicate vectors:', k, kmatch

                   !Now, compare their counterparts in the previous generation to decide which vector will be kept, which will be reverted
                   picksurvivor: if (Xtemp%values(k) .eq. X%values(k)) then        !vector at k was inherited, so keep it & revert kmatch
                      if (verbose) write (*,*) '  Reverting vector...'                     
                      call replace_vector(Xtemp, X, run_params, kmatch)
                      if (verbose) write (*,*) kmatch, Xtemp%vectors(kmatch, :), '->', Xtemp%values(kmatch)

                   else if (Xtemp%values(kmatch) .eq. X%values(kmatch)) then       !vector at kmatch was inherited. Keep it
                      if (verbose) write (*,*) '  Reverting vector ', k
                      call replace_vector(Xtemp, X, run_params, k)
                      if (verbose) write (*,*) k, Xtemp%vectors(k, :), '->', Xtemp%values(k)

                   else if (X%values(k) .lt. X%values(kmatch)) then                !kmatch improved more, so keep it
                      if (verbose) write (*,*) '  Reverting vector ', k
                      call replace_vector(Xtemp, X, run_params, k)
                      if (verbose) write (*,*) k, Xtemp%vectors(k, :), '->', Xtemp%values(k)

                   else                                                            !k improved more, so keep it
                      if (verbose) write (*,*) '  Reverting vector ', kmatch
                      call replace_vector(Xtemp, X, run_params, kmatch)
                      if (verbose) write (*,*) kmatch, Xtemp%vectors(kmatch, :), '->', Xtemp%values(kmatch)

                   end if picksurvivor
                end if

             end do findmatch

          end if
       end do checkpop
    end if

    !replace old population members with new ones
    X%vectors = Xtemp%vectors
    X%values = Xtemp%values
    X%derived = Xtemp%derived
    if (run_params%DE%jDE) then
       X%FjDE = Xtemp%FjDE
       X%CrjDE = Xtemp%CrjDE
    end if

  end subroutine replace_generation


!replace a vector in Xtemp by its counterpart in the previous generation (X)
  subroutine replace_vector(Xtemp, X, run_params, l)
    type(population), intent(inout) :: Xtemp
    type(population), intent(in) :: X
    type(codeparams), intent(in) :: run_params
    integer, intent(in) :: l  !index of vector to replace
    
    Xtemp%vectors(l,:) = X%vectors(l,:)
    Xtemp%values(l) = X%values(l)
    Xtemp%derived(l,:) = X%derived(l,:)
    if (run_params%DE%jDE) then
       Xtemp%FjDE(l) = X%FjDE(l)
       Xtemp%CrjDE(l) = X%CrjDE(l)
    end if
    
  end subroutine replace_vector

end module selection
