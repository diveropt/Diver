module deutils

use detypes

implicit none

#ifdef MPI
  include 'mpif.h'
#endif

private
public int_to_string, quit_de, quit_all_processes, roundvector, newBFs

contains

  function int_to_string(int)
    integer, intent(in) :: int
    character (len=30) :: string, int_to_string

    write (string, *) int
    string = adjustl(string)
    string = trim(string)

    int_to_string = string

  end function int_to_string


  subroutine quit_de(error_message)
    character(LEN=*), intent(in), optional :: error_message
    integer ierror

    if (present(error_message)) write (*,*) error_message
#ifdef MPI
    call MPI_Finalize(ierror)
#endif
    stop

  end subroutine quit_de


  subroutine quit_all_processes(error_message)
    character(LEN=*), intent(in), optional :: error_message
    integer ierror, errorcode

    if (present(error_message)) write (*,*) error_message
#ifdef MPI
    errorcode = 1
    call MPI_Abort(MPI_COMM_WORLD, errorcode, ierror)
#endif
    stop
  end subroutine quit_all_processes


  !rounds vectors to nearest discrete values for all dimensions listed in run_params%discrete
  !all other dimensions are kept the same
  function roundvector(trialvector, run_params)
    real(dp), dimension(:), intent(in) :: trialvector
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D) :: roundvector

    roundvector = trialvector
    roundvector(run_params%discrete) = anint(roundvector(run_params%discrete))

  end function roundvector


  !Updates current best fit
  subroutine newBFs(X,BF)
    type(population), intent(in) :: X     !population of vectors
    type(population), intent(inout) :: BF !best-fit vector
    integer :: bestloc(1)
    real(dp) :: bestvalue

    bestloc = minloc(X%values)
    bestvalue = X%values(bestloc(1))
    if (bestvalue .le. BF%values(1)) then
      BF%values(1) = bestvalue
      BF%vectors(1,:) = X%vectors(bestloc(1),:)
      BF%vectors_and_derived(1,:) = X%vectors_and_derived(bestloc(1),:)
    endif

  end subroutine newBFs

end module deutils
