module deutils

use detypes

implicit none

private
public int_to_string, quit_de, roundvector

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
#ifdef USEMPI
    call MPI_Finalize(ierror)
#endif
    stop

  end subroutine quit_de


  !rounds vectors to nearest discrete values for all dimensions listed in run_params%discrete
  !all other dimensions are kept the same
  function roundvector(trialvector, run_params)
    real(dp), dimension(:), intent(in) :: trialvector
    type(codeparams), intent(in) :: run_params
    real(dp), dimension(run_params%D) :: roundvector

    roundvector = trialvector
    roundvector(run_params%discrete) = anint(roundvector(run_params%discrete))

  end function roundvector


end module deutils
