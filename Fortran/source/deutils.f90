module deutils

implicit none

private
public int_to_string, quit_de

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


end module deutils
