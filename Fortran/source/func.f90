  !function to be minimized.
  function func(X, fcall)
    real, dimension(2), intent(in) :: X
    integer, intent(inout) :: fcall
    real func

    fcall = fcall + 1
    func = (1.0 - X(1))**2 + (5.0 - X(2))**2
  end function func