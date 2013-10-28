! iso_c_binding interface for de::run_de
!
! Allows one to call DEvoPack using the C prototypes
!
! void runde(double (*minusloglike)(double[], const int, int&, bool&, bool), 
!            int nPar, const double lowerbounds[], const double upperbounds[], const char path[], int nDerived, 
!            int nDiscrete, const int discrete[], bool partitionDiscrete, int maxciv, int maxgen, int NP, int nF, 
!            const double F[], double Cr, double lambda, bool current, bool expon, int bndry, bool jDE, bool lambdajDE, 
!            bool removeDuplicates, bool doBayesian, double(*prior)(const double[], const int),
!            double maxNodePop, double Ztolerance, int savecount, bool resume)
!
! double minusloglike(double params[], const int param_dim, int &fcall, bool &quit, bool validvector)
!
! double prior(const double true_params[], const int true_param_dim)
!
! Loosely based on cwrapper.f90 in MultiNest v3.3 by Michele Vallisneri, 2013/09/20


module de_c_interface

implicit none

contains
	
    subroutine runde(minusloglike, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, &
                     partitionDiscrete, maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, &
                     bndry, jDE, lambdajDE, removeDuplicates, doBayesian, prior, maxNodePop, Ztolerance, savecount, &
                     resume) bind(c)

    use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, C_NULL_CHAR
    use de, only: run_de

    type(c_funptr),  intent(in), value :: minusloglike, prior
    integer(c_int),  intent(in), value :: nPar, nDerived, nDiscrete, maxciv, maxgen, NP, nF, bndry, savecount
    integer(c_int),  intent(in), target:: discrete(nDiscrete)
    logical(c_bool), intent(in), value :: partitionDiscrete, current, expon, jDE, lambdajDE, removeDuplicates, doBayesian, resume
    real(c_double),  intent(in), value :: Cr, lambda, maxNodePop, Ztolerance
    real(c_double),  intent(in)        :: lowerbounds(nPar), upperbounds(nPar), F(nF)
    character(kind=c_char,len=1), dimension(1), intent(in) :: path

    integer :: i
    integer, target :: discrete_empty(0)
    integer, pointer :: discrete_f(:)
    integer, parameter :: maxpathlen = 300
    character(len=maxpathlen) :: path_f
    interface
      real(dp) function prior_f_proto(X)
        use detypes, only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: X
      end function prior_f_proto
    end interface
    procedure(prior_f_proto), pointer :: priorPtr

    ! Fix up the string, which is null-terminated in C
    path_f = ' '
    do i = 1, maxpathlen
       if (path(i) == C_NULL_CHAR) then
          exit
       else
	  path_f(i:i) = path(i)
       end if
    end do

    ! Fix up the potential null pointer passed in instead of an illegal zero-element C array if nDiscrete = 0 
    if (nDiscrete .eq. 0) then
      discrete_f => discrete_empty
    else 
      discrete_f => discrete
    endif

    ! Choose a null pointer for the prior function if it isn't needed 
    if (doBayesian) then
      priorPtr => prior_f
    else
      priorPtr => NULL()
    endif
 
    ! Call the actual fortran differential evolution function
    call run_de(minusloglike_f, lowerbounds, upperbounds, path_f, nDerived=nDerived, discrete=discrete_f,      &
                partitionDiscrete=logical(partitionDiscrete), maxciv=maxciv, maxgen=maxgen, NP=NP, F=F, Cr=Cr, &
                lambda=lambda, current=logical(current), expon=logical(expon), bndry=bndry, jDE=logical(jDE),  &
                lambdajDE=logical(lambdajDE), removeDuplicates=logical(removeDuplicates),                      &
                doBayesian=logical(doBayesian), prior = priorPtr, maxNodePop=maxNodePop, Ztolerance=Ztolerance, &
                savecount=savecount, resume=logical(resume))
    
    contains

       ! Wrapper for the likelihood function
       real(dp) function minusloglike_f(params, fcall, quit, validvector)
       use iso_c_binding, only: c_f_procpointer
       use detypes, only: dp
  
       real(dp), intent(inout) :: params(nPar+nDerived)
       integer,  intent(inout) :: fcall 
       logical,  intent(out)   :: quit
       logical,  intent(in)    :: validvector
       logical(c_bool)         :: quit_c, validvector_c

       interface
          real(c_double) function minusloglike_proto(params, param_dim, fcall, quit, validvector)
             use iso_c_binding, only: c_int, c_double, c_bool
             implicit none
             real(c_double),  intent(inout)     :: params(*)
             integer(c_int),  intent(in), value :: param_dim
             integer(c_int),  intent(inout)     :: fcall 
             logical(c_bool), intent(out)       :: quit
             logical(c_bool), intent(in), value :: validvector
          end function minusloglike_proto
       end interface

       procedure(minusloglike_proto), pointer :: minusloglike_c

          ! Cast c_funptr minusloglike to a pointer with signature minusloglike_proto, and assign the result to minusloglike_c
          call c_f_procpointer(minusloglike,minusloglike_c)

          validvector_c = validvector  ! Do boolen type conversion default logical --> c_bool 
          minusloglike_f = minusloglike_c(params, size(params), fcall, quit_c, validvector_c)
          quit = quit_c                ! Do boolen type conversion c_bool --> default logical

       end function minusloglike_f


       ! Wrapper for the prior function
       real(dp) function prior_f(true_params)
       use iso_c_binding, only: c_f_procpointer
       use detypes, only: dp

       real(dp), intent(in) :: true_params(nPar)

       interface
          real(c_double) function prior_proto(true_params, true_param_dim)
             use iso_c_binding, only: c_double, c_int
             implicit none
             real(c_double),  intent(in)        :: true_params(*)
             integer(c_int),  intent(in), value :: true_param_dim
          end function prior_proto
       end interface

       procedure(prior_proto), pointer :: prior_c

          ! Cast c_funptr prior to a pointer with signature prior_proto, and assign the result to prior_c
          call c_f_procpointer(prior,prior_c)

          prior_f = prior_c(true_params, size(true_params))

       end function prior_f


   end subroutine


end module