! iso_c_binding interface for DEvoPack de::run_de
! Pat Scott 2013/10/16
! Loosely based on cwrapper.f90 in MultiNest v3.3 by Michele Vallisneri, 2013/09/20
!
! Allows one to call DEvoPack using the C prototypes
!
! void runDE(double (*minusloglike)(double*, int, bool, bool), double(*prior)(double*), double lowerbounds[], double upperbounds[], 
!            char path[], int nDerived, int discrete[], bool partitionDiscrete, int maxciv, int maxgen, int NP, double F[],
!            double Cr, double lambda, bool current, bool expon, int bndry, bool jDE, bool lambdajDE,
!            bool removeDuplicates, bool doBayesian, double maxNodePop, double Ztolerance, int savecount, bool resume)
!
! double minusloglike(double* params, int fcall, bool quit, bool validvector)
!
! double prior(double* true_params)
!


module de_c_interface

contains
	
    subroutine runDE(minusloglike, prior, lowerbounds, upperbounds, path, nDerived, discrete, partitionDiscrete, &
                    maxciv, maxgen, NP, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE,                    &
                    removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount, resume) bind(c)

    use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, C_NULL_CHAR
    use de, only: run_de
    implicit none

    type(c_funptr),  intent(in), value :: minusloglike, prior
    integer(c_int),  intent(in), value :: nDerived, maxciv, maxgen, NP, bndry, savecount
    integer(c_int),  intent(in)        :: discrete(:)
    logical(c_bool), intent(in), value :: partitionDiscrete, current, expon, jDE, lambdajDE, removeDuplicates, doBayesian, resume
    real(c_double),  intent(in), value :: Cr, lambda, maxNodePop, Ztolerance
    real(c_double),  intent(in)        :: lowerbounds(:), upperbounds(:), F(:)
    character(kind=c_char,len=1), dimension(1), intent(in) :: path

    integer :: i
    integer, parameter :: maxpathlen = 300
    character(len=maxpathlen) :: fortran_path

    fortran_path = ' '
    do i = 1, maxpathlen
       if (path(i) == C_NULL_CHAR) then
          exit
       else
	  fortran_path(i) = path(i)
       end if
    end do

    call run_de(minusloglike_f, prior_f, lowerbounds, upperbounds, path, nDerived=nDerived, discrete=discrete,    &
                partitionDiscrete=logical(partitionDiscrete), maxciv=maxciv, maxgen=maxgen, NP=NP, F=F, Cr=Cr,    &
                lambda=lambda, current=logical(current), expon=logical(expon), bndry=bndry, jDE=logical(jDE),     &
                lambdajDE=logical(lambdajDE), removeDuplicates=logical(removeDuplicates),                         &
                doBayesian=logical(doBayesian), maxNodePop=maxNodePop, Ztolerance=Ztolerance, savecount=savecount,&
                resume=logical(resume))
	
    contains

       real(c_double) function minusloglike_f(params, fcall, quit, validvector)
       use iso_c_binding, only: c_double, c_f_procpointer
	  
       implicit none

       integer          :: fcall
       double precision :: params(:)
       logical          :: quit, validvector

       interface
          real(c_double) function minusloglike_proto(params, fcall, quit, validvector)
             use iso_c_binding, only: c_int, c_double, c_bool
             implicit none
             real(c_double),  intent(inout) :: params(:)
             integer(c_int),  intent(inout) :: fcall 
             logical(c_bool), intent(out)   :: quit
             logical(c_bool), intent(in)    :: validvector
          end function minusloglike_proto
       end interface

          procedure(minusloglike_proto), pointer :: minusloglike_c
          call c_f_procpointer(minusloglike,minusloglike_c)

          minusloglike_f = minusloglike_c(params, fcall, quit, validvector)

       end function minusloglike_f


       real(c_double) function prior_f(realparams)
       use iso_c_binding, only: c_double, c_f_procpointer
	  
       implicit none

       double precision :: realparams(:)

       interface
          real(c_double) function prior_proto(realparams)
             use iso_c_binding, only: c_double
             implicit none
             real(c_double),  intent(inout) :: realparams(:)
          end function prior_proto
       end interface

          procedure(prior_proto), pointer :: prior_c
          call c_f_procpointer(prior,prior_c)

          prior_f = prior_c(realparams)

       end function prior_f


   end subroutine


end module
