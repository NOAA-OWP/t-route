!*---------------------------------------------------------------------------------------------
!*         Numerical Recipes Types (nrtype), Appendix C1.1,NR F90

!The file supplied as nrtype.f90 contains a single module named nrtype, which in turn contains
!definitions for a number of named constants (that is, PARAMETERs), and a couple of elementary
!derived data types used by the sparse matrix routines in this book. Of the named constants,
!by far the most important are those that define the KIND types of virtually all the variables
!used in this book: I4B, I2B, and I1B for integer variables, SP and DP for
!real variables (and SPC and DPC for the corresponding complex cases), and LGT for the default
!logical type
!-----------------------------------------------------------------------------------------------
module nrtype
    !* symbolic names for kind types of 4-, 2-, and 1-byte integers:
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: i2b = selected_int_kind(4)
    integer, parameter :: i1b = selected_int_kind(2)
    !* symbolic names for kind types of single- and double-precision reals:
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
    !* symbolic names for kind types of single- and double-precision complex:
    integer, parameter :: spc = kind((1.0,1.0))
    integer, parameter :: dpc = kind((1.0d0,1.0d0))
    !* symbolic name for kind type of default logical:
    integer, parameter :: lgt = kind(.true.)
    !* frequently used mathematical constants (with precision to spare):
    real(sp), parameter :: pi=3.141592653589793238462643383279502884197_sp
    real(sp), parameter :: pio2=1.57079632679489661923132169163975144209858_sp
    real(sp), parameter :: twopi=6.283185307179586476925286766559005768394_sp
    real(sp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_sp
    real(sp), parameter :: euler=0.5772156649015328606065120900824024310422_sp
    real(dp), parameter :: pi_d=3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: pio2_d=1.57079632679489661923132169163975144209858_dp
    real(dp), parameter :: twopi_d=6.283185307179586476925286766559005768394_dp
    !* derived data types for sparse matrices, single and double precision (see use in chapter b2):
    type sprs2_sp
        integer(i4b) :: n,len
        real(sp), dimension(:), pointer :: val
        integer(i4b), dimension(:), pointer :: irow
        integer(i4b), dimension(:), pointer :: jcol
    end type sprs2_sp
    type sprs2_dp
        integer(i4b) :: n,len
        real(dp), dimension(:), pointer :: val
        integer(i4b), dimension(:), pointer :: irow
        integer(i4b), dimension(:), pointer :: jcol
    end type sprs2_dp
end module nrtype

