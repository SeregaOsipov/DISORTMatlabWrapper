#include "fintrf.h"
!c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!c Author: Sergey Osipov, KAUST, May 2015
!c matlab wrapper for disort routine
!c for question email to Serega.Osipov@gmail.com
!c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!C     Gateway routine
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    !C     Declarations
    implicit none

    !C     mexFunction arguments:
    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs

    !C     Function declarations:
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxCreateNumericArray
    integer*4 mxClassIDFromClassName
    integer mxIsNumeric
    mwPointer mxGetM, mxGetN

    !C     Pointers to input/output mxArrays:r

    !C     Array information:
    mwPointer mrows, ncols
    mwSize size

    !C     Arguments for mxCreateNumericArray
    integer*4 classid
    integer*4 complexflag
    mwSize ndim
    mwSize dims(3)

    !C     Arguments for computational routine:
    !C     .. scalar arguments ..
    CHARACTER HEADER*127
    INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
    LOGICAL  USRANG, USRTAU, ONLYFL, PRNT(5), &
            PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
    real(kind = 4) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
            PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0
    real(kind = 4), parameter :: EARTH_RADIUS = 6371.0

    !     .. array arguments ..
    real(kind = 4), dimension(:), allocatable :: DTAUC, PHI, SSALB, TEMPER, UMU, UTAU
    real(kind = 4), dimension(:, :), allocatable :: PMOM
    real(kind = 4), dimension(:, :, :), allocatable :: RHOQ, RHOU
    real(kind = 4), dimension(:), allocatable :: EMUST, BEMST
    real(kind = 4), dimension(:, :), allocatable :: RHO_ACCURATE
    real(kind = 4), dimension(:), allocatable :: RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED
    real(kind = 4), dimension(:, :, :), allocatable :: UU
    real(kind = 4), dimension(:), allocatable :: H_LYR

    !C     temproray storage for MATLAB logical items
    integer*1 logicalItem
    character*200 shortLine  ! for debugging it is nice to have compact prints
    character*1000 line
    character*15000 longLine
    integer k
    integer mexPrintf
    integer x

    !C-----------------------------------------------------------------------
    !C     Check for proper number of arguments.
    if(nrhs .ne. 34) then
        call mexErrMsgIdAndTxt ('MATLAB:disortMatlab:nInput', 'exactly 34 input parameters required.')
    endif

    !C     Validate inputs
    !C     TODO:// check that input is single and not is double!
    !C     Check that the input is a number.
    !C      if(mxIsNumeric(prhs(1)) .eq. 0) then
    !C         call mexErrMsgIdAndTxt ('MATLAB:disortMatlab:NonNumeric', 'Input must be a number.')
    !C      endif

    !C     Get the size of the input array.
    mrows = mxGetM(prhs(1))
    ncols = mxGetN(prhs(1))
    size = mrows * ncols

    !C     Create Fortran array from the input argument.
    call mxCopyPtrToInteger8(mxGetPr(prhs(1)), NLYR, 1)
!    write(shortLine, *) '1, NLYR = ', NLYR
!    k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger8(mxGetPr(prhs(4)), NMOM, 1)
!   write(shortLine, *) '4, NMOM = ', NMOM
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal8(mxGetPr(prhs(7)), WVNMLO, 1)
!   write(shortLine, *) '7, WVNMLO = ', WVNMLO
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal8(mxGetPr(prhs(8)), WVNMHI, 1)
!   write(shortLine, *) '8, WVNMHI = ', WVNMHI
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(9)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, USRTAU)
!   write(shortLine, *) '9, USRTAU = ', USRTAU
!   k = mexPrintf(shortLine // achar(13))

    if (USRTAU) then  ! NTAU could be input/output variable
        call mxCopyPtrToInteger8(mxGetPr(prhs(10)), NTAU, 1)
    else
        NTAU = NLYR + 1
    end if
!   write(shortLine, *) '10, NTAU = ', NTAU
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger8(mxGetPr(prhs(12)), NSTR, 1)
!   write(shortLine, *) '12, NSTR = ', NSTR
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(13)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, USRANG)
!   write(shortLine, *) '13, USRANG = ', USRANG
!   k = mexPrintf(shortLine // achar(13))

    if (USRANG) then  ! NUMU and UMU could be input/output variable
        call mxCopyPtrToInteger8(mxGetPr(prhs(14)), NUMU, 1)
    else
        NUMU = NSTR
    endif
!   write(shortLine, *) '14, NUMU = ', NUMU
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger8(mxGetPr(prhs(16)), NPHI, 1)
!   write(shortLine, *) '16, NPHI = ', NPHI
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger8(mxGetPr(prhs(18)), IBCND, 1)
!   write(shortLine, *) '18, IBCND = ', IBCND
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(19)), FBEAM, 1)
!   write(shortLine, *) '19, FBEAM = ', FBEAM
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(20)), UMU0, 1)
!   write(shortLine, *) '20, UMU0 = ', UMU0
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(21)), PHI0, 1)
!   write(shortLine, *) '21, PHI0 = ', PHI0
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(22)), FISOT, 1)
!   write(shortLine, *) '22, FISOT = ', FISOT
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(23)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, LAMBER)
!   write(shortLine, *) '23, LAMBER = ', LAMBER
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(24)), ALBEDO, 1)
!   write(shortLine, *) '24, ALBEDO = ', ALBEDO
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(25)), BTEMP, 1)
!   write(shortLine, *) '25, BTEMP = ', BTEMP
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(26)), TTEMP, 1)
!   write(shortLine, *) '26, TTEMP = ', TTEMP
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(27)), TEMIS, 1)
!   write(shortLine, *) '27, TEMIS = ', TEMIS
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(28)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, PLANK)
!   write(shortLine, *) '28, PLANK = ', PLANK
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(29)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, ONLYFL)
!   write(shortLine, *) '29, ONLYFL = ', ONLYFL
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(30)), ACCUR, 1)
!   write(shortLine, *) '30, ACCUR = ', ACCUR
!   k = mexPrintf(shortLine // achar(13))

!    call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )
    allocate(DTAUC(NLYR), SSALB(NLYR), PMOM(0:NMOM, NLYR), &
            TEMPER(0:NLYR), UTAU(NTAU), UMU(NUMU), PHI(NPHI), H_LYR(0:NLYR))
    allocate(RHOQ(NSTR / 2, 0:NSTR / 2, 0:(NSTR - 1)), RHOU(NUMU, 0:NSTR / 2, 0:(NSTR - 1)), &
            EMUST(NUMU), BEMST(NSTR / 2), RHO_ACCURATE(NUMU, NPHI))
    allocate(RFLDIR(NTAU), RFLDN(NTAU), FLUP(NTAU), DFDT(NTAU), UAVG(NTAU), &
            ALBMED(NUMU), TRNMED(NUMU), UU(NUMU, NTAU, NPHI))
    DTAUC = 0.0; SSALB = 0.0; PMOM = 0.0; TEMPER = 0.0; UTAU = 0.0; UMU = 0.0; PHI = 0.0;
    H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
    RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; DFDT = 0.0; UAVG = 0.0; UU = 0.0;
    ALBMED = 0.0; TRNMED = 0.0;

    mrows = mxGetM(prhs(2))
    ncols = mxGetN(prhs(2))
    size = mrows * ncols
    call mxCopyPtrToReal4(mxGetPr(prhs(2)), DTAUC, size)
!   write(line, *) '2, DTAUC = ', DTAUC
!   k = mexPrintf(line // achar(13))

    call mxCopyPtrToReal4(mxGetPr(prhs(3)), SSALB, size)
!   write(line, *) '3, SSALB = ', SSALB
!   k = mexPrintf(line // achar(13))

    mrows = mxGetM(prhs(6))
    ncols = mxGetN(prhs(6))
    size = mrows * ncols
    call mxCopyPtrToReal4(mxGetPr(prhs(6)), TEMPER, size)
!   write(line, *) '6, TEMPER = ', TEMPER
!   k = mexPrintf(line // achar(13))

    mrows = mxGetM(prhs(5))
    ncols = mxGetN(prhs(5))
    size = mrows * ncols
    call mxCopyPtrToReal4(mxGetPr(prhs(5)), PMOM, size)
!   write(longLine, *) '5, PMOM = ', PMOM
!   k = mexPrintf(longLine // achar(13))

    if (USRTAU) then  ! NTAU could be input/output variable
        mrows = mxGetM(prhs(11))
        ncols = mxGetN(prhs(11))
        size = mrows * ncols
        call mxCopyPtrToReal4(mxGetPr(prhs(11)), UTAU, size)  ! UTAU could be input/output variable
!       write(line, *) '11, UTAU = ', UTAU
!       k = mexPrintf(line // achar(13))
    end if

    if (USRANG) then  !C     NUMU and UMU could be input/output variable
        mrows = mxGetM(prhs(15))
        ncols = mxGetN(prhs(15))
        size = mrows * ncols
        call mxCopyPtrToReal4(mxGetPr(prhs(15)), UMU, size)
!        write(line, *) '15, UMU = ', UMU
!       k = mexPrintf(line // achar(13))
    endif

    mrows = mxGetM(prhs(17))
    ncols = mxGetN(prhs(17))
    size = mrows * ncols
    call mxCopyPtrToReal4(mxGetPr(prhs(17)), PHI, size)
!   write(line, *) '17, PHI = ', PHI
!   k = mexPrintf(line // achar(13))

    HEADER = 'dummy header'  ! TODO: header is not read properly but I don't need it
    !C      call mxCopyPtrToCharacter(mxGetPr(prhs(32)),HEADER,127)
    !C      k=mexPrintf(HEADER//achar(13))

    ! TODO: hardcoded boolean values just yet
    PRNT(1) = .false.  ! input
    PRNT(2) = .false.  ! fluxes
    PRNT(3) = .false.  ! intensities
    PRNT(4) = .false.  ! transmissivity & albedo
    PRNT(5) = .false.  ! phase function moments

!    PRNT(1) = .true.
!    PRNT(2) = .true.
!    PRNT(3) = .true.
!    PRNT(4) = .false.
!    PRNT(5) = .false.

    call mxCopyPtrToInteger1(mxGetPr(prhs(33)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, DELTAMPLUS)
!   write(shortLine, *) '33, DELTAMPLUS = ', DELTAMPLUS
!   k = mexPrintf(shortLine // achar(13))

    call mxCopyPtrToInteger1(mxGetPr(prhs(34)), logicalItem, 1)
    call mxCopyInteger1ToLogical(logicalItem, DO_PSEUDO_SPHERE)
!   write(shortLine, *) '34, DO_PSEUDO_SPHERE = ', DO_PSEUDO_SPHERE
!   k = mexPrintf(shortLine // achar(13))

!    k=mexPrintf("in Wrapper: before DISORT"//achar(13))

    CALL DISORT(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU, &
            USRANG, USRTAU, IBCND, ONLYFL, PRNT, &
            PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE, &
            DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI, &
            UTAU, UMU0, PHI0, UMU, PHI, FBEAM, &
            FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, &
            EARTH_RADIUS, H_LYR, &
            RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST, &
            ACCUR, HEADER, &
            RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, &
            ALBMED, TRNMED)

!    k=mexPrintf("in Wrapper: after disort"//achar(13))

    !C     Prepare MATLAB output variables
    !C     Create matrix for the return argument.
    size = NTAU * 1
    classid = mxClassIDFromClassName('single')  ! double
    complexflag = 0
    ndim = 2
    dims(1) = NTAU
    dims(2) = 1

    plhs(1) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(RFLDIR, mxGetPr(plhs(1)), size)

    plhs(2) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(RFLDN, mxGetPr(plhs(2)), size)

    plhs(3) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(FLUP, mxGetPr(plhs(3)), size)

    plhs(4) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(DFDT, mxGetPr(plhs(4)), size)

    plhs(5) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(UAVG, mxGetPr(plhs(5)), size)

    size = NUMU * NTAU * NPHI
    ndim = 3
    dims(1) = NUMU
    dims(2) = NTAU
    dims(3) = NPHI

    plhs(6) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(UU, mxGetPr(plhs(6)), size)

    size = NUMU * 1
    ndim = 2
    dims(1) = NUMU
    dims(2) = 1
    plhs(7) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(ALBMED, mxGetPr(plhs(7)), size)

    plhs(8) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(TRNMED, mxGetPr(plhs(8)), size)

    !C     independent of USRANG always return UMU
    !C     it will the copy of the provided UMU or from the DISORT
    plhs(9) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(UMU, mxGetPr(plhs(9)), size)

    !C     independent of USRTAU always return UTAU
    !C     it will the copy of the provided UTAU or from the DISORT
    size = NTAU * 1
    ndim = 2
    dims(1) = NTAU
    dims(2) = 1
    plhs(10) = mxCreateNumericArray(ndim, dims, classid, complexflag)
    call mxCopyReal4ToPtr(UTAU, mxGetPr(plhs(10)), size)

!    k=mexPrintf("in Wrapper: before deallocate"//achar(13))
!    call deallocate_disort_allocatable_arrays()
    deallocate( DTAUC, SSALB, PMOM, TEMPER, UTAU, UMU, PHI, H_LYR )
    deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )
    deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED, UU )

    return
end

!C     this routine converts a matlab logical to fortran logical, matrices are not implemented yet
subroutine mxCopyInteger1ToLogical(logicaldata, fortran)
    implicit none
    integer*1, intent(in) :: logicaldata
    logical, intent(out) :: fortran
    fortran = (logicaldata /= 0)
    return
end subroutine

!subroutine allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)
!implicit none
!integer,intent(in) :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU

!allocate( DTAUC( NLYR ), SSALB( NLYR ), PMOM( 0:NMOM, NLYR ), &
!          TEMPER( 0:NLYR ), UTAU( NTAU ), UMU( NUMU ), PHI( NPHI ), H_LYR( 0:NLYR ) )
!allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
!          EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )
!allocate( RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ), DFDT( NTAU ), UAVG( NTAU ),&
!          ALBMED( NUMU ), TRNMED( NUMU ), UU( NUMU, NTAU, NPHI ) )
!DTAUC = 0.0; SSALB = 0.0; PMOM = 0.0; TEMPER = 0.0; UTAU = 0.0; UMU = 0.0; PHI = 0.0;
!H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
!RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; DFDT = 0.0; UAVG = 0.0; UU = 0.0;
!ALBMED = 0.0; TRNMED = 0.0;

!end subroutine allocate_disort_allocatable_arrays

!subroutine deallocate_disort_allocatable_arrays()

!deallocate( DTAUC, SSALB, PMOM, TEMPER, UTAU, UMU, PHI, H_LYR )
!deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )
!deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED, UU )

!end subroutine deallocate_disort_allocatable_arrays
