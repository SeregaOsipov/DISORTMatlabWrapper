#include "fintrf.h"
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c Author: Sergey Osipov, KAUST, May 2015
c matlab wrapper for disort routine
c for question email to Serega.Osipov@gmail.com
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

C     Pointers to input/output mxArrays:r

C     Array information:
      mwPointer mrows, ncols
      mwSize size

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(3)

C     Arguments for computational routine:
C     .. scalar arguments ..
      CHARACTER HEADER*127
      integer :: NLYR, NMOM, NTAU, NSTR, NUMU, NPHI, IBCND,
     &          MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM
      logical USRTAU, USRANG, LAMBER, PLANK, ONLYFL
      real  WVNMLO, WVNMHI, FBEAM, UMU0, PHI0, FISOT, ALBEDO, BTEMP,
     &      TTEMP, TEMIS, ACCUR

c     .. array arguments ..
      logical PRNT(5)
      real, dimension(:), allocatable::  DTAUC,SSALB,TEMPER,UTAU,UMU,
     &          PHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED
      real, dimension(:,:), allocatable:: PMOM
      real, dimension(:,:,:), allocatable:: UU

c     temproray storage for MATLAB logical items
      integer*1 logicalItem   

      character*20000 line
c     for debugging it is nice to have compact prints
c      character*200 line
      integer k
      integer mexPrintf
      integer x

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 37) then
         call mexErrMsgIdAndTxt ('MATLAB:disortMatlab:nInput',
     +   'exactly 37 input parameters required.')     
      endif

C     Validate inputs
C     TODO:// check that input is single and not is double!
C     Check that the input is a number.
C      if(mxIsNumeric(prhs(1)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:disortMatlab:NonNumeric',
C     +                           'Input must be a number.')
C      endif

C     Get the size of the input array.
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      size = mrows*ncols

C     Create Fortran array from the input argument.
      call mxCopyPtrToInteger8(mxGetPr(prhs(1)),NLYR,1)
c      write(line,*) '1, NLYR = ',NLYR
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(4)),NMOM,1)
c      write(line,*) '4, NMOM = ',NMOM
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(33)),MAXCLY,1)
c      write(line,*) '33, MAXCLY = ',MAXCLY
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(34)),MAXULV,1)
c      write(line,*) '34, MAXULV = ',MAXULV
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(35)),MAXUMU,1)
c      write(line,*) '35, MAXUMU = ',MAXUMU
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(36)),MAXPHI,1)
c      write(line,*) '36, MAXPHI = ',MAXPHI
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(37)),MAXMOM,1)
c      write(line,*) '37, MAXMOM = ',MAXMOM
c      k=mexPrintf(line//achar(13))

      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      size = mrows*ncols

      allocate(DTAUC(NLYR))      
      call mxCopyPtrToReal8(mxGetPr(prhs(2)),DTAUC,size)
c      write(line,*) '2, DTAUC = ',DTAUC
c      k=mexPrintf(line//achar(13))

      allocate(SSALB(NLYR))      
      call mxCopyPtrToReal8(mxGetPr(prhs(3)),SSALB,size)
c      write(line,*) '3, SSALB = ',SSALB
c      k=mexPrintf(line//achar(13))

      mrows = mxGetM(prhs(6))
      ncols = mxGetN(prhs(6))
      size = mrows*ncols

      allocate(TEMPER(0:NLYR))      
      call mxCopyPtrToReal8(mxGetPr(prhs(6)),TEMPER,size)
c      write(line,*) '6, TEMPER = ',TEMPER
c      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(5))
      ncols = mxGetN(prhs(5))
      size = mrows*ncols
      
      allocate(PMOM(0:NMOM, NLYR))      
      call mxCopyPtrToReal8(mxGetPr(prhs(5)),PMOM,size)
c      write(line,*) '5, PMOM = ',PMOM
c      k=mexPrintf(line//achar(13))
      
      call mxCopyPtrToReal8(mxGetPr(prhs(7)),WVNMLO,1)
c      write(line,*) '7, WVNMLO = ',WVNMLO
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(8)),WVNMHI,1)
c      write(line,*) '8, WVNMHI = ',WVNMHI
c      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger1(mxGetPr(prhs(9)),logicalItem,1)
      call mxCopyInteger1ToLogical( logicalItem, USRTAU)
c      write(line,*) '9, USRTAU = ',USRTAU
c      k=mexPrintf(line//achar(13))

c     NTAU could be input/output variable
      if (USRTAU) then
          call mxCopyPtrToInteger8(mxGetPr(prhs(10)),NTAU,1)
c          write(line,*) '10, NTAU = ',NTAU
c          k=mexPrintf(line//achar(13))

          mrows = mxGetM(prhs(11))
          ncols = mxGetN(prhs(11))
          size = mrows*ncols

c         UTAU could be input/output variable
          allocate(UTAU(NTAU))    
          call mxCopyPtrToReal8(mxGetPr(prhs(11)),UTAU,size)
c          write(line,*) '11, UTAU = ',UTAU
c          k=mexPrintf(line//achar(13))    
      else    
c     both NTAU and UTAU are output variables
          NTAU = NLYR+1
          allocate(UTAU(NTAU))
      end if

      call mxCopyPtrToInteger8(mxGetPr(prhs(12)),NSTR,1)
c      write(line,*) '12, NSTR = ',NSTR
c      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger1(mxGetPr(prhs(13)),logicalItem,1)
      call mxCopyInteger1ToLogical( logicalItem, USRANG)
c      write(line,*) '13, USRANG = ',USRANG
c      k=mexPrintf(line//achar(13))
    
c     NUMU and UMU could be input/output variable
      if ( USRANG) then
          call mxCopyPtrToInteger8(mxGetPr(prhs(14)),NUMU,1)
c          write(line,*) '14, NUMU = ',NUMU
c          k=mexPrintf(line//achar(13))

          mrows = mxGetM(prhs(15))
          ncols = mxGetN(prhs(15))
          size = mrows*ncols

          allocate(UMU(NUMU))    
          call mxCopyPtrToReal8(mxGetPr(prhs(15)),UMU,size)
c          write(line,*) '15, UMU = ',UMU
c          k=mexPrintf(line//achar(13))
      else
c     both NUMU and UMU are output variables
          NUMU = MAXUMU
          allocate(UMU(NUMU))  
      endif  

      call mxCopyPtrToInteger8(mxGetPr(prhs(16)),NPHI,1)
c      write(line,*) '16, NPHI = ',NPHI
c      k=mexPrintf(line//achar(13))

      mrows = mxGetM(prhs(17))
      ncols = mxGetN(prhs(17))
      size = mrows*ncols

      allocate(PHI(NPHI))      
      call mxCopyPtrToReal8(mxGetPr(prhs(17)),PHI,size)
c      write(line,*) '17, PHI = ',PHI
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger8(mxGetPr(prhs(18)),IBCND,1)
c      write(line,*) '18, IBCND = ',IBCND
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(19)),FBEAM,1)
c      write(line,*) '19, FBEAM = ',FBEAM
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(20)),UMU0,1)
c      write(line,*) '20, UMU0 = ',UMU0
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(21)),PHI0,1)
c      write(line,*) '21, PHI0 = ',PHI0
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(22)),FISOT,1)
c      write(line,*) '22, FISOT = ',FISOT
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger1(mxGetPr(prhs(23)),logicalItem,1)
      call mxCopyInteger1ToLogical( logicalItem, LAMBER)
c      write(line,*) '23, LAMBER = ',LAMBER
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(24)),ALBEDO,1)
c      write(line,*) '24, ALBEDO = ',ALBEDO
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(25)),BTEMP,1)
c      write(line,*) '25, BTEMP = ',BTEMP
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(26)),TTEMP,1)
c      write(line,*) '26, TTEMP = ',TTEMP
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(27)),TEMIS,1)
c      write(line,*) '27, TEMIS = ',TEMIS
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger1(mxGetPr(prhs(28)),logicalItem,1)
      call mxCopyInteger1ToLogical( logicalItem, PLANK)
c      write(line,*) '28, PLANK = ',PLANK
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger1(mxGetPr(prhs(29)),logicalItem,1)
      call mxCopyInteger1ToLogical( logicalItem, ONLYFL)
c      write(line,*) '29, ONLYFL = ',ONLYFL
c      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal8(mxGetPr(prhs(30)),ACCUR,1)
c      write(line,*) '30, ACCUR = ',ACCUR
c      k=mexPrintf(line//achar(13))

      HEADER = 'dummy header'
c      call mxCopyPtrToCharacter(mxGetPr(prhs(32)),HEADER,127)
C     //TODO: header is not read properly but I don't need it
c      k=mexPrintf(HEADER//achar(13)) 

c     allocate output variables
      allocate(RFLDIR(NTAU))      
      allocate(RFLDN(NTAU))            
      allocate(FLUP(NTAU))            
      allocate(DFDT(NTAU))            
      allocate(UAVG(NTAU))            
      allocate(UU(MAXUMU, MAXULV, MAXPHI))            
      allocate(ALBMED(MAXUMU))            
      allocate(TRNMED(MAXUMU))            

c     TODO:// fix the boolean values, hardcode boolean values just yet
      PRNT(1) = .false.
      PRNT(2) = .false.
      PRNT(3) = .false.
      PRNT(4) = .false.
      PRNT(5) = .false.

C     Call the computational subroutine.
c      k=mexPrintf("in Wrapper: before DISORT"//achar(13))
c      return
      call DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
     &  WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI,
     &  PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO,
     &  BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, HEADER,
     &  MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
     &  FLUP, DFDT, UAVG, UU, ALBMED, TRNMED)

c      k=mexPrintf("in Wrapper: after disort"//achar(13))
C     Prepare MATLAB output variables
C     Create matrix for the return argument.
      size = NTAU*1
      classid = mxClassIDFromClassName('double')
      complexflag = 0
      ndim = 2
      dims(1) = NTAU
      dims(2) = 1

      plhs(1) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(RFLDIR,mxGetPr(plhs(1)),size)

      plhs(2) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(RFLDN,mxGetPr(plhs(2)),size)

      plhs(3) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(FLUP,mxGetPr(plhs(3)),size)

      plhs(4) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(DFDT,mxGetPr(plhs(4)),size)

      plhs(5) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(UAVG,mxGetPr(plhs(5)),size)

      size = MAXUMU*MAXULV*MAXPHI
      ndim = 3
      dims(1) = MAXUMU
      dims(2) = MAXULV
      dims(3) = MAXPHI

      plhs(6) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(UU,mxGetPr(plhs(6)),size)

      size = NUMU*1
      ndim = 2
      dims(1) = NUMU
      dims(2) = 1

      plhs(7) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(ALBMED,mxGetPr(plhs(7)),size)

      plhs(8) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal8ToPtr(TRNMED,mxGetPr(plhs(8)),size)

c     independent of USRANG always return UMU
c     it will the copy of the provided UMU or from the DISORT
      plhs(9) = mxCreateNumericArray(ndim,dims, classid,
     &                                       complexflag)
      call mxCopyReal8ToPtr(UMU,mxGetPr(plhs(9)),size)

c     independent of USRTAU always return UTAU
c     it will the copy of the provided UTAU or from the DISORT
      size = NTAU*1
      ndim = 2
      dims(1) = NTAU
      dims(2) = 1
      plhs(10) = mxCreateNumericArray(ndim,dims, classid,
     &                                       complexflag)
      call mxCopyReal8ToPtr(UTAU,mxGetPr(plhs(10)),size)
 
      return
      end

c     this routine converts a matlab logical to fortran logical, matrices are not implemented yet
      subroutine mxCopyInteger1ToLogical( logicaldata, fortran)
      implicit none
      integer*1, intent(in) :: logicaldata
      logical, intent(out) :: fortran
      fortran = (logicaldata /= 0)
      return
      end subroutine
