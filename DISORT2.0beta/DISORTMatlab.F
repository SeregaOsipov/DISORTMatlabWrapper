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
      integer*4 ::  NLYR, NMOM, NTAU, NSTR, NUMU, NPHI, IBCND,
     &          MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM
      logical USRTAU, USRANG, LAMBER, PLANK, ONLYFL
      real*4  WVNMLO, WVNMHI, FBEAM, UMU0, PHI0, FISOT, ALBEDO, BTEMP, 
     &      TTEMP, TEMIS, ACCUR

c     .. array arguments ..
      logical PRNT(5)
      real*4, dimension(:), allocatable::  DTAUC,SSALB,TEMPER,UTAU,UMU,
     &          PHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED
      real*4, dimension(:,:), allocatable:: PMOM
      real*4, dimension(:,:,:), allocatable:: UU

c     temproray storage for MATLAB logical items      
      integer*1 logicalItem


      character*20000 line
      integer*4 k
      integer*4 mexPrintf
      integer x

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
c     TODO:// uncomment the check back
c      if(nrhs .ne. 45) then
c         call mexErrMsgIdAndTxt ('MATLAB:disortMatlab:nInput',
c     +                          '45 input parameters required.')     
c      endif

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
      call mxCopyPtrToInteger4(mxGetPr(prhs(1)),NLYR,1)
      write(line,*) '1, NLYR = ',NLYR
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(4)),NMOM,1)
      write(line,*) '4, NMOM = ',NMOM
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(33)),MAXCLY,1)
      write(line,*) '33, MAXCLY = ',MAXCLY
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(34)),MAXULV,1)
      write(line,*) '34, MAXULV = ',MAXULV
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(35)),MAXUMU,1)
      write(line,*) '35, MAXUMU = ',MAXUMU
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(36)),MAXPHI,1)
      write(line,*) '36, MAXPHI = ',MAXPHI
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger4(mxGetPr(prhs(37)),MAXMOM,1)
      write(line,*) '37, MAXMOM = ',MAXMOM
      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      size = mrows*ncols

      allocate(DTAUC(NLYR))      
      call mxCopyPtrToReal4(mxGetPr(prhs(2)),DTAUC,size)
      write(line,*) '2, DTAUC = ',DTAUC
      k=mexPrintf(line//achar(13)) 

      allocate(SSALB(NLYR))      
      call mxCopyPtrToReal4(mxGetPr(prhs(3)),SSALB,size)
      write(line,*) '3, SSALB = ',SSALB
      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(6))
      ncols = mxGetN(prhs(6))
      size = mrows*ncols

      allocate(TEMPER(0:NLYR))      
      call mxCopyPtrToReal4(mxGetPr(prhs(6)),TEMPER,size)
      write(line,*) '6, TEMPER = ',TEMPER
      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(5))
      ncols = mxGetN(prhs(5))
      size = mrows*ncols
      
      allocate(PMOM(0:NMOM, NLYR))      
      call mxCopyPtrToReal4(mxGetPr(prhs(5)),PMOM,size)
      write(line,*) '5, PMOM = ',PMOM
      k=mexPrintf(line//achar(13))
      
      call mxCopyPtrToReal4(mxGetPr(prhs(7)),WVNMLO,1)
      write(line,*) '7, WVNMLO = ',WVNMLO
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(8)),WVNMHI,1)
      write(line,*) '8, WVNMHI = ',WVNMHI
      k=mexPrintf(line//achar(13)) 

      call mxCopyPtrToInteger1(mxGetPr(prhs(9)),logicalItem,1)
      write(line,*) '9, logicalItem = ',logicalItem
      k=mexPrintf(line//achar(13)) 

c     NTAU could be input/output variable
      call mxCopyPtrToInteger4(mxGetPr(prhs(10)),NTAU,1)
      write(line,*) '10, NTAU = ',NTAU
      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(11))
      ncols = mxGetN(prhs(11))
      size = mrows*ncols

c     UTAU could be input/output variable
      allocate(UTAU(NTAU))    
      call mxCopyPtrToReal4(mxGetPr(prhs(11)),UTAU,size)
      write(line,*) '11, UTAU = ',UTAU
      k=mexPrintf(line//achar(13))        

      call mxCopyPtrToInteger4(mxGetPr(prhs(12)),NSTR,1)
      write(line,*) '12, NSTR = ',NSTR
      k=mexPrintf(line//achar(13)) 
    
c     NUMU could be input/output variable
      call mxCopyPtrToInteger4(mxGetPr(prhs(14)),NUMU,1)
      write(line,*) '14, NUMU = ',NUMU
      k=mexPrintf(line//achar(13)) 

      mrows = mxGetM(prhs(15))
      ncols = mxGetN(prhs(15))
      size = mrows*ncols

c     UMU is output variable
      allocate(UMU(NUMU))    
      call mxCopyPtrToReal4(mxGetPr(prhs(15)),UMU,size)
      write(line,*) '15, UMU = ',UMU
      k=mexPrintf(line//achar(13))     

      call mxCopyPtrToInteger4(mxGetPr(prhs(16)),NPHI,1)
      write(line,*) '16, NPHI = ',NPHI
      k=mexPrintf(line//achar(13))

      mrows = mxGetM(prhs(17))
      ncols = mxGetN(prhs(17))
      size = mrows*ncols

      allocate(PHI(NPHI))      
      call mxCopyPtrToReal4(mxGetPr(prhs(17)),PHI,size)
      write(line,*) '17, PHI = ',PHI
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToInteger4(mxGetPr(prhs(18)),IBCND,1)
      write(line,*) '18, IBCND = ',IBCND
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(19)),FBEAM,1)
      write(line,*) '19, FBEAM = ',FBEAM
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(20)),UMU0,1)
      write(line,*) '20, UMU0 = ',UMU0
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(21)),PHI0,1)
      write(line,*) '21, PHI0 = ',PHI0
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(22)),FISOT,1)
      write(line,*) '22, FISOT = ',FISOT
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(24)),ALBEDO,1)
      write(line,*) '24, ALBEDO = ',ALBEDO
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(25)),BTEMP,1)
      write(line,*) '25, BTEMP = ',BTEMP
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(26)),TTEMP,1)
      write(line,*) '26, TTEMP = ',TTEMP
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(27)),TEMIS,1)
      write(line,*) '27, TEMIS = ',TEMIS
      k=mexPrintf(line//achar(13))

      call mxCopyPtrToReal4(mxGetPr(prhs(30)),ACCUR,1)
      write(line,*) '30, ACCUR = ',ACCUR
      k=mexPrintf(line//achar(13))

      HEADER = 'fortran header'
c      call mxCopyPtrToCharacter(mxGetPr(prhs(32)),HEADER,127)
C     //TODO: header is not read properly but I don't think I need it      
      k=mexPrintf(HEADER//achar(13)) 

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
      USRTAU = .true.
      USRANG = .true.
      LAMBER = .true.
      PLANK = .false.
      ONLYFL = .false.

      PRNT(1) = .false.
      PRNT(2) = .false.
      PRNT(3) = .false.
      PRNT(4) = .false.
      PRNT(5) = .false.

C     Call the computational subroutine.
      k=mexPrintf("before DISORT"//achar(13)) 

      call DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
     &  WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI,
     &  PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO,
     &  BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, HEADER,
     &  MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
     &  FLUP, DFDT, UAVG, UU, ALBMED, TRNMED)

      k=mexPrintf("after disort")

C     Prepare MATLAB output variables
C     Create matrix for the return argument.
      size = NTAU*1
      classid = mxClassIDFromClassName('single')
      complexflag = 0
      ndim = 2
      dims(1) = NTAU
      dims(2) = 1

      plhs(1) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(RFLDIR,mxGetPr(plhs(1)),size)    

      plhs(2) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(RFLDN,mxGetPr(plhs(2)),size)    

      plhs(3) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(FLUP,mxGetPr(plhs(3)),size)    

      plhs(4) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(DFDT,mxGetPr(plhs(4)),size)    

      plhs(5) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(UAVG,mxGetPr(plhs(5)),size)    

      size = MAXUMU*MAXULV*MAXPHI
      ndim = 3
      dims(1) = MAXUMU
      dims(2) = MAXULV
      dims(3) = MAXPHI

      plhs(6) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(UU,mxGetPr(plhs(6)),size)   

      size = NUMU*1
      ndim = 2
      dims(1) = NUMU
      dims(2) = 1

      plhs(7) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(ALBMED,mxGetPr(plhs(7)),size)   

      plhs(8) = mxCreateNumericArray(ndim,dims, classid, complexflag)
      call mxCopyReal4ToPtr(TRNMED,mxGetPr(plhs(8)),size)   
 
      return

      end