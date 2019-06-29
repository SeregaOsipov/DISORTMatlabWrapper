c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: ErrPack.f,v 2.1 2000/03/27 21:40:49 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /

      character*20000 line
      integer*4 k
      integer*4 mexPrintf
      integer x


      IF ( FATAL )  THEN
         WRITE ( line, * )  ' ******* ERROR >>>>>>  ', MESSAG         
         call mexErrMsgIdAndTxt('DisortMatlab:ErrMsg',line)     
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE (line, * )  ' ******* WARNING >>>>>>  ', MESSAG
         k=mexPrintf(line//achar(13))
      ELSE
         WRITE ( *,99 )
         k=mexPrintf('99'//achar(13))
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     &   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /

      character*20000 line
      integer*4 k
      integer*4 mexPrintf
      integer x


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( line, '(3A)' )  ' ****  Input variable  ', VarNam,
     &                     '  in error  ****'
      k=mexPrintf(line//achar(13))
      IF ( NumMsg.EQ.MaxMsg )
     &   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal

      character*20000 line
      integer*4 k
      integer*4 mexPrintf
      integer x

      write(line,*) DimNam,' should be increased to at least ',MinVal      
c      call mexErrMsgIdAndTxt('SymbolicDimension:Size',line)     
c      WRITE (line, '(/,3A,I7)' )  '****Symbolic dimension ', DimNam,
c     &                  '  should be increased to at least ', MinVal
      k=mexPrintf(line//achar(13))
      WrtDim = .TRUE.

      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr

      character*20000 line
      integer*4 k
      integer*4 mexPrintf
      integer x


      TstBad = .FALSE.
      WRITE( line, '(/,3A,1P,E11.2,A)' )
     &       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     &       ' per cent from correct value.  Self-test failed.'
      k=mexPrintf(line//achar(13))
      RETURN
      END

