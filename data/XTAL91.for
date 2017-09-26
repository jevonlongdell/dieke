C
C        ***************************************************
C        *             XTAL91                              *
C        * MODIFIED FROM XTAL84.1, FOR VAX/VMS SYSTEM      *
C        * ADDITIONAL SUBROUTINS: XTAL91SUB.FOR            *
C        * ADDITIONAL UNIT NEEDED: 5 XTAL91.DAT            *
C        *                         6 XTAL91.OUT            *
C        *            JUNE 10,91                           *
C        *                                                 * 
C        ***************************************************
C        ***************************************************
C        *                                                 *
C        *              XTAL84.1                           *
C        *     PROGRAM TO WRITE COMPLETE CRYSTAL TAPES     *
C        *         OR TRUNCATED TAPES                      *
C        *                                                 *
C        *    NO MATRIX SIZE CHECK--XTLCF5                 *
C        *                                                 *
C        *    13 FREE-ION SUBMATRICES CAN BE USED          *
C        *    23 FREE ION PARAMETERS MAY BE USED 6/83      *
C        *    20 LEVELS PER SUBMATRIX MAY BE KEPT          *
C        *    BKQ AND BK-Q HAVE THE SAME PAR NO   6/83     *
C        *    PROVISION FOR Q=3 PROVIDED  2/84             *
C        *    STATE NAMES WRITTEN ON UNIT 16   4/84        *
C        *    FORMATTED OUTPUT  4(I9,F11.6)    4/84        *
C        *    OPTION TO USE V1K (WITH S-FACTOR INCLUDED    *
C        *       IN MATRIX ELEMENT)            4/84        *
C        *                                                 *
C        *  STACKING                                       *
C        *                                                 *
C        *    CARD 1  TITLE  (10A8)                        *
C        *    CARD 2  CODE CARD (FREE FORM INPUT--8 NOS.)  *
C        *         NF = NO. OF F-ELECTRONS                 *
C        *         NPAR=NO. FREE-ION PARAMETERS            *
C        *         NPXT = TOTAL NUMBER OF PARAMETERS       *
C        *          (NOT INCLUDING V1K)                    *
C        *         NTRC = 0 IF COMPLETE CRYSTAL            *
C        *              MATRIX IS CALCULATED               *
C        *              .GT.0 IF MATRIX IS TO BE           *
C        *               TRANSFORMED                       *
C        *         MXNO = NO. CRYSTAL SUBMATRICES          *
C        *         NQ = Q MINIMUM                          *
C        *         NLOOP=NO. OF SUBMATRICES FOR EACH MU    *
C        *            (UP TO 3)                            *
C        *         NV1KOP = 0 (UK ONLY) OR 1 (V1K ALSO)    *
C        *                                                 *
C        * * *        IF NTRC IS NOT 0                 * * *
C        *    CARD 3  NUMBER TO KEEP PER SUBMATRIX         *
C        *    CARD 4  NUMBER TO SKIP PER SUBMTRIX          *
C        *    CARD GROUP 5                                 *
C        *         FREE ION PARAMETERS (I3,12X,F12.3,A8)   *
C        *                                                 *
C        *         UNITS USED                              *
C        *    8  SCRATCH                                   *
C        *   15  FREE-ION MATRIX ELEMENTS                  *
C        *   14  FN U(K)                                   *
C        *   16  NEW CRYSTAL TAPE  (INCLUDES STATE NAMES)  *
C        *                                                 *
C        * * *       IF NTRC IS NOT 0                 * *  *
C        *   20  SCRATCH FOR UNTRANSFORMED UK (LARGE)      *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(50,50),Z(50,20),NDIM(13),TITLE(10),NOKEEP(3,13)
      REAL*4 STC(13,50),COEFFG(13,12,21,21)
      DIMENSION D(50),E1(50),E(50),IND(20),SCR1(50),SCR2(50),
     1SCR3(50),SCR4(50),SCR5(50),EVAL(20),NOSKIP(3,13)
      COMMON/COM1/RVEC(13,50,20),COEFFG
      COMMON/COM2/MXNO,NQ,NPAR,NPXT
      COMMON/COM3/NOM(13,50),MUL(119),LST(119)
C
      IL = 20
      NEIGEN = 50
      NKPM = 20
C
      READ (5,1000) TITLE
      WRITE (6,1000) TITLE
C      PRINT 1000, TITLE
 1000 FORMAT (10A8)
 1001 FORMAT (24I3)
  101 FORMAT (1H )
      READ(5, *) NF,NPAR,NPXT,NTRC,MXNO,NQ,NLOOP,NV1KOP
C      PRINT *, 'NLOOP =',NLOOP
      IF (NLOOP.EQ.0) NLOOP = 1
      MMM = 0
      XADD = 0.0
      NFN = NF/2
      NF2 = 2*NFN
      IF (NF2.NE.NF) XADD = 0.5
C
C               READ IN TERMS OF CONFIGURATION, UK (AND V1K)
C
C             CALCULATE CRYSTAL REDUCED MATRIX ELEMENTS FOR COMPLETE
C                  CONFIGURATION, WRITE ON 8
C                   OR ON 20 IF NTRC.GT.0
C
      CALL RME84(H,MMM,NEIGEN,NF,NDIM,STC,NTRC,NSUB,NV1KOP,XADD)
      IF (MMM.NE.0) GO TO 345
C
C      PRINT *, 'PASS 1: AFTER CALL RME84'
C      PRINT *, 'NSUB =',NSUB
      DO 3 K=1,NSUB
      NOSKIP(1,K) = 0
    3 NOKEEP(1,K) = NDIM(K)
C
C    READ NUMBER OF LEVELS TO KEEP IF NTRC IS GREATER THAN 0
C
      IF (NTRC.EQ.0) GO TO 25
      DO 5 NL=1,NLOOP
      READ (5,1001) (NOKEEP(NL,MT),MT=1,NSUB)
      READ (5,1001) (NOSKIP(NL,MT),MT=1,NSUB)
C      PRINT 1001, (NOKEEP(NL,MT),MT=1,NSUB)
C      PRINT 1001, (NOSKIP(NL,MT),MT=1,NSUB)
    5 CONTINUE
C
C             READ IN FREE-ION MATRIX ELEMENTS
C                  TRANSFORM IF NTRC IS NOT 0
C
   25 CONTINUE
      NP2 = NPAR/2
      N2 = NP2*2
      IF (N2.NE.NPAR) NP2 = NP2 + 1
      DO 35 MT=1,NSUB
      NKPT = NOKEEP(NLOOP,MT) + 1
      DO 30 I=1,NKPT
      DO 28 J=1,NKPT
      DO 26 K=1,NP2
   26 COEFFG(MT,K,I,J) = 0.0
   28 CONTINUE
   30 CONTINUE
   35 CONTINUE
C
C  FREE ION MATRIX ELEMENTS ARE READ FROM TAPE 15
C
      CALL FRICF(NEIGEN,NKPM,NTRC,H,Z,MMM,D,E,E1,IND,SCR1,SCR2,
     1SCR3,SCR4,SCR5,EVAL,STC,NSUB,NOKEEP,NLOOP)
      IF (MMM.NE.0) GO TO 345
C      PRINT *, 'PASS 2: AFTER FRICF'
C
C             *****************************************
C             *                                       *
C             *    MAKE CRYSTAL TAPE                  *
C             *                                       *
C             *****************************************
C
      IF (NTRC.EQ.0) GO TO 55
C
C                  TRANSFORM CRYSTAL REDUCED MATRIX ELEMENTS
      CALL RMTR84(H,Z,NEIGEN,NKPM,NOKEEP,NLOOP)
C      PRINT *, 'PASS 3: AFTER RMTR84'
C
   55 CALL XLTP84(MMM,H,NEIGEN,NKPM,STC,XADD,NSUB,NOKEEP,NOSKIP,
     1NLOOP,NTRC,NV1KOP)
C
C      PRINT *, 'PASS 4: AFTER XLTP84'
      IF (MMM.EQ.0) GO TO 60
C  345 PRINT 1345, MMM
  345 WRITE (6,1345) MMM
 1345 FORMAT (1X,' ERROR AT ',I6)
   60 CONTINUE
      REWIND 15
      REWIND 16
   65 CONTINUE
      STOP
      END
C
C ------
C
      SUBROUTINE RME84(H,MMM,NEIGEN,NF,NDIM,STC,NTRC,NSUB,
     1NV1KOP,XADD)
C
C        ***************************************************
C        *                                                 *
C        *    CALCULATE CRYSTAL REDUCED MATRIX ELEMENTS    *
C        *    FOR COMPLETE CONFIGURATIONS                  *
C        *                                                 *
C        *    USE WITH XTAL84 AND XINT84                   *
C        *    READ IN STATE NAMES FOR FREE ION SUBMATRICES *
C        *      AND SUBMATRIX DIMENSIONS FROM UNIT 14      *
C        *    OPTION TO PUT IN V1K (IN NV1KOP)             *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 STC(13,50),ASTC(119),U4(6,119,119)
      DIMENSION H(NEIGEN,NEIGEN),NDIM(13)
      COMMON/COM1/HZ(50,50),U4,ASTC,NSCR(4711)
      COMMON/COM3/NOM(13,50),MUL(119),LST(119)
C
      RTFAC = 1.
      IF (NF.GT.7) RTFAC = -1.
      IK1 = 8
      IF (NTRC.NE.0) IK1=20
      REWIND IK1
      REWIND 14
      READ (14,1001) NF,NOST,NSUB,(NDIM(MT),MT=1,NSUB)
 1001 FORMAT (20I4)
C
C   READ IN TERMS OF THE CONFIGURATION
C
      READ (14,1002) (ASTC(I),I=1,NOST)
 1002 FORMAT (20(1X,A3))
      READ (14,1003) (MUL(I),LST(I),I=1,NOST)
 1003 FORMAT (20(1X,I1,I2))
C
C READ IN STATE NAMES FOR J SUBMATRICES
C
      DO 10 MT=1,NSUB
      ND = NDIM(MT)
      READ (14,1004) (STC(MT,I),I=1,ND)
   10 CONTINUE
 1004 FORMAT (6(A3,9X))
C
C    READ IN NON-ZERO U(K) (AND V(1K) IF NV1KOP>0) K=2,4,6
C
      KEND = 3
      IF (NV1KOP.NE.0) KEND = 6
      DO 16 J=1,NOST
      DO 15 I=1,NOST
      DO 14 K=1,KEND
   14 U4(K,I,J) = 0.0
   15 CONTINUE
   16 CONTINUE
   20 READ (14,1005,END=22) I,J,(U4(K,I,J),K=1,KEND)
      GO TO 20
 1005 FORMAT (2I4,6F12.7)
   22 CONTINUE
      REWIND 14
      DO 30 MT=1,NSUB
      ND1 = NDIM(MT)
      DO 28 I=1,ND1
      DO 25 J=1,NOST
      IF (STC(MT,I).EQ.ASTC(J)) GO TO 27
   25 CONTINUE
      WRITE (6, 2000) I,MT,STC(MT,I)
C      PRINT 2000, I,MT,STC(MT,I)
 2000 FORMAT (2I3,1X,A3)
      MMM = 25
      WRITE (6, 2001) MMM
C      PRINT 2001, MMM
 2001 FORMAT (1X,' ERROR IN RME, CODE = ',I6)
      RETURN
   27 NOM(MT,I) = J
   28 CONTINUE
   30 CONTINUE
C
C      MAIN CALCULATION
C
      DO 80 MT=1,NSUB
      AJ = MT-1
      AJ = AJ + XADD
      NJ = 2.*AJ
      ND1 = NDIM(MT)
      DO 78 MT1 = MT,NSUB
      IF((MT1-MT).GT.6) GO TO 80
      AJ1 = MT1-1
      AJ1 = AJ1 + XADD
      NJ1 = 2.*AJ1
      FACJ = DSQRT((2.*AJ+1.)*(2.*AJ1+1.))
      ND2 = NDIM(MT1)
      NJX = AJ + AJ1
      DO 75 K=2,6,2
      IF (NJX.LT.K) GO TO 75
      IF ((MT1-MT).GT.K) GO TO 75
      K1 = K/2
      NK = 2*K
      K1V = K1+3
      K10 = 10 + K
      DO 70 I=1,ND1
      I1 = NOM(MT,I)
      NA = 2*LST(I1)
      NB = MUL(I1)-1
      S = NB
      S = 0.5*S
      SFAC = DSQRT(S*(S+1.D0)/(2.*S + 1.D0))
      DO 65 J=1,ND2
      H(I,J) = 0.0
      HZ(I,J) = 0.0
      J1 = NOM(MT1,J)
      IF (MUL(I1).NE.MUL(J1)) GO TO 65
      NA1 = 2*LST(J1)
      IF1 = (NB+NA1+NJ+NK)/2
      RL = SIXJ(NJ,NJ1,NK,NA1,NA,NB)
      H(I,J) = (-1.)**IF1*RL*U4(K1,I1,J1)*FACJ*RTFAC
      IF (NV1KOP.EQ.0) GO TO 65
      IF (NB.EQ.0) GO TO 65
      HZ(I,J) =(-1.)**IF1*RL*U4(K1V,I1,J1)*FACJ*SFAC
C   NOTE THAT SFACTOR IS INCLUDED
   65 CONTINUE
   70 CONTINUE
      WRITE (IK1) MT,ND1,MT1,ND2,K,((H(I,J),J=1,ND2),I=1,ND1)
      IF (NV1KOP.NE.0)
     1WRITE (IK1) MT,ND1,MT1,ND2,K10,((HZ(I,J),J=1,ND2),I=1,ND1)
C     IF (MT1.GT.3) GO TO 75
C     PRINT 101
  101 FORMAT (1H )
C     PRINT 1001, MT,ND1,MT1,ND2,K,K10
C  PRINT CHECK FOR REDUCED UK
C     N1 = ND1
C     N2 = ND2
C     DO 68 I=1,N1
C  68 PRINT 1006,(H(I,J),J=1,N2)
C     PRINT 101
C     DO 69 I=1,N1
C  69 PRINT 1006,(HZ(I,J),J=1,N2)
C     PRINT 101
 1006 FORMAT (10F12.6)
   75 CONTINUE
   78 CONTINUE
   80 CONTINUE
      END FILE IK1
      REWIND IK1
      RETURN
      END
C
C ------
C
      SUBROUTINE FRICF(NEIGEN,NKPM,NTRC,H,Z,MMM,D,E,E1,IND,SCR1,
     1SCR2,SCR3,SCR4,SCR5,EVAL,STC,NSUB,NOKEEP,NLOOP)
C
C        ***************************************************
C        *                                                 *
C        *    READ FREE-ION MATRIX ELEMENTS,               *
C        *         TRANSFORM IF NECESSARY                  *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NEIGEN,NEIGEN),Z(NEIGEN,NKPM),EVAL(NKPM)
      DIMENSION PARM(23),NOKEEP(3,13)
      REAL*4 AME(12,51,51),COEFFG(13,12,21,21),STC(13,50)
      DIMENSION D(NEIGEN),E(NEIGEN),E1(NEIGEN),IND(NKPM),SCR1(NEIGEN),
     1SCR2(NEIGEN),SCR3(NEIGEN),SCR4(NEIGEN),SCR5(NEIGEN)
      COMMON/COM1/RVEC(13,50,20),COEFFG
      COMMON/COM2/MXNO,NQ,NPAR,NPXT
C
C             *****************************************
C             *                                       *
C             *    FREE-ION ME ARE READ FROM UNIT 15  *
C             *    ODD-NUMBERED PARAMETERS ARE PUT    *
C             *      IN UPPER RIGHT OF COEFFG AND AME *
C             *    EVEN NUMBERED PARAMETERS ARE IN    *
C             *      LOWER LEFT OF COEFFG AND AME     *
C             *                                        *
C             *****************************************
C
      REWIND 15
      IF (NTRC.EQ.0) GO TO 15
C
C      READ PARAMETERS
C
      WRITE (6, 2000)
C      PRINT 2000
 2000 FORMAT (1X,' PARAMETERS FOR TRUNCATION')
      WRITE (6, 101)
C      PRINT 101
  101 FORMAT (1H )
      DO 10 I=1,NPAR
C      PRINT *,'PASS DO 1'
      READ (5,1002) K,PARM(K),PARN
   10 WRITE (6, 2001) K,PARM(K),PARN
C   10 PRINT 2001, K,PARM(K),PARN
      WRITE (6,101)
C      PRINT 101
   15 READ (15,1009) MT,NDT
      NKPT = NOKEEP(NLOOP,MT)
      IF (NKPT.LE.NKPM) GO TO 17
      MMM = 15
      WRITE (6,1345) MMM,NKPT,NKPM
C      PRINT 1345, MMM,NKPT,NKPM
 1345 FORMAT (1X,' ERROR IN FRICF, CODE = 15, NKPT =',I6,' MAX '
     1,'ALLOWED IS',I6)
      RETURN
   17 IF (NTRC.EQ.0) GO TO 30
      NP2P = NPAR/2
      N2 = 2*NP2P
      IF (N2.NE.NPAR) NP2P = NP2P + 1
      NDT1 = NDT+1
      DO 22 I=1,NDT1
      DO 21 J=1,NDT1
      DO 20 K=1,NP2P
   20 AME(K,I,J) = 0.0
C      PRINT *,'PASS DO 2'
   21 CONTINUE
   22 CONTINUE
      DO 25 I=1,NDT
      DO 24 J=1,NDT
      H(I,J) = 0.0
C      PRINT *,'PASS DO 3'
   24 CONTINUE
   25 CONTINUE
   30 CONTINUE
      WRITE (6, 1009) MT,NDT
C      PRINT 1009,MT,NDT
   50 READ (15,1009) M1,I,J,NP,X
C      PRINT *, M1,I,J,NP,X
      IF (M1.EQ.0) GO TO 55
      IF (NKPT.EQ.0) GO TO 50
      I1 = I
      J1 = J
      NP2 = NP/2
      NP22 = 2*NP2
      NP2 = NP2 + 1
      IF (NP22.NE.NP) GO TO 52
      J1 = I
      I1 = J+1
      NP2 = NP2 -1
   52 IF (NTRC.EQ.0) GO TO 54
      AME(NP2,I1,J1) = X
      H(I,J) = H(I,J) + PARM(NP)*X
      GO TO 50
   54 COEFFG(MT,NP2,I1,J1) = X
      GO TO 50
   55 CONTINUE
      IF (NTRC.EQ.0) GO TO 92
C
C      TRANSFORM AND TRUNCATE
      IF (NKPT.EQ.0) GO TO 92
C
      NDT1 = NDT-1
      DO 65 I=1,NDT1
      I1 = I + 1
      DO 64 J=I1,NDT
   64 H(J,I) = H(I,J)
C      PRINT *, 'PASS DO 4'
   65 CONTINUE
      NVPT = 0
C      PRINT *, 'PASS DO 44'
      CALL VECTR(H,Z,EVAL,NDT,NKPT,NVPT,MMM,NEIGEN,NKPM,D,E,E1,IND,
     1SCR1,SCR2,SCR3,SCR4,SCR5)
C      PRINT *, 'PASS DO 444'
      IF (MMM.GT.0) GO TO 95
C
C      PRINT VECTORS
      CALL VCPR(Z,EVAL,NDT,NKPT,MT,NEIGEN,NKPM,STC)
      DO 70 I=1,NDT
      DO 68 J=1,NKPT
      RVEC(MT,I,J) = Z(I,J)
C      PRINT *, 'PASS DO 5'
   68 CONTINUE
   70 CONTINUE
      DO 90 L=1,NP2P
      NP2 = L
C  NP = 2*L-1
      N1 = NDT-1
      NN = 1
      DO 81 I=1,NDT
      DO 80 J=I,NDT
   80 H(I,J) = AME(L,I,J)
C      PRINT *, 'PASS DO 6'
   81 CONTINUE
   82 DO 84 I=1,N1
      I1 = I+1
      DO 83 J=I1,NDT
   83 H(J,I) = H(I,J)
   84 CONTINUE
      CALL COEFTR(H,Z,MT,NDT,NKPT,MT,NDT,NKPT,NEIGEN,NKPM)
      GO TO (85,75),NN
   85 DO 87 I=1,NKPT
      DO 86 J=I,NKPT
   86 COEFFG(MT,NP2,I,J) = H(I,J)
   87 CONTINUE
      NP = 2*L
      IF (NP.GT.NPAR) GO TO 92
      NN = 2
      DO 89 I=1,NDT
      J1 = I
      DO 88 J=I,NDT
      I1 = J + 1
   88 H(I,J) = AME(L,I1,J1)
   89 CONTINUE
      GO TO 82
   75 DO 78 I=1,NKPT
      J1 = I
      DO 77 J=I,NKPT
      I1 = J + 1
   77 COEFFG(MT,NP2,I1,J1) = H(I,J)
   78 CONTINUE
   90 CONTINUE
   92 CONTINUE
      IF (MT.LT.NSUB) GO TO 15
   95 CONTINUE
 1002 FORMAT (I3,12X,F12.3,A8)
 1009 FORMAT (4I6,F12.9)
 2001 FORMAT (I3,12X,F12.3,A8)
      RETURN
      END
C
C ------
C
      SUBROUTINE VCPR(Z,EVAL,NDT,NKPT,MT,NEIGEN,NKPM,STC)
C
C        ***************************************************
C        *                                                 *
C        *    PRINT VECTORS                                *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(NEIGEN,NKPM),EVAL(NKPM)
      REAL*4 STC(13,50)
      WRITE (6, 101)
C      PRINT 101
  101 FORMAT (1H )
      WRITE (6, 2002) MT
C      PRINT 2002, MT
 2002 FORMAT (1X,' VECTORS FOR MATRIX',I3)
      WRITE (6, 101)
C      PRINT 101
      DO 10 K=1,NKPT,10
      KK = K + 9
      IF (KK.GT.NKPT) KK = NKPT
      WRITE (6, 2003) (EVAL(I),I=K,KK)
C      PRINT 2003, (EVAL(I),I=K,KK)
      WRITE (6, 101)
C      PRINT 101
 2003 FORMAT (12X,10(1X,F9.2,1X))
      DO 12 J=1,NDT
      WRITE (6, 2004) STC(MT,J),(Z(J,I1),I1=K,KK)
C      PRINT 2004, STC(MT,J),(Z(J,I1),I1=K,KK)
   12 CONTINUE
   10 CONTINUE
 2004 FORMAT (1X,A3,6X,10F11.7)
      WRITE (6, 101)
C      PRINT 101
      RETURN
      END
C
C ------
C
      SUBROUTINE VECTR(H,Z,EVAL,ND,NKP,NVPT,MMM,NGN,NPM,D,E,E1,IND,
     1SCR1,SCR2,SCR3,SCR4,SCR5)
C
C        ***************************************************
C        *                                                 *
C        *    CALCULATE EIGENVECTORS AND EIGENVALUES       *
C        *         ND=DIMENSION                            *
C        *         NKP = NUMBER DESIRED                    *
C        *         NVPT = 1 FOR PRINTING LARGEST COMPONENT *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SCR1(NGN),SCR2(NGN),SCR3(NGN),SCR4(NGN),SCR5(NGN)
      DIMENSION E(NGN),E1(NGN),D(NGN),IND(NPM),EVAL(NPM),H(NGN,NGN),
     1Z(NGN,NPM)
      NEIGEN=NGN
      EPS1 = -1
      M11 = 1
      CALL TRED1(NEIGEN,ND,H,D,E,E1)
C      PRINT *, 'PASS CALL 1'
C      PRINT *, 'NKP =', NKP
      CALL TRIDIB(ND,EPS1,D,E,E1,BL,UB,M11,NKP,EVAL,IND,IERR,SCR4,SCR5)
C      PRINT *, 'PASS CALL 2'
      IF (IERR.EQ.0) GO TO 10
    5 WRITE (6,1012) IERR
C    5 PRINT 1012, IERR
      MMM = 5
      RETURN
   10 CALL TINVIT(NEIGEN,ND,D,E,E1,NKP,EVAL,IND,Z,IERR,SCR1,SCR2,SCR3,
     1SCR4,SCR5)
C      PRINT *, 'PASS CALL 3'
      IF (IERR.EQ.0) GO TO 15
      WRITE (6, 1014) IERR
C      PRINT 1014, IERR
      MMM = 10
      RETURN
   15 CALL TRBAK1(NEIGEN,ND,H,E,NKP,Z)
C      PRINT *, 'PASS CALL 4'
      IF (NVPT.EQ.0) GO TO 30
      WRITE (6, 1010)
C      PRINT 1010
 1010 FORMAT (1X,' VECTORS CALCULATED')
      DO 25 I=1,NKP
      MM = 1
      TM = DABS(Z(1,I))
      DO 20 J=2,ND
      IF ((DABS(Z(J,I))-TM).LE.0.01) GO TO 20
      MM = J
      TM = DABS(Z(J,I))
   20 CONTINUE
      WRITE (6, 1020) I,EVAL(I),MM,Z(MM,I)
C      PRINT 1020, I,EVAL(I),MM,Z(MM,I)
   25 CONTINUE
      WRITE (6, 101)
C      PRINT 101
   30 CONTINUE
  101 FORMAT (1H )
 1012 FORMAT (1X,' ERROR IN TRIDIB, CODE = ',I5)
 1014 FORMAT (1X,' ERROR IN TINVIT, CODE = ',I5)
 1020 FORMAT (I6,F12.2,I6,F12.8)
      RETURN
      END
C
C ------
C
      SUBROUTINE WRT84(COEFF,INDC,NCOEFF,NNTT)
C
C        ***************************************************
C        *                                                 *
C        *    ROUTINE FOR WRITING COEFFICIENTS TAPES       *
C        *    IN FORMATTED FORM (4(I9,F11.6))              *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION COEFF(4),INDC(4)
      NNTT = NNTT + NCOEFF
      WRITE (16,1015) (INDC(I),COEFF(I),I=1,4)
 1015 FORMAT (4(I9,F11.6))
      DO 10 I=1,4
      INDC(I) = 0
   10 COEFF(I) = 0.D0
      NCOEFF = 0
      RETURN
      END
C
C -----
C
      SUBROUTINE COEFTR(H,Z,MT1,ND1,NKP1,MT2,ND2,NKP2,NEIGEN,NKPM)
C
C        ***************************************************
C        *                                                 *
C        *    TRANSFORM COEFFICIENTS USING RVEC FROM COM1  *
C        *    PUT INTERMEDIATE RESULTS IN Z, OUTPUT IN H   *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 COEFFG(13,12,21,21)
      DIMENSION H(NEIGEN,NEIGEN),Z(NEIGEN,NKPM)
      COMMON/COM1/RVEC(13,50,20),COEFFG
      DO 20 I=1,ND1
      DO 15 J=1,NKP2
      Z(I,J) = 0.0
      DO 10 K=1,ND2
      Z(I,J) = Z(I,J) + H(I,K)*RVEC(MT2,K,J)
   10 CONTINUE
   15 CONTINUE
   20 CONTINUE
      DO 35 I=1,NKP1
      DO 30 J=1,NKP2
      H(I,J) = 0.0
      DO 25 K=1,ND1
      H(I,J) = H(I,J) +RVEC(MT1,K,I)*Z(K,J)
   25 CONTINUE
   30 CONTINUE
   35 CONTINUE
      RETURN
      END
C
C -----
C
      SUBROUTINE RMTR84(H,Z,NEIGEN,NKPM,NOKEEP,NLOOP)
C
C        ***************************************************
C        *                                                 *
C        *    TRANSFORM CRYSTAL REDUCED MATRIX ELEMENTS    *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NEIGEN,NEIGEN),Z(NEIGEN,NKPM)
      REAL*4 COEFFG(13,12,21,21)
      DIMENSION NOKEEP(3,13)
      COMMON/COM1/RVEC(13,50,20),COEFFG
   10 READ (20,END=25) MT,ND1,MT1,ND2,K,((H(I,J),J=1,ND2),I=1,ND1)
      NKP1 = NOKEEP(NLOOP,MT)
      NKP2 = NOKEEP(NLOOP,MT1)
C     PRINT 1001, MT,NKP1,MT1,NKP2,K
      IF (NKP1.EQ.0.OR.NKP2.EQ.0) GO TO 10
       CALL COEFTR(H,Z,MT,ND1,NKP1,MT1,ND2,NKP2,NEIGEN,NKPM)
      WRITE (8) MT,NKP1,MT1,NKP2,K,((H(I,J),J=1,NKP2),I=1,NKP1)
      GO TO 10
   25 CONTINUE
      WRITE (6,1001) MT,NKP1,MT1,NKP2,K
C      PRINT 1001, MT,NKP1,MT1,NKP2,K
 1001 FORMAT (24I3)
      WRITE (6,2002)
C      PRINT 2002
 2002 FORMAT (1X,' REDUCED MATRIX ELEMENTS TRANSFORMED ON UNIT  8' )
      RETURN
      END
C
C ------
C
      SUBROUTINE XLTP84(MMM,H,NEIGEN,NKPM,STC,XADD,NSUB,NOKEEP,
     1NOSKIP,NLOOP,NTRC,NV1KOP)
C
C        ***************************************************
C        *                                                 *
C        *    MAKE CRYSTAL TAPE FROM PREPARED FREE-ION AND *
C        *      REDUCED CRYSTAL MATRIX ELEMENTS            *
C        *       FORMATTED OUTPUT                          *
C        *                                                 *
C        ***************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NEIGEN,NEIGEN),COEFF(4),INDC(4),C(3),
     1NOKEEP(3,13),NOSKIP(3,13),HZ(50,50)
      DIMENSION NUM(13),NXM(6),NLL(6),NAJ1(6),KNQ(3)
      REAL*4 COEFFG(13,12,21,21),SLL(6),STC(13,50)
      COMMON/COM1/RVEC(13,50,20),COEFFG
      COMMON/COM2/MXNO,NQ,NPAR,NPXT
      EQUIVALENCE (HZ(1),RVEC(1))
C
      KNQ(1) = 2/NQ
      KNQ(2) = 4/NQ
      KNQ(3) = 6/NQ
      C(1) = -1.366260
      C(2) = 1.128152
      C(3) = -1.277381
      XMLP = IABS(NQ)
      XDEL = 0.0
      IF (DABS(XADD).LT.0.0001) XDEL=-0.4
      REWIND 14
      DO 700 NOLP=1,NLOOP
      MNO = MXNO *(NOLP-1)
      DO 500 MATNO=1,MXNO
      MTN = MNO + MATNO
      NTOT = 0
      AMU = MATNO-1
      AMU = AMU + XADD
C
C     CALCULATE DIMENSION
      NNTT = 0
      NCOEFF = 0
      ND = 0
      DO 25 I=1,NSUB
      ANO = 1.
      NUM(I) = 0
      IF (I.LT.MATNO) GO TO 25
      AJ1 = I-1
      AJ1 = AJ1 + XADD
      MXM = NOKEEP(NOLP,I)-NOSKIP(NOLP,I)
      IF (MXM.EQ.0) GO TO 25
   20 XM = DABS(AMU-ANO*XMLP)
      IF ((AJ1-XM).LT.XDEL) GO TO 22
      ND = ND + MXM
      NUM(I) = NUM(I) + 1
      XM = AMU + ANO*XMLP
      IF ((AJ1-XM).LT.XDEL) GO TO 22
      ND = ND + MXM
      NUM(I) = NUM(I) + 1
      ANO = ANO + 1.
      GO TO 20
   22 NUM(I) = NUM(I) + 1
      ND = ND + MXM
      WRITE (6, 1001) I,NUM(I),ND
C      PRINT 1001, I,NUM(I),ND
 1001 FORMAT (24I3)
   25 CONTINUE
      NDDD = ND
C     IF (NDDD.LE.200) GO TO 30
C     PRINT 7777, NDDD
C     MMM = 25
C7777 FORMAT (1X,' MATRIX IS TOO LARGE')
C     RETURN
   30 CONTINUE
      NTOT = NDDD
C
C      WRITE FIRST RECORD ON COEFFICIENT TAPE 16
C
      WRITE (16,1010) MTN,NTOT
 1010 FORMAT (3I8)
C
C     WRITE STATE LABELS ON UNIT 16  4/84
C
      NP = 1.
      X = 1.
      ND = 0
      LPRT = 0
      DO 65 K=1,NSUB
      IF (K.LT.MATNO) GO TO 65
      NDA = NOSKIP(NOLP,K)
      NM = NOKEEP(NOLP,K)-NDA
      IF (NM.EQ.0) GO TO 65
      ANO = 1.
      AJ1 = K-1
      AJ1 = AJ1 + XADD
      XM = AMU
      NK = 1
      ND1 = ND
   35 DO 45 L=1,NM
      LL = L + NDA
      I = ND1 + L
      LPRT = LPRT + 1
      IF (NTRC.EQ.0) SLL(LPRT) = STC(K,LL)
      IF (NTRC.NE.0) NLL(LPRT) = LL
      NAJ1(LPRT) = AJ1
      NXM(LPRT) = XM
      IF (NTRC.EQ.0) GO TO 41
      WRITE (6, 1026) MTN,I,I,NP,X,LL,AJ1,XM
C      PRINT 1026, MTN,I,I,NP,X,LL,AJ1,XM
      GO TO 42
   41 CONTINUE
      WRITE (6, 1025) MTN  ,I,I,NP,X,SLL(LPRT),AJ1,XM
C      PRINT 1025, MTN  ,I,I,NP,X,SLL(LPRT),AJ1,XM
   42 CONTINUE
      IF (LPRT.LT.6) GO TO 40
      IF (NTRC.EQ.0) GO TO 36
      WRITE (16,6001) (NLL(LX),NAJ1(LX),NXM(LX),LX=1,6)
      GO TO 37
   36 WRITE (16,6002) (SLL(LX),NAJ1(LX),NXM(LX),LX=1,6)
   37 DO 38 LX=1,6
      NLL(LX) = 0
      SLL(LX) = 0
      NXM(LX) = 0
      NAJ1(LX) = 0
   38 CONTINUE
      LPRT = 0
   40 CONTINUE
   45 CONTINUE
      IF (NK.EQ.NUM(K) ) GO TO 60
      NK2 = NK/2
      NKK2 = 2*NK2
      IF (NKK2.EQ.NK) GO TO 55
      XM = AMU -ANO*XMLP
      GO TO 50
   55 XM = AMU + ANO*XMLP
      ANO = ANO + 1
   50 NK = NK + 1
      ND1 = ND1 + NM
      GO TO 35
   60 ND = ND + NUM(K)*NM
   65 CONTINUE
      IF (LPRT.EQ.0) GO TO 72
      IF (NTRC.EQ.0) WRITE (16, 6002) (SLL(LX),NAJ1(LX),NXM(LX),
     1LX=1,6)
      IF (NTRC.NE.0)
     1WRITE (16,6001) (NLL(LX),NAJ1(LX),NXM(LX),LX=1,6)
   72 DO 68 I=1,NTOT
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF) = X
      INDC(NCOEFF) = I + 1000*(I + 1000*NP)
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
   68 CONTINUE
C
C      PUT IN FREE-ION MATRIX ELEMENTS
C
   70 ND = 0
      NDP = 0
      NK = 1
      MT = MATNO
      NP = NP + 1
      IF (NP.GT.NPAR) GO TO 110
      KC = 0
      NP2 = NP/2
      NP22 = 2*NP2
      IF (NP22.EQ.NP) GO TO 75
      KC = 1
      NP2 = NP2 + 1
   75 CONTINUE
      NDA = NOSKIP(NOLP,MT)
      ND1 = NOKEEP(NOLP,MT)-NDA
      IF (ND1.EQ.0) GO TO 91
      IF (MT.EQ.MATNO) GO TO 80
      MTX = MT-1
      ND = 0
      DO 76 MM=MATNO,MTX
   76 ND = ND + NUM(MM)*(NOKEEP(NOLP,MM)-NOSKIP(NOLP,MM))
   80 CONTINUE
   85 DO 90 I1=1,ND1
      I = I1 + ND
      IF (I.GT.NDDD) GO TO 92
      DO 88 J1=I1,ND1
      J = J1 + ND
      IF (J1.GT.NDDD) GO TO 92
      J2 = J1 + NDA
      I2 = I1 +NDA
      IF (KC.NE.0) GO TO 86
      I2 = J1 + NDA + 1
      J2 = I1 + NDA
   86 CONTINUE
      IF (ABS(COEFFG(MT,NP2,I2,J2)).LE.0.0001) GO TO 88
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF)= COEFFG(MT,NP2,I2,J2)
      INDC(NCOEFF) = 1000*(I + 1000*NP) + J
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
   88 CONTINUE
   90 CONTINUE
      IF (NK.LT.NUM(MT)) GO TO 95
   91 MT = MT + 1
      NK = 1
      IF (MT.LE.NSUB) GO TO 75
      GO TO 70
   92 MMM = 92
      WRITE (6, 2001) MMM,NK,NUM(MT),MT,I1,J1
C      PRINT 2001, MMM,NK,NUM(MT),MT,I1,J1
 2001 FORMAT (1X,' ERROR IN XTLTP AT',I6,5I3)
      RETURN
   95 NK = NK + 1
      ND = ND + ND1
      GO TO 85
C
C      PUT IN CRYSTAL MATRIX ELEMENTS
C
  110 CONTINUE
      REWIND 8
  120 READ (8,END=495) MT,ND1,MTT,ND2,KU,((H(I,J),J=1,ND2),I=1,ND1)
      IF (NV1KOP.GT.0)
     1READ (8) MT,ND1,MTT,ND2,KV,((HZ(I,J),J=1,ND2),I=1,ND1)
  121 IF (MT.LT.MATNO) GO TO 120
      IF (NOKEEP(NOLP,MT).GT.0) GO TO 123
      IF (MT.GE.NSUB) GO TO  495
      GO TO 120
  123 IF (NOKEEP(NOLP,MTT).EQ.0) GO TO 120
      K2 = KU/2
      NP = NPAR + K2
      NPV = NPXT + K2
      J32 = 2*KU
      ND = 0
      NDA = NOSKIP(NOLP,MT)
      NDB = NOSKIP(NOLP,MTT)
      ND1X = NOKEEP(NOLP,MT)-NDA
      ND2X = NOKEEP(NOLP,MTT)-NDB
C
C   CYLINDRICAL CRYSTAL FIELD MATRIX ELEMENTS
      NDP = 0
      IF (MT.EQ.MATNO) GO TO 125
      MT1 = MT-1
      DO 122 MM=MATNO,MT1
  122 ND = ND +NUM(MM)*(NOKEEP(NOLP,MM)-NOSKIP(NOLP,MM))
  125 IF (MTT.EQ.MATNO) GO TO 130
      MT1 = MTT-1
      DO 126 MM=MATNO,MT1
  126 NDP=NDP +NUM(MM)*(NOKEEP(NOLP,MM)-NOSKIP(NOLP,MM))
  130 AJ1 = MT-1
      AJ1 = AJ1 + XADD
      J31 = 2.*AJ1
      AJ2 = MTT-1
      AJ2 = AJ2 + XADD
      J33 = 2.*AJ2
      ANO = 1.
      NDD = ND
      NDDP = NDP
      NK = 1
      XM = AMU
      MJ5 = 0
  135 MJ4 = 2.*XM*(-1.)
      NFAC = (J31+MJ4)/2
      MJ6 = 2.*XM
      RS = THRJ(J31,J32,J33,MJ4,MJ5,MJ6)
      IF (DABS(RS).LT.0.0001) GO TO 155
      RS = (-1.)**NFAC*RS
      DO 150 I1=1,ND1X
      I = I1 + NDD
      I2 = I1 + NDA
      IF (I.LE.NDDD) GO TO 136
      MMM = 150
  134 WRITE (6, 5005) MMM,KU,MT,MTT,I1,J1
C      PRINT 5005, MMM,KU,MT,MTT,I1,J1
      RETURN
  136 N1 = 1
      IF (MT.EQ.MTT) N1 = I1
      DO 140 J1=N1,ND2X
      J2 = J1 + NDB
      J = J1 + NDDP
      IF (J.LE.NDDD) GO TO 138
      MMM = 138
      GO TO 134
  138 IF (DABS(H(I2,J2)).LE.0.0001) GO TO 139
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF) = RS*H(I2,J2)*C(K2)
      INDC(NCOEFF) = 1000*(I + 1000*NP) + J
      IF (MT.EQ.1) WRITE (6, 1016) MT,I,MTT,J,NP,COEFF(NCOEFF)
C      IF (MT.EQ.1) PRINT 1016, MT,I,MTT,J,NP,COEFF(NCOEFF)
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
  139 IF (NV1KOP.EQ.0) GO TO 140
      IF (DABS(HZ(I2,J2)).LE.0.0001) GO TO 140
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF) = RS*HZ(I2,J2)*C(K2)
      INDC(NCOEFF) = 1000*(I + 1000*NPV) + J
      IF (MT.EQ.1) WRITE (6,1016) MT,I,MTT,J,NPV,COEFF(NCOEFF)
C      IF (MT.EQ.1) PRINT 1016, MT,I,MTT,J,NPV,COEFF(NCOEFF)
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
  140 CONTINUE
  150 CONTINUE
  155 CONTINUE
      IF (NK.EQ.NUM(MT) ) GO TO 160
            NKKX = NK/2
      NKXX2 = 2*NKKX
      IF (NK-NKXX2) 157,158,157
  157 XM = AMU-ANO*XMLP
      GO TO 159
  158 XM = AMU + ANO*XMLP
      ANO = ANO + 1
  159 NK = NK + 1
      NDDP = NDDP + ND2X
      NDD = NDD + ND1X
      GO TO 135
  160 CONTINUE
C
C      MATRIX ELEMENTS FOR Q.NE.0
      IF (KU.LT.NQ) GO TO 210
      NK = 0
      NKP = 0
      NO1 = 0
      NO2 = 0
  165 XM1 = NO1*NQ
      XM1 = XM1+AMU
      NDD =  ND + NK*(NOKEEP(NOLP,MT)-NDA)
      MJ1 = -2.*XM1
  170 XMP = NO2*NQ
      XMP = XMP + AMU
      MJ3 = 2*XMP
      NDDP = NDP + NKP*(NOKEEP(NOLP,MTT)-NDB)
      MJ2 =(MJ1 + MJ3)*(-1)
      IF (MJ2.EQ.0) GO TO 172
      MQ = MJ2/2
      NQ2 = IABS(MQ)
      IF (IABS(MJ2).LE.J32) GO TO 180
  172 NKP = NKP + 1
      IF ((NKP+1).GT.NUM(MTT)) GO TO 175
      IF (NO2) 174,173,173
  173 NO2=NO2 + 1
  174 NO2 = -NO2
      GO TO 170
  175 NK = NK + 1
      IF ((NK+1).GT.NUM(MT)) GO TO 210
      NKP = 0
      IF (MTT.EQ.MT) NKP = NK
      NO2 = 0
      IF (NO1) 178,176,176
  176 NO1 = NO1 + 1
  178 NO1 = -NO1
      IF (MT.EQ.MTT) NO2=NO1
      GO TO 165
  180 RS = THRJ(J31,J32,J33,MJ1,MJ2,MJ3)
      IF (DABS(RS).LT.0.0001) GO TO 172
      NPXD = NPAR+3
      NPXDV = NPXT + 3
      IF (K2.EQ.1) GO TO 182
      K21 = K2-1
      DO 181 I=1,K21
      NPXD = NPXD + KNQ(I)
      NPXDV = NPXDV + KNQ(I)
 181  CONTINUE
  182 NX = 1
C               PAR. NO CHANGE  6/83
C   B(KQ) AND B(K-Q) HAVE THE SAME PAR. NO
C
C          PHASE CORRECTION FOR Q=3
C
  183 NP1 = NPXD + NX
      NPV = NPXDV + NX
      IF (NQ2.EQ.NX*NQ) GO TO 184
      NX = NX + 1
      GO TO 183
  184 CONTINUE
      NPZ = (J31+MJ1)/2
      IF (MQ.GT.0) NPZ=NPZ+MQ
      RS = RS*(-1.)**NPZ
      DO 200 I1 = 1,ND1X
      I = I1 + NDD
      I2 = I1 + NDA
      DO 195 J1=1,ND2X
      J = J1 + NDDP
      J2 = J1 + NDB
      IF (DABS(H(I2,J2)).LT.0.0001) GO TO 194
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF) = RS*H(I2,J2)*C(K2)
      INDC(NCOEFF) = 1000*(I+1000*NP1) + J
      IF (MT.EQ.1) WRITE (6,1016) MT,I,MTT,J,NP1,COEFF(NCOEFF)
C      IF (MT.EQ.1) PRINT 1016, MT,I,MTT,J,NP1,COEFF(NCOEFF)
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
  194 IF (NV1KOP.EQ.0) GO TO 195
      IF (DABS(HZ(I2,J2)).LT.0.0001) GO TO 195
      NCOEFF = NCOEFF + 1
      COEFF(NCOEFF) = RS*HZ(I2,J2)*C(K2)
      INDC(NCOEFF) = 1000*(I+1000*NPV) + J
      IF (MT.EQ.1) WRITE (6, 1016) MT,I,MTT,J,NPV,COEFF(NCOEFF)
C      IF (MT.EQ.1) PRINT 1016, MT,I,MTT,J,NPV,COEFF(NCOEFF)
      IF (NCOEFF.EQ.4) CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
  195 CONTINUE
  200 CONTINUE
      GO TO 172
  210 IF (MT.LT.NSUB) GO TO 120
      IF (KU.LT.6) GO TO 120
  495 CONTINUE
      CALL WRT84(COEFF,INDC,NCOEFF,NNTT)
      WRITE (6, 6003) MTN,NNTT
  500 CONTINUE
C      PRINT 6003, MTN,NNTT
 6003 FORMAT (/1X,' MATRIX',I3,' HAS ',I12,' COEFFICIENTS'/)
  700 CONTINUE
      WRITE(6,9000)
      DATA ISE,ISV/' B  ','BV  '/
      NP = NPAR
      ISYM = ISE
  800 DO 850 I = 1,3
      K = 2*I
      MQ = 0
      NP = NP + 1
  850 WRITE(6,9001)NP,ISYM,K,MQ
      DO 860 I = 1,3
      K = 2*I
      LMI = KNQ(I)
      DO 860 J = 1,LMI
      MQ = NQ*J
      NP = NP + 1
  860 WRITE(6,9001) NP,ISYM,K,MQ
      ISYM = ISV
      IF(NV1KOP.NE.0.AND.NP.LE.NPXT) GOTO 800
 9000 FORMAT(5X,'FOR CRYSTAL FIELD PARAMETERS:'/
     & 5X,'INDEX PARAMETER')
 9001 FORMAT(5X,I3,7X,A2,2I1)
 1025 FORMAT (4I6,F12.8,1X,A3,2F5.1)
 1026 FORMAT (4I6,F12.8,1X,I3,2F5.1)
 1016 FORMAT (5I6,F12.8)
 5005 FORMAT (7I6,F10.8)
  101 FORMAT (1H )
 6001 FORMAT (6(I3,I2,I3,4X),2I3)
 6002 FORMAT (6(A3,I2,I3,4X),2I3)
      RETURN
      END
