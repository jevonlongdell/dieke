C
C    ==================================================
C                  FIT91.FOR
C       MODIFIED PROGRAM FROM FIT84.FOR FOR VAX/VMS SYSTEM
C	ADDITIONAL SUBROUTINS IN THE FILE: FIT91SUB.FOR
C       ADDITIONAL UNIT NEEDED 5 FIT91.DAT
C                              6 FIT91.OUT
C       	
C                 JUNE 10, 1991
C    ==================================================
C    ==================================================
C
C                  FIT84
C
C      LEAST SQUARES FITTING OF PARAMETERS TO ENERGY
C                        LEVELS
C               USES OUTPUT FROM XTAL84
C
C       LOWEST  LEVEL ADJUSTMENT TO 0 OPTIONAL
C
C                RANK 300-CRYSTAL VERSION
C          300 - MAXIMUM RANK (NEIGEN)
C          300 - MAXIMUM NUMBER OF LEVELS KEPT
C               PER SUBMATRIX (MAXKP)
C          40 - MAXIMUM NUMBER OF PARAMETERS (MPAR)
C          40 - MAXIMUM NUMBER OF VARIABLES (NPARAM)
C          6 - MAXIMUM NUMBER OF SUBMATRICES (MR)
C
C       INPUT:
C          1 - TITLE 2(10A8)
C          3-  6 NUMBERS EXPECTED
C              NUMBER OF SUBMATRICES (MAXSM)
C              NUMBER OF PARAMETERS (MAXPAR)
C              NUMBER OF ITERATIONS (MAXITR), IF = 0,
C                 NOT FITTED, NO LEVELS NEEDED
C              KEVEC GT 0 FOR VECTORS  UNIT 4 IS NEEDED
C              NZERO = 1 TO ZERO LOWEST LEVEL
C              NOLEV=NO. LEVELS TO PRINT FOR ITER=0
C               (IF 0, PRINTS THE SMALLER OF NDIM OR 100)
C          4 - PARAMETERS
C          9 - LEVELS IF MAXITR GT 0
C         10 - -1
C    ==================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 H( 300, 300),EVNAM( 300,6),DERIV( 300,40)
      REAL*8 TITLE(20),COEFFN(1000),NEWERR,DETDUM(2)
      REAL*8 LEVEL( 300,6),DIFF( 300),EVAL(300),E(300)
      REAL*8 A(40,40),ANAM(40,3),PARAM(40),PARLMP(40)
      REAL*8 B(40),BB(40),ZW(40),OLDPAR(40)
      REAL*4 RATIO(40),COEFFX(2000),BIGEV( 300),SECEV( 300)
      REAL*4 RMS(40)
      INTEGER MRANK(6),NEVAL(6),INDX(2000),I1EV( 300)
      INTEGER I2EV( 300),INERT(3),NEWNO(40),LUMPAR(40)
      INTEGER IPVT(40)
C
      COMMON/CWRT/COEFFX,INDX,NNT,NCOEFF,KKKDX
C
      NEIGEN =  300
      MAXKP =  300
      MR = 6
      MPAR = 40
      NPARAM = 40
      MLMP = 40
      MMM = 0
      READ (5,1000) (TITLE(K),K=1,10)
      READ (5,1000) (TITLE(K),K=11,20)
 1000 FORMAT (10A8)
      WRITE (6, 1003), (TITLE(K),K=1,20)
 1003 FORMAT (1H1,10A8/1H ,10A8)
      READ (5,*) MAXSM,MAXPAR,MAXITR,KEVEC,NZERO,NOLEV
 1001 FORMAT (24I3)
      WRITE (6, 1004) MAXSM,MAXPAR,MAXITR
 1004 FORMAT (1H0,I3,' SUBMATRICES, ',I6,' PARAMETERS, ',
     1' UP TO ',I2,' ITERATIONS'/)
C
C    ====================================================
C                  READ PARAMETERS
C    ====================================================
C
      CALL RDPAR(MPAR,NPARAM,MLMP,PARAM,PARLMP,
     1ANAM,RATIO,NEWNO,LUMPAR,NOVAR,MAXPAR,MAXLMP,MMM)
      IF (MMM.GT.0) GO TO 345
      NDERIV = NOVAR
      IF (MAXLMP.GT.0) NDERIV=NOVAR+1
C
C    ====================================================
C           READ STATE NAMES AND COEFFICIENTS
C    ====================================================
C
      CALL RDNMCF(NEIGEN,MPAR,MLMP,H,PARLMP,RATIO,NEWNO,
     1LUMPAR,MRANK,MMM,COEFFN,MR,MAXLMP,EVNAM,MAXSM)
      IF (MMM.NE.0) GO TO 345
C
C     ==================================================
C                   READ LEVELS
C     ===================================================
C
      CALL RDLEV4(NEVAL,LEVEL,DEG,NOVAR,MAXKP,
     1MRANK,INDLEV,MR,MAXSM,MAXITR,NOLEV)
      IF (DEG.GT.0.) GO TO 45
      NUMITR =  0
      MAXITR = 0
      GO TO 58
   45 NUMITR = 0
      NEWERR = 0.0
      Q = NOVAR
      Q = Q/DEG
      OLDRES = 0.0
   50 CONTINUE
      NUMITR = NUMITR + 1
      IF (NUMITR.GT.MAXITR) GO TO 550
      IF (KEVEC.GT.0) REWIND 4
C
C     =================================================
C           ZERO DERIVATIVE MATRIX
C     ============================================
C
      RESID = 0.0
      DO 55 NP=1,NDERIV
      DO 54 NP1=1,NDERIV
   54 A(NP,NP1) = 0.0
      B(NP) = 0.0
   55 BB(NP) = 0.0
   58 CONTINUE

C
C     =================================================
C                MAIN LOOP
C     =================================================
C
      ADJ = 1000000.
      WRITE (6, 1003) (TITLE(I),I=1,20)
      IF (NZERO.NE.0) GO TO 59
      WRITE (6, 1008) NUMITR
 1008 FORMAT (//,25X,'ITERATION NO.',I3/)
      WRITE (6, 1009)
 1009 FORMAT (1H0,'  M  N  OBS. LEVEL(N)   EIGENVALUE(N)',
     1'   T-W    LARGEST EV.COMP.   NEXT LGST  EV     '/1H )
   59 DO 400 M=1,MAXSM
      WRITE (6, 101)
  101 FORMAT (1H )
      NEVALM = NEVAL(M)
      IF (NEVALM.EQ.0) GO TO 400
C     IF NO LEVELS AND MAXITR.GT.0 , OMIT THE SUBMATRIX
C
      MDIM  = MRANK(M)
C
C     =================================================
C              FORM ENERGY MATRIX
C     =================================================
C
      CALL CALC(NEIGEN,NPARAM,COEFFN,PARAM,H,MDIM,MAXLMP)
      IF (MDIM.NE.1) GO TO 142
      EVAL(1) = H(1,1)
      H(1,1) = 1.0D0
      H(2,1) = 0.0
      GO TO 145
C
C      ================================================
C                        DIAGONALIZE
C             CAN KEEP ALL VECTORS
C      ================================================
C
  142 NVPT = 0
      CALL VECTR2(H,EVAL,MDIM,NEVALM,NVPT,MMM,NEIGEN,E)
      IF (MMM.EQ.0) GO TO 145
      MMM = 142
      GO TO 345
  145 CONTINUE
      IF (KEVEC.NE.0) WRITE (4) ((H(K,I),K=1,MDIM),
     1I=1,NEVALM),(EVAL(J),J=1,NEVALM)
C
C      FIND TWO LARGEST VECTOR COMPONENTS
C
      CALL SORTVC(H,MDIM,NEVALM,I1EV,I2EV,NEIGEN,MAXKP,
     1BIGEV,SECEV)
C
      RSUM = 0.0
      DO 65 I=1,NEVALM
      DIFF(I) = 1.
      IF (LEVEL(I,M).EQ.1.0) GO TO 65
      DIFF(I) = LEVEL(I,M)-EVAL(I)
      RSUM = DIFF(I)**2 + RSUM
   65 CONTINUE
      RESID = RESID + RSUM
      IF (NZERO.EQ.0) GO TO 68
      IF (EVAL(1).LT.ADJ) ADJ = EVAL(1)
      WRITE (19) (EVAL(I),I1EV(I),BIGEV(I),I2EV(I),SECEV(I),
     1I=1,NEVALM)
C
C     ==================================================
C             ADD TO A MATRIX
C     ==================================================
C
   68 CONTINUE
      IF (MAXITR.NE.0) CALL ABCAL4(M,NEIGEN,NPARAM,MAXKP,
     1NOVAR ,NEVALM,DERIV,H,A,B,BB,DIFF,LEVEL,MR)
      IF (NZERO.NE.0) GO TO 400
      DO 70 I=1,NEVALM
      J1EV = I1EV(I)
      J2EV = I2EV(I)
      WRITE (6, 1010)M,I,LEVEL(I,M),EVAL(I),DIFF(I),I1EV(I),
     1BIGEV(I),EVNAM(J1EV,M),J2EV,SECEV(I),EVNAM(J2EV,M)
   70 CONTINUE
  400 CONTINUE
      IF (NZERO.EQ.0) GO TO 72
C
      REWIND 19
      WRITE (6, 1018) NUMITR, ADJ
 1018 FORMAT (//50X,'ITERATION NO.',I3/5X,
     1'LOWEST LEVEL MADE 0.0 IN PRINTED OUTPUT BY',
     21X,'SUBTRACTING ',F12.3//)
      WRITE (6, 1009)
      DO 85 M=1,MAXSM
      NEVALM = NEVAL(M)
      READ (19)(EVAL(I),I1EV(I),BIGEV(I),I2EV(I),SECEV(I),
     1I=1,NEVALM)
      DO 170 I=1,NEVALM
      J1EV = I1EV(I)
      J2EV = I2EV(I)
      CALCEV = EVAL(I) -ADJ
      DIFFX = 1.
      IF (LEVEL(I,M).NE.1.) DIFFX = LEVEL(I,M)-CALCEV
      WRITE (6, 1010)M,I,LEVEL(I,M),CALCEV,DIFFX,I1EV(I),
     1BIGEV(I),EVNAM(J1EV,M),J2EV,SECEV(I),EVNAM(J2EV,M)
  170 CONTINUE
      WRITE (6, 101)
   85 CONTINUE
      REWIND 19
   72 IF (MAXITR.EQ.0) GO TO 470
      REWIND 2
      REWIND 31
      IF (MAXLMP.GT.0) A(NDERIV,NDERIV) = 1.D0
C
C     ================================================
C              CALCULATE PARAMETER CHANGES
C    =================================================
C
      SIGMA = DSQRT(RESID/DEG)
      CALL DSICO(A,NPARAM,NDERIV,IPVT,RCOND,ZW)
      T = 1.0 + RCOND
      IF (T.NE.1.0) GO TO 405
      MMM = 405
      WRITE (6, 1015) T
      GO TO 345
 1010 FORMAT (1X,2I3,F12.3,F12.3,F9.3 ,2(1X,I3,1H(,
     1F6.3,1H),A8))
 1015 FORMAT (1X,' A-MATRIX IS ILL-CONDITIONED, T=',
     1E13.6)
  405 CALL DSISL(A,NPARAM,NDERIV,IPVT,B)
      CALL DSIDI(A,NPARAM,NDERIV,IPVT,DETDUM,INERT,ZW,1)
 1050 FORMAT (10E12.4)
      OLDERR = NEWERR
      DRESID = OLDRES-RESID
      OLDRES = RESID
      NEWERR = 0.0001 * OLDRES
      DO 410 NP=1,NOVAR
      OLDPAR(NP) = PARAM(NP)
      RMS(NP) = SIGMA*DSQRT(DABS(A(NP,NP )))
      PARAM(NP) = OLDPAR(NP) + B(NP)
  410 CONTINUE
  450 WRITE (6, 1020)NUMITR
 1020 FORMAT (///30X,'ITERATION NO. ',I3//)
      IF (NUMITR.EQ.1) GO TO 460
      IF( DABS(DRESID).LE.NEWERR) WRITE (6, 1021)
 1021 FORMAT (32X,'CONVERGED')
      IF (NEWERR.GT.OLDERR) WRITE (6, 1022)
 1022 FORMAT (32X,' DIVERGED'/)
  460 WRITE (6, 1025)OLDRES,DRESID,SIGMA,NOVAR,INDLEV,DEG
 1025 FORMAT (1H0,'RESIDUAL=',E13.3,1X,'CHANGE IN RESIDUAL ='
     1,E13.3,/1H ,12X,' SIGMA=',E13.3,/1H ,'NO. OF INDEPENDENT ',
     2'PARAMETERS =',I4,' NO. OF LEVELS FITTED =',I4,/
     3'   DEGREES OF FREEDOM =', F5.1)
  470 CONTINUE
      WRITE (6, 1028)
 1028 FORMAT (1H0,'  P  OLD PARAMETER   ERROR(RMS) NAME',
     114X,' CHANGE   NEW PARAMETER DR/DP'/)
      DO 480 NP=1,NOVAR
      IF (MAXITR.EQ.0) GO TO 475
      WRITE (6, 1030)NP,OLDPAR(NP),RMS(NP),(ANAM(NP,K),K=1,3),
     1B(NP),PARAM(NP),BB(NP)
      GO TO 480
  475 WRITE (6, 1032)NP,PARAM(NP),(ANAM(NP,K1),K1=1,3)
  480 CONTINUE
 1032 FORMAT (1H ,I4,3X,F12.3,1X,3A6)
 1030 FORMAT (1H ,I4,3X,2F12.3,1X,3A6,F9.3,F12.3,3F9.3)
      IF (MAXLMP.EQ.0) GO TO 492
      WRITE (6, 1034)
 1034 FORMAT (1H0,' LUMPED PARAMETERS')
      DO 490 I=1,MAXLMP
  490 WRITE (6, 1032)LUMPAR(I),PARLMP(I)
  492 CONTINUE
      IF (KEVEC.EQ.0) GO TO 500
C
C     =================================================
C             PRINT ALL VECTOR COMPONENTS
C      ================================================
C
      CALL PRTV4(MAXSM,NEVAL,MRANK,NEIGEN,H,EVAL,EVNAM,MR,ADJ)
  500 CONTINUE
      IF (MAXITR.EQ.0) GO TO 550
      IF (DABS(DRESID).LT.NEWERR) GO TO 550
      GO TO 50
  345 WRITE (6, 1345) MMM
 1345 FORMAT (1X, 'ERROR AT ',I6)
  550  REWIND 16
      STOP
      END
C
C ------
C
      SUBROUTINE RDPAR(MPAR,NPARAM,MLMP,PARAM,PARLMP,
     1ANAM,RATIO,NEWNO,LUMPAR,NDERIV,MAXPAR,MAXLMP,MMM)
C
C     ===============================================
C          READ PARAMETERS--JULY 1982
C         IF TOO MANY LUMPED PARAMETERS, THE
C           EXTRA ARE VARIED
C         IF TOO MANY VARIABLES, THE EXTRA
C            ARE FIXED
C         IF PARAMETER IS IN RATIO TO A FIXED
C           PARAMETER IT IS FIXED
C        RATIO MAY BE READ IN IN COL 10-15 AS WELL
C            AS IN COLUMN 46-54
C
C    ================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PARAM(NPARAM),PARLMP(MLMP),ANAM(NPARAM,3)
      REAL*4 RATIO(MPAR),RTO,RT1
      INTEGER NEWNO(MPAR),LUMPAR(MLMP)
      I = 1
      J = 1
      MAXLMP = 0
      DO 20 IP=1,MAXPAR
      READ (5,1200) NP,IVAR,NP1,RT1,PARM,TNAM1,TNAM2,TNAM3,RTO
      WRITE (6, 1201) NP,IVAR,NP1,RT1, PARM,TNAM1,TNAM2,TNAM3,RTO
 1200 FORMAT (3I3,F6.3,F11.3,1X,3A6,F8.5)
 1201 FORMAT (1X,3I3,F6.3,F11.3,1X,3A6,F8.5)
      RATIO(NP) = 1.
      IF (IVAR.EQ.0) GO TO 18
      IF (IVAR.EQ.2) GO TO 15
      IF (IVAR.EQ.1) GO TO 10
      IVAR = 0
      GO TO 18
   10 IF (I.GT.MPAR) GO TO 18
      NEWNO(NP) = I
      PARAM(I) = PARM
      ANAM(I,1) = TNAM1
      ANAM(I,2) = TNAM2
      ANAM(I,3) = TNAM3
      I = I + 1
      GO TO 20
   15 NEWNO(NP) = NEWNO(NP1)
      IF (NEWNO(NP1).EQ.0) GO TO 18
      RATIO(NP) = RT1
      IF (RT1.EQ.0.) RATIO(NP) = RTO
      IF (RATIO(NP).EQ.0.) RATIO(NP) = 1.
      GO TO 20
   18 IF (J.LE.MLMP) GO TO 19
      IF (I.LE.MPAR) GO TO 10
      WRITE (6, 1202)
 1202 FORMAT (1X,' TOO MANY PARAMETERS')
      MMM = 200
      RETURN
   19 CONTINUE
      LUMPAR(J) = NP
      PARLMP(J) = PARM
      NEWNO(NP) = 0
      J = J + 1
   20 CONTINUE
      MAXLMP = J-1
      NDERIV = I-1
      WRITE (6, 1205)  MAXPAR,NDERIV
 1205 FORMAT (1H0,I6,' PARAMETERS ARE READ',I6,' ARE FREE')
      WRITE (6, 1210)
 1210 FORMAT (1H0,1X,' VARIABLE PARAMETERS'/1X,
     1' NO. NEWNO PARAMETER   RATIO')
      DO 35 I=1,MAXPAR
      IF (NEWNO(I).EQ.0) GO TO 35
      J = NEWNO(I)
      WRITE (6, 1215) I,NEWNO(I),PARAM(J),RATIO(I)
 1215 FORMAT (1X,I3,I5,F12.3,F6.3)
   35 CONTINUE
      IF (MAXLMP.EQ.0) GO TO 45
      WRITE (6, 1220)
 1220 FORMAT (1H0,1X,' LUMPED PARAMETERS')
      DO 40 I=1,MAXLMP
      WRITE (6, 1225)  LUMPAR(I),I,PARLMP(I)
 1225 FORMAT (1X,2I3,F12.3)
   40 CONTINUE
   45 CONTINUE
      RETURN
      END
C
C ------
C
      SUBROUTINE RDNMCF(NEIGEN,MPAR,MLMP,H,PARLMP,RATIO,
     1NEWNO,LUMPAR,MRANK,MMM,COEFFN,MR,MAXLMP,EVNAM,MAXSM)
C
C     ====================================================
C                 VERSION RD16
C
C           READ STATE NAMES
C           READ COEFFICIENTS, LUMP AS NOTED
C           RD16 TO READ COEFF IN 4(I9,F11.6) FORMAT
C           WRITE ON 2 FOR ENERGY MATRIX
C            WRITE ON 31 FOR DERIVATIVES
C     ====================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 H(NEIGEN,NEIGEN),PARLMP(MLMP),COEFFN(1000)
      REAL*4 RATIO(MPAR),COEFFX(2000)
      INTEGER INDX(2000),NEWNO(MPAR),MRANK(MR)
      DIMENSION LUMPAR(MLMP),EVNAM(NEIGEN,MR),AME(4),NAJ(4)
      DIMENSION NAI(4),NANP(4)
      COMMON/CWRT/COEFFX,INDX,NNT,NCOEFF,KKKDX
      NNT = 0
      DO 200 M=1,MAXSM
      NUMBER = 0
      NCOEFF = 0
      WRITE (6, 1010)
 1010 FORMAT (1H0,60(1H-))
      KKKDX = 1
       READ (16,1030) M1,MDIM
 1030 FORMAT (3I8)
      MRANK(M) = MDIM
      DO 20 I=1,MDIM
      DO 18 J=I,MDIM
   18 H(I,J) = 0.0
   20 CONTINUE
C
C  READ STATE NAMES
C
      READ (16,1110) (EVNAM(I,M),I=1,MDIM)
 1110 FORMAT (6(A8,4X))
   22 READ (16 ,1035)(NANP(LX),NAI(LX),NAJ(LX),AME(LX),LX=1,4)
 1035 FORMAT (4(3I3,F11.6))
      DO 48 LX=1,4
      NP = NANP(LX)
      IF (NP.EQ.0) GO TO 60
      I = NAI(LX)
      J = NAJ(LX)
      IF (NEWNO(NP).GT.0) GO TO 45
      DO 42 K1=1,MAXLMP
      IF (NP.EQ.LUMPAR(K1)) GO TO 44
   42 CONTINUE
      MMM = 2042
      WRITE (6, 1024) NP
 1024 FORMAT (1X,'ERROR IN COEFFICIENT ',I6)
      RETURN
   44 H(I,J) = H(I,J) + PARLMP(K1)*AME(LX)
      GO TO 48
   45 NCOEFF = NCOEFF + 1
      COEFFX(NCOEFF) = AME(LX)*RATIO(NP)
      INDX(NCOEFF) = (I + 1000*NEWNO(NP))*1000+J
      IF (NCOEFF.EQ.2000) CALL WRIT1(M,MDIM)
   48 CONTINUE
      GO TO 22
   60 KKKDX = 0
      CALL WRIT1(M,MDIM)
C  PUT LUMPED VALUES ON UNIT 2
      NTX = 0
      IF (MAXLMP.EQ.0) GO TO 199
      KKKDX = 1
      DO 70 I=1,MDIM
      DO 68 J=I,MDIM
      IF (DABS(H(I,J)).LT.0.00001) GO TO 68
      NCOEFF = NCOEFF + 1
      COEFFN(NCOEFF) = H(I,J)
      INDX(NCOEFF) = J + 1000*I
      IF (NCOEFF.LT.1000) GO TO 68
      WRITE (2) M,MDIM,NCOEFF,KKKDX,(COEFFN(LX),INDX(LX),
     1LX=1,NCOEFF)
      NTX = NTX + NCOEFF
      NCOEFF = 0
   68 CONTINUE
   70 CONTINUE
      IF (NCOEFF.NE.0) GO TO 75
      NCOEFF = 1
      COEFFN(1) = 0.0
      INDX(1) = 1001*MDIM
   75 CONTINUE
      KKKDX = 0
      WRITE (2) M,MDIM,NCOEFF,KKKDX,(COEFFN(LX),INDX(LX),
     1LX=1,NCOEFF)
  199 CONTINUE
      NTX = NTX + NCOEFF
      WRITE (6, 1015) M,MDIM,NNT,NTX
 1015 FORMAT (/1X,' MATRIX ',I3,' (RANK ',I4,'), VARIABLE COEFF=',
     1I8,' CONSTANT COEFF=',I8/)
  200 CONTINUE
      WRITE (6, 2000)
 2000 FORMAT (1H0,' COEFFICIENTS READ')
      REWIND 2
      REWIND 31
      RETURN
      END
C
C ------
C
      SUBROUTINE RDLEV4(NEVAL,LEVEL,DEG,NDERIV,MAXKP,
     1MRANK,INDLEV,MR,MAXSM,MAXITR,NOLEV)
C
C    ======================================================
C                READ LEVELS
C    ======================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LEVEL(MAXKP,MR),LEV
      INTEGER NEVAL(MR),MRANK(MR)
      DEG = 0.0
      INDLEV = 0
      IF (NOLEV.EQ.0) NOLEV=100
      IF (MAXITR.EQ.0) GO TO 40
      DO 20 M=1,MAXSM
   20 NEVAL(M) = 0
      M = 0
   25 READ (5,1006) M1,LEV
 1006 FORMAT (I3,3X,F12.1,4X,2A8)
      IF (M1.LE.-1) GO TO 30
C  NOTE CHANGE TO PERMIT INCORRECT SUBMATRIX NUMBER IF IT IS
C     LESS THAN PRESENT
      IF (M1.LE.M) GO TO 28
      M = M1
      I = 0
   28 NEVAL(M) = NEVAL(M) + 1
      I = I + 1
      LEVEL(I,M) = LEV
      GO TO 25
   30 DO 35 M=1,MAXSM
      NEVALM = NEVAL(M)
      IF (NEVALM.EQ.0) GO TO 35
      INDLEV = INDLEV + NEVALM
      DO 34 I=1,NEVALM
      IF (LEVEL(I,M).EQ.1.0) INDLEV = INDLEV-1
   34 CONTINUE
   35 CONTINUE
      DEG = INDLEV-NDERIV
      IF (DEG.LE.0) MAXITR = 0
      RETURN
   40 DO 45 M=1,MAXSM
      NEVAL(M) = MRANK(M)
      IF (MRANK(M).GT.NOLEV) NEVAL(M) =NOLEV
      ND = NEVAL(M)
      DO 42 I=1,ND
      LEVEL(I,M) = 1.
   42 CONTINUE
   45 CONTINUE
      RETURN
      END
C
C ------
C
      SUBROUTINE CALC(NEIGEN,NPARAM,COEFFN,PARAM,H,
     1MDIM,MAXLMP)
C
C     ================================================
C         FORM ENERGY MATRIX
C     ================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 H(NEIGEN,NEIGEN),PARAM(NPARAM),COEFFN(1000)
      REAL*4 COEFFX(2000)
      INTEGER INDX(2000)
      COMMON/CWRT/COEFFX,INDX,NNT,NCOEFF,KKKDX
      DO 20 I=1,MDIM
      DO 19 J=1,I
   19 H(I,J) = 0.0
   20 CONTINUE
   25 READ (2) M1,MD,NCOEFF,KKKD,(COEFFX(LX),INDX(LX),
     1LX=1,NCOEFF)
      DO 30 IAT=1,NCOEFF
      NP = INDX(IAT)/1000000
      K1 = INDX(IAT)-1000000*NP
      J = K1/1000
      I = K1-1000*J
      H(I,J) = H(I,J) + COEFFX(IAT)*PARAM(NP)
   30 CONTINUE
      IF (KKKD.GT.0) GO TO 25
C  PUT IN CONSTANT TERMS
      IF (MAXLMP.EQ.0) GO TO 36
   33 READ (2) M1,MD,NCOEFF,KKKD,(COEFFN(LX),INDX(LX),
     1LX=1,NCOEFF)
      DO 35 LX=1,NCOEFF
      J = INDX(LX)/1000
      I = INDX(LX)-1000*J
      H(I,J) = H(I,J) + COEFFN(LX)
   35 CONTINUE
      IF (KKKD.NE.0) GO TO 33
   36 CONTINUE
      RETURN
      END
C
C ------
C
      SUBROUTINE VECTR2(H,EVAL,ND,NKP,NVPT,MMM,NEIGEN,E)
C
C     ==================================================
C          CALCULATE EIGENVECTORS AND EIGENVALUES
C            ALL VECTORS  CALCULATED
C               ND = DIMENSION
C               NKP=NUMBER DESIRED TO KEEP
C               NVPT=1 FOR PRINTING LARGEST VECTOR
C                      COMPONENT
C     ===================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NEIGEN,NEIGEN),EVAL(NEIGEN),E(NEIGEN)
      CALL TRED2(NEIGEN,ND,H,EVAL,E,H)
      CALL TQL2(NEIGEN,ND,EVAL,E,H,IERR)
      IF (IERR.EQ.0) GO TO 10
    5 WRITE (6, 1012) IERR
      MMM = 5
      RETURN
   10 IF (NVPT.EQ.0) GO TO 30
      WRITE (6, 1010)
      DO 25 I=1,NKP
      MM = 1
      TM = DABS(H(1,I))
      DO 20 J=2,ND
      IF ((DABS(H(J,I))-TM).LE.0.01) GO TO 20
      MM = J
      TM = DABS(H(J,I))
   20 CONTINUE
      WRITE (6, 1020) I,EVAL(I),MM,H(MM,I)
   25 CONTINUE
      WRITE (6, 101)
   30 CONTINUE
  101 FORMAT (1H )
 1010 FORMAT (1X,' VECTORS CALCULAATED')
 1012 FORMAT (1X,' ERROR IN TQL2, CODE = ',I5)
 1020 FORMAT (I6,F12.2,I6,F12.8)
      RETURN
      END
C
C ------
C
      SUBROUTINE SORTVC(H,MDIM,NEVALM,I1EV,I2EV,
     1NEIGEN,MAXKP,BIGEV,SECEV)
C
C     ===============================================
C              FIND TWO LARGEST VECTOR
C                 COMPONENTS
C     ===============================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(NEIGEN,NEIGEN),I1EV(MAXKP),I2EV(MAXKP)
      REAL*4 BIGEV(MAXKP),SECEV(MAXKP)
      DO 150 I=1,NEVALM
      BIGGST = 0.0
      BIGGA = DABS(BIGGST)
      SECOND = 0.0
      SECA = DABS(SECOND)
      I1EV(I) = 0
      DO 148 K=1,MDIM
      TESTEV = DABS(H(K,I))
      IF (BIGGA.GE.TESTEV) GO TO 146
      SECOND = BIGGST
      SECA = DABS(SECOND)
      I2EV(I) = I1EV(I)
      BIGGST = H(K,I)
      BIGGA = DABS(BIGGST)
      I1EV(I) = K
      GO TO 148
  146 IF (SECA.GE.TESTEV) GO TO 148
      SECOND = H(K,I)
      SECA = DABS(SECOND)
      I2EV(I) = K
  148 CONTINUE
      BIGEV(I) = BIGGST
      SECEV(I) = SECOND
  150 CONTINUE
      RETURN
      END
C
C ------
C
      SUBROUTINE ABCAL4(M,NEIGEN,NPARAM,MAXKP,NDERIV,NEVALM,
     1DERIV,Z,A,B,BB,DIFF,LEVEL,MR)
C
C     ======================================================
C          ACCUMULATE DERIV**2 IN 'A' MATRIX
C          ACCUMULATE B AND BB VECTORS
C     ======================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LEVEL(MAXKP,MR),DERIV(NEIGEN,NPARAM),Z(NEIGEN,
     1NEIGEN),A(NPARAM,NPARAM),B(NPARAM),BB(NPARAM),
     2DIFF(MAXKP)
      REAL*4 COEFFX(2000)
      INTEGER INDX(2000)
      COMMON/CWRT/COEFFX,INDX,NNT,NCOEFF,KKKDX
      DO 20 I=1,NEVALM
      DO 18 J=1,NDERIV
   18 DERIV(I,J) = 0.0
   20 CONTINUE
   75 READ (31) M1,MD,NCOEFF,KKKD,(COEFFX(LX),
     1INDX(LX),LX=1,NCOEFF)
      DO 78 LX=1,NCOEFF
      NP = INDX(LX)/1000000
      K1 = INDX(LX) -1000000*NP
      I = K1/1000
      K = K1-1000*I
      FAC = 1.
      IF (I.NE.K) FAC = 2.
      DO 76 J=1,NEVALM
   76 DERIV (J,NP) = DERIV(J,NP) + Z(I,J)*Z(K,J)*
     1COEFFX(LX)*FAC
   78 CONTINUE
      IF (KKKD.NE.0) GO TO 75
      DO 82 NP1=1,NDERIV
      DO 85 NP=1,NDERIV
      ATEMP = 0.0
      DO 80 J=1,NEVALM
      IF (LEVEL(J,M).EQ.1.0) GOTO 80
      ATEMP = ATEMP + DERIV(J,NP)*DERIV(J,NP1)
   80 CONTINUE
      A(NP,NP1) = A(NP,NP1) + ATEMP
   85 CONTINUE
   82 CONTINUE
      DO 90 NP=1,NDERIV
      BTEMP = 0.0
      DO 88 J=1,NEVALM
      IF (LEVEL(J,M).EQ.1.0) GO TO 88
      BTEMP = DERIV(J,NP)*DIFF(J) + BTEMP
   88 CONTINUE
      B(NP) = B(NP) + BTEMP
      BB(NP) = B(NP)
   90 CONTINUE
      RETURN
      END
C
C ------
C
      SUBROUTINE WRIT1(M,MDIM)
C
C     ===================================================
C              WRITE MATRIX ELEMENTS ON UNITS 2 AND 31
C     ===================================================
C
      REAL*4 COEFFX(2000)
      DIMENSION INDX(2000)
      COMMON/CWRT/COEFFX,INDX,NNT,NCOEFF,KKKD
      NNT = NNT + NCOEFF
      WRITE (2) M,MDIM,NCOEFF,KKKD,(COEFFX(LX),INDX(LX),
     1LX=1,NCOEFF)
      WRITE (31) M,MDIM,NCOEFF,KKKD,(COEFFX(LX),INDX(LX),
     1LX=1,NCOEFF)
      NCOEFF = 0
      RETURN
      END
C
C ------
C
      SUBROUTINE PRTV4(MAXSM,NEVAL,MRANK,NEIGEN,Z,EVAL,
     1EVNAM,MR,ADJ)
C
C     ====================================================
C              PRINT COMPLETE VECTORS
C     ====================================================
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(NEIGEN,NEIGEN),EVAL(NEIGEN),EVNAM(NEIGEN,MR)
      DIMENSION MRANK(MR),NEVAL(MR)
      REWIND 4
      DO 50 M=1,MAXSM
      WRITE (6, 1005) M
      NEVALM=NEVAL(M)
      MDIM = MRANK(M)
      READ (4) ((Z(K,I),K=1,MDIM),I=1,NEVALM),(EVAL(I1),
     1I1=1,NEVALM)
      DO 10 I=1,NEVALM
   10 EVAL(I) = EVAL(I) -ADJ
      DO 40 K=1,NEVALM,10
      KK = K + 9
      IF (KK.GT.NEVALM) KK = NEVALM
      WRITE (6, 1006) (EVAL(I),I=K,KK)
      DO 35 J=1,MDIM
      WRITE (6, 1007) EVNAM(J,M),(Z(J,I),I=K,KK)
   35 CONTINUE
   40 CONTINUE
   50 CONTINUE
      REWIND 4
 1005 FORMAT (1H0,132(1H-)/1H0,40X,'EIGENVECTORS, SUBMATRIX',
     1I3/)
 1006 FORMAT (/11X,10(1X,F9.2,1X)/)
 1007 FORMAT (1X,A8,1X,10F11.7)
      RETURN
      END
