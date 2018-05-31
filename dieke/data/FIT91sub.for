      SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)                            00000010
      INTEGER LDA,N,KPVT(1)                                             00000020
      DOUBLE PRECISION A(LDA,1),Z(1)                                    00000030
      DOUBLE PRECISION RCOND                                            00000040
C                                                                       00000050
C     DSICO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION  00000060
C     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE        00000070
C     MATRIX.                                                           00000080
C                                                                       00000090
C     IF  RCOND  IS NOT NEEDED, DSIFA IS SLIGHTLY FASTER.               00000100
C     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL.                        00000110
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL.                 00000120
C     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI.                   00000130
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI.               00000140
C     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI.                    00000150
C                                                                       00000160
C     ON ENTRY                                                          00000170
C                                                                       00000180
C        A       DOUBLE PRECISION(LDA, N)                               00000190
C                THE SYMMETRIC MATRIX TO BE FACTORED.                   00000200
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         00000210
C                                                                       00000220
C        LDA     INTEGER                                                00000230
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000240
C                                                                       00000250
C        N       INTEGER                                                00000260
C                THE ORDER OF THE MATRIX  A .                           00000270
C                                                                       00000280
C     OUTPUT                                                            00000290
C                                                                       00000300
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      00000310
C                WERE USED TO OBTAIN IT.                                00000320
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     00000330
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         00000340
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            00000350
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            00000360
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         00000370
C                                                                       00000380
C        KPVT    INTEGER(N)                                             00000390
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000400
C                                                                       00000410
C        RCOND   DOUBLE PRECISION                                       00000420
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        00000430
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       00000440
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             00000450
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 00000460
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     00000470
C                           1.0 + RCOND .EQ. 1.0                        00000480
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           00000490
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         00000500
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          00000510
C                UNDERFLOWS.                                            00000520
C                                                                       00000530
C        Z       DOUBLE PRECISION(N)                                    00000540
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  00000550
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      00000560
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           00000570
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    00000580
C                                                                       00000590
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000600
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000610
C                                                                       00000620
C     SUBROUTINES AND FUNCTIONS                                         00000630
C                                                                       00000640
C     LINPACK DSIFA                                                     00000650
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       00000660
C     FORTRAN DABS,DMAX1,IABS,DSIGN                                     00000670
C                                                                       00000680
C     INTERNAL VARIABLES                                                00000690
C                                                                       00000700
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T                  00000710
      DOUBLE PRECISION ANORM,S,DASUM,YNORM                              00000720
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS                                  00000730
C                                                                       00000740
C                                                                       00000750
C     FIND NORM OF A USING ONLY UPPER HALF                              00000760
C                                                                       00000770
      DO 30 J = 1, N                                                    00000780
         Z(J) = DASUM(J,A(1,J),1)                                       00000790
         JM1 = J - 1                                                    00000800
         IF (JM1 .LT. 1) GO TO 20                                       00000810
         DO 10 I = 1, JM1                                               00000820
            Z(I) = Z(I) + DABS(A(I,J))                                  00000830
   10    CONTINUE                                                       00000840
   20    CONTINUE                                                       00000850
   30 CONTINUE                                                          00000860
      ANORM = 0.0D0                                                     00000870
      DO 40 J = 1, N                                                    00000880
         ANORM = DMAX1(ANORM,Z(J))                                      00000890
   40 CONTINUE                                                          00000900
C                                                                       00000910
C     FACTOR                                                            00000920
C                                                                       00000930
      CALL DSIFA(A,LDA,N,KPVT,INFO)                                     00000940
C                                                                       00000950
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              00000960
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .         00000970
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           00000980
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .                   00000990
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            00001000
C                                                                       00001010
C     SOLVE U*D*W = E                                                   00001020
C                                                                       00001030
      EK = 1.0D0                                                        00001040
      DO 50 J = 1, N                                                    00001050
         Z(J) = 0.0D0                                                   00001060
   50 CONTINUE                                                          00001070
      K = N                                                             00001080
   60 IF (K .EQ. 0) GO TO 120                                           00001090
         KS = 1                                                         00001100
         IF (KPVT(K) .LT. 0) KS = 2                                     00001110
         KP = IABS(KPVT(K))                                             00001120
         KPS = K + 1 - KS                                               00001130
         IF (KP .EQ. KPS) GO TO 70                                      00001140
            T = Z(KPS)                                                  00001150
            Z(KPS) = Z(KP)                                              00001160
            Z(KP) = T                                                   00001170
   70    CONTINUE                                                       00001180
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))                       00001190
         Z(K) = Z(K) + EK                                               00001200
         CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                          00001210
         IF (KS .EQ. 1) GO TO 80                                        00001220
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))                00001230
            Z(K-1) = Z(K-1) + EK                                        00001240
            CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)                   00001250
   80    CONTINUE                                                       00001260
         IF (KS .EQ. 2) GO TO 100                                       00001270
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 90                  00001280
               S = DABS(A(K,K))/DABS(Z(K))                              00001290
               CALL DSCAL(N,S,Z,1)                                      00001300
               EK = S*EK                                                00001310
   90       CONTINUE                                                    00001320
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   00001330
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         00001340
         GO TO 110                                                      00001350
  100    CONTINUE                                                       00001360
            AK = A(K,K)/A(K-1,K)                                        00001370
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  00001380
            BK = Z(K)/A(K-1,K)                                          00001390
            BKM1 = Z(K-1)/A(K-1,K)                                      00001400
            DENOM = AK*AKM1 - 1.0D0                                     00001410
            Z(K) = (AKM1*BK - BKM1)/DENOM                               00001420
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               00001430
  110    CONTINUE                                                       00001440
         K = K - KS                                                     00001450
      GO TO 60                                                          00001460
  120 CONTINUE                                                          00001470
      S = 1.0D0/DASUM(N,Z,1)                                            00001480
      CALL DSCAL(N,S,Z,1)                                               00001490
C                                                                       00001500
C     SOLVE TRANS(U)*Y = W                                              00001510
C                                                                       00001520
      K = 1                                                             00001530
  130 IF (K .GT. N) GO TO 160                                           00001540
         KS = 1                                                         00001550
         IF (KPVT(K) .LT. 0) KS = 2                                     00001560
         IF (K .EQ. 1) GO TO 150                                        00001570
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)                     00001580
            IF (KS .EQ. 2)                                              00001590
     *         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)            00001600
            KP = IABS(KPVT(K))                                          00001610
            IF (KP .EQ. K) GO TO 140                                    00001620
               T = Z(K)                                                 00001630
               Z(K) = Z(KP)                                             00001640
               Z(KP) = T                                                00001650
  140       CONTINUE                                                    00001660
  150    CONTINUE                                                       00001670
         K = K + KS                                                     00001680
      GO TO 130                                                         00001690
  160 CONTINUE                                                          00001700
      S = 1.0D0/DASUM(N,Z,1)                                            00001710
      CALL DSCAL(N,S,Z,1)                                               00001720
C                                                                       00001730
      YNORM = 1.0D0                                                     00001740
C                                                                       00001750
C     SOLVE U*D*V = Y                                                   00001760
C                                                                       00001770
      K = N                                                             00001780
  170 IF (K .EQ. 0) GO TO 230                                           00001790
         KS = 1                                                         00001800
         IF (KPVT(K) .LT. 0) KS = 2                                     00001810
         IF (K .EQ. KS) GO TO 190                                       00001820
            KP = IABS(KPVT(K))                                          00001830
            KPS = K + 1 - KS                                            00001840
            IF (KP .EQ. KPS) GO TO 180                                  00001850
               T = Z(KPS)                                               00001860
               Z(KPS) = Z(KP)                                           00001870
               Z(KP) = T                                                00001880
  180       CONTINUE                                                    00001890
            CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                       00001900
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)    00001910
  190    CONTINUE                                                       00001920
         IF (KS .EQ. 2) GO TO 210                                       00001930
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 200                 00001940
               S = DABS(A(K,K))/DABS(Z(K))                              00001950
               CALL DSCAL(N,S,Z,1)                                      00001960
               YNORM = S*YNORM                                          00001970
  200       CONTINUE                                                    00001980
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   00001990
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         00002000
         GO TO 220                                                      00002010
  210    CONTINUE                                                       00002020
            AK = A(K,K)/A(K-1,K)                                        00002030
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  00002040
            BK = Z(K)/A(K-1,K)                                          00002050
            BKM1 = Z(K-1)/A(K-1,K)                                      00002060
            DENOM = AK*AKM1 - 1.0D0                                     00002070
            Z(K) = (AKM1*BK - BKM1)/DENOM                               00002080
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               00002090
  220    CONTINUE                                                       00002100
         K = K - KS                                                     00002110
      GO TO 170                                                         00002120
  230 CONTINUE                                                          00002130
      S = 1.0D0/DASUM(N,Z,1)                                            00002140
      CALL DSCAL(N,S,Z,1)                                               00002150
      YNORM = S*YNORM                                                   00002160
C                                                                       00002170
C     SOLVE TRANS(U)*Z = V                                              00002180
C                                                                       00002190
      K = 1                                                             00002200
  240 IF (K .GT. N) GO TO 270                                           00002210
         KS = 1                                                         00002220
         IF (KPVT(K) .LT. 0) KS = 2                                     00002230
         IF (K .EQ. 1) GO TO 260                                        00002240
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)                     00002250
            IF (KS .EQ. 2)                                              00002260
     *         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)            00002270
            KP = IABS(KPVT(K))                                          00002280
            IF (KP .EQ. K) GO TO 250                                    00002290
               T = Z(K)                                                 00002300
               Z(K) = Z(KP)                                             00002310
               Z(KP) = T                                                00002320
  250       CONTINUE                                                    00002330
  260    CONTINUE                                                       00002340
         K = K + KS                                                     00002350
      GO TO 240                                                         00002360
  270 CONTINUE                                                          00002370
C     MAKE ZNORM = 1.0                                                  00002380
      S = 1.0D0/DASUM(N,Z,1)                                            00002390
      CALL DSCAL(N,S,Z,1)                                               00002400
      YNORM = S*YNORM                                                   00002410
C                                                                       00002420
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         00002430
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               00002440
      RETURN                                                            00002450
      END                                                               00002460

      SUBROUTINE DSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)                 00000010
      INTEGER LDA,N,JOB                                                 00000020
      DOUBLE PRECISION A(LDA,1),WORK(1)                                 00000030
      DOUBLE PRECISION DET(2)                                           00000040
      INTEGER KPVT(1),INERT(3)                                          00000050
C                                                                       00000060
C     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE               00000070
C     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM     00000080
C     DSIFA.                                                            00000090
C                                                                       00000100
C     ON ENTRY                                                          00000110
C                                                                       00000120
C        A       DOUBLE PRECISION(LDA,N)                                00000130
C                THE OUTPUT FROM DSIFA.                                 00000140
C                                                                       00000150
C        LDA     INTEGER                                                00000160
C                THE LEADING DIMENSION OF THE ARRAY A.                  00000170
C                                                                       00000180
C        N       INTEGER                                                00000190
C                THE ORDER OF THE MATRIX A.                             00000200
C                                                                       00000210
C        KPVT    INTEGER(N)                                             00000220
C                THE PIVOT VECTOR FROM DSIFA.                           00000230
C                                                                       00000240
C        WORK    DOUBLE PRECISION(N)                                    00000250
C                WORK VECTOR.  CONTENTS DESTROYED.                      00000260
C                                                                       00000270
C        JOB     INTEGER                                                00000280
C                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE              00000290
C                   IF  C .NE. 0, THE INVERSE IS COMPUTED,              00000300
C                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,          00000310
C                   IF  A .NE. 0, THE INERTIA IS COMPUTED.              00000320
C                                                                       00000330
C                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.               00000340
C                                                                       00000350
C     ON RETURN                                                         00000360
C                                                                       00000370
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.                   00000380
C                                                                       00000390
C        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF           00000400
C               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE         00000410
C               IS NEVER REFERENCED.                                    00000420
C                                                                       00000430
C        DET    DOUBLE PRECISION(2)                                     00000440
C               DETERMINANT OF ORIGINAL MATRIX.                         00000450
C               DETERMINANT = DET(1) * 10.0**DET(2)                     00000460
C               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0                    00000470
C               OR DET(1) = 0.0.                                        00000480
C                                                                       00000490
C        INERT  INTEGER(3)                                              00000500
C               THE INERTIA OF THE ORIGINAL MATRIX.                     00000510
C               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.            00000520
C               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.            00000530
C               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.                00000540
C                                                                       00000550
C     ERROR CONDITION                                                   00000560
C                                                                       00000570
C        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED       00000580
C        AND  DSICO  HAS SET RCOND .EQ. 0.0                             00000590
C        OR  DSIFA  HAS SET  INFO .NE. 0 .                              00000600
C                                                                       00000610
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000620
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB             00000630
C                                                                       00000640
C     SUBROUTINES AND FUNCTIONS                                         00000650
C                                                                       00000660
C     BLAS DAXPY,DCOPY,DDOT,DSWAP                                       00000670
C     FORTRAN DABS,IABS,MOD                                             00000680
C                                                                       00000690
C     INTERNAL VARIABLES.                                               00000700
C                                                                       00000710
      DOUBLE PRECISION AKKP1,DDOT,TEMP                                  00000720
      DOUBLE PRECISION TEN,D,T,AK,AKP1                                  00000730
      INTEGER J,JB,K,KM1,KS,KSTEP                                       00000740
      LOGICAL NOINV,NODET,NOERT                                         00000750
C                                                                       00000760
      NOINV = MOD(JOB,10) .EQ. 0                                        00000770
      NODET = MOD(JOB,100)/10 .EQ. 0                                    00000780
      NOERT = MOD(JOB,1000)/100 .EQ. 0                                  00000790
C                                                                       00000800
      IF (NODET .AND. NOERT) GO TO 140                                  00000810
         IF (NOERT) GO TO 10                                            00000820
            INERT(1) = 0                                                00000830
            INERT(2) = 0                                                00000840
            INERT(3) = 0                                                00000850
   10    CONTINUE                                                       00000860
         IF (NODET) GO TO 20                                            00000870
            DET(1) = 1.0D0                                              00000880
            DET(2) = 0.0D0                                              00000890
            TEN = 10.0D0                                                00000900
   20    CONTINUE                                                       00000910
         T = 0.0D0                                                      00000920
         DO 130 K = 1, N                                                00000930
            D = A(K,K)                                                  00000940
C                                                                       00000950
C           CHECK IF 1 BY 1                                             00000960
C                                                                       00000970
            IF (KPVT(K) .GT. 0) GO TO 50                                00000980
C                                                                       00000990
C              2 BY 2 BLOCK                                             00001000
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)     00001010
C                      (S  C)                                           00001020
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.                    00001030
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.      00001040
C                                                                       00001050
               IF (T .NE. 0.0D0) GO TO 30                               00001060
                  T = DABS(A(K,K+1))                                    00001070
                  D = (D/T)*A(K+1,K+1) - T                              00001080
               GO TO 40                                                 00001090
   30          CONTINUE                                                 00001100
                  D = T                                                 00001110
                  T = 0.0D0                                             00001120
   40          CONTINUE                                                 00001130
   50       CONTINUE                                                    00001140
C                                                                       00001150
            IF (NOERT) GO TO 60                                         00001160
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1                00001170
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1                00001180
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1                00001190
   60       CONTINUE                                                    00001200
C                                                                       00001210
            IF (NODET) GO TO 120                                        00001220
               DET(1) = D*DET(1)                                        00001230
               IF (DET(1) .EQ. 0.0D0) GO TO 110                         00001240
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80                 00001250
                     DET(1) = TEN*DET(1)                                00001260
                     DET(2) = DET(2) - 1.0D0                            00001270
                  GO TO 70                                              00001280
   80             CONTINUE                                              00001290
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100                  00001300
                     DET(1) = DET(1)/TEN                                00001310
                     DET(2) = DET(2) + 1.0D0                            00001320
                  GO TO 90                                              00001330
  100             CONTINUE                                              00001340
  110          CONTINUE                                                 00001350
  120       CONTINUE                                                    00001360
  130    CONTINUE                                                       00001370
  140 CONTINUE                                                          00001380
C                                                                       00001390
C     COMPUTE INVERSE(A)                                                00001400
C                                                                       00001410
      IF (NOINV) GO TO 270                                              00001420
         K = 1                                                          00001430
  150    IF (K .GT. N) GO TO 260                                        00001440
            KM1 = K - 1                                                 00001450
            IF (KPVT(K) .LT. 0) GO TO 180                               00001460
C                                                                       00001470
C              1 BY 1                                                   00001480
C                                                                       00001490
               A(K,K) = 1.0D0/A(K,K)                                    00001500
               IF (KM1 .LT. 1) GO TO 170                                00001510
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)                       00001520
                  DO 160 J = 1, KM1                                     00001530
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)                   00001540
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)          00001550
  160             CONTINUE                                              00001560
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)           00001570
  170          CONTINUE                                                 00001580
               KSTEP = 1                                                00001590
            GO TO 220                                                   00001600
  180       CONTINUE                                                    00001610
C                                                                       00001620
C              2 BY 2                                                   00001630
C                                                                       00001640
               T = DABS(A(K,K+1))                                       00001650
               AK = A(K,K)/T                                            00001660
               AKP1 = A(K+1,K+1)/T                                      00001670
               AKKP1 = A(K,K+1)/T                                       00001680
               D = T*(AK*AKP1 - 1.0D0)                                  00001690
               A(K,K) = AKP1/D                                          00001700
               A(K+1,K+1) = AK/D                                        00001710
               A(K,K+1) = -AKKP1/D                                      00001720
               IF (KM1 .LT. 1) GO TO 210                                00001730
                  CALL DCOPY(KM1,A(1,K+1),1,WORK,1)                     00001740
                  DO 190 J = 1, KM1                                     00001750
                     A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)                 00001760
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)        00001770
  190             CONTINUE                                              00001780
                  A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1) 00001790
                  A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)   00001800
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)                       00001810
                  DO 200 J = 1, KM1                                     00001820
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)                   00001830
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)          00001840
  200             CONTINUE                                              00001850
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)           00001860
  210          CONTINUE                                                 00001870
               KSTEP = 2                                                00001880
  220       CONTINUE                                                    00001890
C                                                                       00001900
C           SWAP                                                        00001910
C                                                                       00001920
            KS = IABS(KPVT(K))                                          00001930
            IF (KS .EQ. K) GO TO 250                                    00001940
               CALL DSWAP(KS,A(1,KS),1,A(1,K),1)                        00001950
               DO 230 JB = KS, K                                        00001960
                  J = K + KS - JB                                       00001970
                  TEMP = A(J,K)                                         00001980
                  A(J,K) = A(KS,J)                                      00001990
                  A(KS,J) = TEMP                                        00002000
  230          CONTINUE                                                 00002010
               IF (KSTEP .EQ. 1) GO TO 240                              00002020
                  TEMP = A(KS,K+1)                                      00002030
                  A(KS,K+1) = A(K,K+1)                                  00002040
                  A(K,K+1) = TEMP                                       00002050
  240          CONTINUE                                                 00002060
  250       CONTINUE                                                    00002070
            K = K + KSTEP                                               00002080
         GO TO 150                                                      00002090
  260    CONTINUE                                                       00002100
  270 CONTINUE                                                          00002110
      RETURN                                                            00002120
      END                                                               00002130

      SUBROUTINE DSISL(A,LDA,N,KPVT,B)                                  00000010
      INTEGER LDA,N,KPVT(1)                                             00000020
      DOUBLE PRECISION A(LDA,1),B(1)                                    00000030
C                                                                       00000040
C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM                00000050
C     A * X = B                                                         00000060
C     USING THE FACTORS COMPUTED BY DSIFA.                              00000070
C                                                                       00000080
C     ON ENTRY                                                          00000090
C                                                                       00000100
C        A       DOUBLE PRECISION(LDA,N)                                00000110
C                THE OUTPUT FROM DSIFA.                                 00000120
C                                                                       00000130
C        LDA     INTEGER                                                00000140
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000150
C                                                                       00000160
C        N       INTEGER                                                00000170
C                THE ORDER OF THE MATRIX  A .                           00000180
C                                                                       00000190
C        KPVT    INTEGER(N)                                             00000200
C                THE PIVOT VECTOR FROM DSIFA.                           00000210
C                                                                       00000220
C        B       DOUBLE PRECISION(N)                                    00000230
C                THE RIGHT HAND SIDE VECTOR.                            00000240
C                                                                       00000250
C     ON RETURN                                                         00000260
C                                                                       00000270
C        B       THE SOLUTION VECTOR  X .                               00000280
C                                                                       00000290
C     ERROR CONDITION                                                   00000300
C                                                                       00000310
C        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0 00000320
C        OR  DSIFA  HAS SET INFO .NE. 0  .                              00000330
C                                                                       00000340
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 00000350
C     WITH  P  COLUMNS                                                  00000360
C           CALL DSIFA(A,LDA,N,KPVT,INFO)                               00000370
C           IF (INFO .NE. 0) GO TO ...                                  00000380
C           DO 10 J = 1, P                                              00000390
C              CALL DSISL(A,LDA,N,KPVT,C(1,J))                          00000400
C        10 CONTINUE                                                    00000410
C                                                                       00000420
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000430
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            00000440
C                                                                       00000450
C     SUBROUTINES AND FUNCTIONS                                         00000460
C                                                                       00000470
C     BLAS DAXPY,DDOT                                                   00000480
C     FORTRAN IABS                                                      00000490
C                                                                       00000500
C     INTERNAL VARIABLES.                                               00000510
C                                                                       00000520
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP                  00000530
      INTEGER K,KP                                                      00000540
C                                                                       00000550
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND                    00000560
C     D INVERSE TO B.                                                   00000570
C                                                                       00000580
      K = N                                                             00000590
   10 IF (K .EQ. 0) GO TO 80                                            00000600
         IF (KPVT(K) .LT. 0) GO TO 40                                   00000610
C                                                                       00000620
C           1 X 1 PIVOT BLOCK.                                          00000630
C                                                                       00000640
            IF (K .EQ. 1) GO TO 30                                      00000650
               KP = KPVT(K)                                             00000660
               IF (KP .EQ. K) GO TO 20                                  00000670
C                                                                       00000680
C                 INTERCHANGE.                                          00000690
C                                                                       00000700
                  TEMP = B(K)                                           00000710
                  B(K) = B(KP)                                          00000720
                  B(KP) = TEMP                                          00000730
   20          CONTINUE                                                 00000740
C                                                                       00000750
C              APPLY THE TRANSFORMATION.                                00000760
C                                                                       00000770
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)                     00000780
   30       CONTINUE                                                    00000790
C                                                                       00000800
C           APPLY D INVERSE.                                            00000810
C                                                                       00000820
            B(K) = B(K)/A(K,K)                                          00000830
            K = K - 1                                                   00000840
         GO TO 70                                                       00000850
   40    CONTINUE                                                       00000860
C                                                                       00000870
C           2 X 2 PIVOT BLOCK.                                          00000880
C                                                                       00000890
            IF (K .EQ. 2) GO TO 60                                      00000900
               KP = IABS(KPVT(K))                                       00000910
               IF (KP .EQ. K - 1) GO TO 50                              00000920
C                                                                       00000930
C                 INTERCHANGE.                                          00000940
C                                                                       00000950
                  TEMP = B(K-1)                                         00000960
                  B(K-1) = B(KP)                                        00000970
                  B(KP) = TEMP                                          00000980
   50          CONTINUE                                                 00000990
C                                                                       00001000
C              APPLY THE TRANSFORMATION.                                00001010
C                                                                       00001020
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)                     00001030
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)                 00001040
   60       CONTINUE                                                    00001050
C                                                                       00001060
C           APPLY D INVERSE.                                            00001070
C                                                                       00001080
            AK = A(K,K)/A(K-1,K)                                        00001090
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  00001100
            BK = B(K)/A(K-1,K)                                          00001110
            BKM1 = B(K-1)/A(K-1,K)                                      00001120
            DENOM = AK*AKM1 - 1.0D0                                     00001130
            B(K) = (AKM1*BK - BKM1)/DENOM                               00001140
            B(K-1) = (AK*BKM1 - BK)/DENOM                               00001150
            K = K - 2                                                   00001160
   70    CONTINUE                                                       00001170
      GO TO 10                                                          00001180
   80 CONTINUE                                                          00001190
C                                                                       00001200
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.                        00001210
C                                                                       00001220
      K = 1                                                             00001230
   90 IF (K .GT. N) GO TO 160                                           00001240
         IF (KPVT(K) .LT. 0) GO TO 120                                  00001250
C                                                                       00001260
C           1 X 1 PIVOT BLOCK.                                          00001270
C                                                                       00001280
            IF (K .EQ. 1) GO TO 110                                     00001290
C                                                                       00001300
C              APPLY THE TRANSFORMATION.                                00001310
C                                                                       00001320
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)                  00001330
               KP = KPVT(K)                                             00001340
               IF (KP .EQ. K) GO TO 100                                 00001350
C                                                                       00001360
C                 INTERCHANGE.                                          00001370
C                                                                       00001380
                  TEMP = B(K)                                           00001390
                  B(K) = B(KP)                                          00001400
                  B(KP) = TEMP                                          00001410
  100          CONTINUE                                                 00001420
  110       CONTINUE                                                    00001430
            K = K + 1                                                   00001440
         GO TO 150                                                      00001450
  120    CONTINUE                                                       00001460
C                                                                       00001470
C           2 X 2 PIVOT BLOCK.                                          00001480
C                                                                       00001490
            IF (K .EQ. 1) GO TO 140                                     00001500
C                                                                       00001510
C              APPLY THE TRANSFORMATION.                                00001520
C                                                                       00001530
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)                  00001540
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)            00001550
               KP = IABS(KPVT(K))                                       00001560
               IF (KP .EQ. K) GO TO 130                                 00001570
C                                                                       00001580
C                 INTERCHANGE.                                          00001590
C                                                                       00001600
                  TEMP = B(K)                                           00001610
                  B(K) = B(KP)                                          00001620
                  B(KP) = TEMP                                          00001630
  130          CONTINUE                                                 00001640
  140       CONTINUE                                                    00001650
            K = K + 2                                                   00001660
  150    CONTINUE                                                       00001670
      GO TO 90                                                          00001680
  160 CONTINUE                                                          00001690
      RETURN                                                            00001700
      END                                                               00001710

C                                                                       90220001
C     ------------------------------------------------------------------90220002
C                                                                       90220003
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  90220004
C                                                                       90220005
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                             90220006
      REAL*8 D(N),E(N),Z(NM,N)                                          90220007
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     90220008
      REAL*8 DSQRT,DABS,DSIGN                                           90220009
C                                                                       90220010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     90220011
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     90220012
C     WILKINSON.                                                        90220013
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   90220014
C                                                                       90220015
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            90220016
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               90220017
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              90220018
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  90220019
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  90220020
C                                                                       90220021
C     ON INPUT:                                                         90220022
C                                                                       90220023
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         90220024
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          90220025
C          DIMENSION STATEMENT;                                         90220026
C                                                                       90220027
C        N IS THE ORDER OF THE MATRIX;                                  90220028
C                                                                       90220029
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          90220030
C                                                                       90220031
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        90220032
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;               90220033
C                                                                       90220034
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           90220035
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      90220036
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        90220037
C          THE IDENTITY MATRIX.                                         90220038
C                                                                       90220039
C      ON OUTPUT:                                                       90220040
C                                                                       90220041
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          90220042
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          90220043
C          UNORDERED FOR INDICES 1,2,...,IERR-1;                        90220044
C                                                                       90220045
C        E HAS BEEN DESTROYED;                                          90220046
C                                                                       90220047
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           90220048
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     90220049
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       90220050
C          EIGENVALUES;                                                 90220051
C                                                                       90220052
C        IERR IS SET TO                                                 90220053
C          ZERO       FOR NORMAL RETURN,                                90220054
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               90220055
C                     DETERMINED AFTER 30 ITERATIONS.                   90220056
C                                                                       90220057
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        90220058
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         90220059
C                                                                       90220060
C     ------------------------------------------------------------------90220061
C                                                                       90220062
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     90220063
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   90220064
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        90220065
C                ON S360 ::::::::::                                     90220066
C      DATA MACHEP/Z3410000000000000/                                    90220067
       MACHEP = 16.0D0**(-13)
C                                                                       90220068
      IERR = 0                                                          90220069
      IF (N .EQ. 1) GO TO 1001                                          90220070
C                                                                       90220071
      DO 100 I = 2, N                                                   90220072
  100 E(I-1) = E(I)                                                     90220073
C                                                                       90220074
      F = 0.0D0                                                         90220075
      B = 0.0D0                                                         90220076
      E(N) = 0.0D0                                                      90220077
C                                                                       90220078
      DO 240 L = 1, N                                                   90220079
         J = 0                                                          90220080
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))                         90220081
         IF (B .LT. H) B = H                                            90220082
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::         90220083
         DO 110 M = L, N                                                90220084
            IF (DABS(E(M)) .LE. B) GO TO 120                            90220085
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               90220086
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::              90220087
  110    CONTINUE                                                       90220088
C                                                                       90220089
  120    IF (M .EQ. L) GO TO 220                                        90220090
  130    IF (J .EQ. 30) GO TO 1000                                      90220091
         J = J + 1                                                      90220092
C     :::::::::: FORM SHIFT ::::::::::                                  90220093
         L1 = L + 1                                                     90220094
         G = D(L)                                                       90220095
         P = (D(L1) - G) / (2.0D0 * E(L))                               90220096
         R = DSQRT(P*P+1.0D0)                                           90220097
         D(L) = E(L) / (P + DSIGN(R,P))                                 90220098
         H = G - D(L)                                                   90220099
C                                                                       90220100
         DO 140 I = L1, N                                               90220101
  140    D(I) = D(I) - H                                                90220102
C                                                                       90220103
         F = F + H                                                      90220104
C     :::::::::: QL TRANSFORMATION ::::::::::                           90220105
         P = D(M)                                                       90220106
         C = 1.0D0                                                      90220107
         S = 0.0D0                                                      90220108
         MML = M - L                                                    90220109
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::             90220110
         DO 200 II = 1, MML                                             90220111
            I = M - II                                                  90220112
            G = C * E(I)                                                90220113
            H = C * P                                                   90220114
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                      90220115
            C = E(I) / P                                                90220116
            R = DSQRT(C*C+1.0D0)                                        90220117
            E(I+1) = S * P * R                                          90220118
            S = C / R                                                   90220119
            C = 1.0D0 / R                                               90220120
            GO TO 160                                                   90220121
  150       C = P / E(I)                                                90220122
            R = DSQRT(C*C+1.0D0)                                        90220123
            E(I+1) = S * E(I) * R                                       90220124
            S = 1.0D0 / R                                               90220125
            C = C * S                                                   90220126
  160       P = C * D(I) - S * G                                        90220127
            D(I+1) = H + S * (C * G + S * D(I))                         90220128
C     :::::::::: FORM VECTOR ::::::::::                                 90220129
            DO 180 K = 1, N                                             90220130
               H = Z(K,I+1)                                             90220131
               Z(K,I+1) = S * Z(K,I) + C * H                            90220132
               Z(K,I) = C * Z(K,I) - S * H                              90220133
  180       CONTINUE                                                    90220134
C                                                                       90220135
  200    CONTINUE                                                       90220136
C                                                                       90220137
         E(L) = S * P                                                   90220138
         D(L) = C * P                                                   90220139
         IF (DABS(E(L)) .GT. B) GO TO 130                               90220140
  220    D(L) = D(L) + F                                                90220141
  240 CONTINUE                                                          90220142
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::          90220143
      DO 300 II = 2, N                                                  90220144
         I = II - 1                                                     90220145
         K = I                                                          90220146
         P = D(I)                                                       90220147
C                                                                       90220148
         DO 260 J = II, N                                               90220149
            IF (D(J) .GE. P) GO TO 260                                  90220150
            K = J                                                       90220151
            P = D(J)                                                    90220152
  260    CONTINUE                                                       90220153
C                                                                       90220154
         IF (K .EQ. I) GO TO 300                                        90220155
         D(K) = D(I)                                                    90220156
         D(I) = P                                                       90220157
C                                                                       90220158
         DO 280 J = 1, N                                                90220159
            P = Z(J,I)                                                  90220160
            Z(J,I) = Z(J,K)                                             90220161
            Z(J,K) = P                                                  90220162
  280    CONTINUE                                                       90220163
C                                                                       90220164
  300 CONTINUE                                                          90220165
C                                                                       90220166
      GO TO 1001                                                        90220167
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN                      90220168
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::              90220169
 1000 IERR = L                                                          90220170
 1001 RETURN                                                            90220171
C     :::::::::: LAST CARD OF TQL2 ::::::::::                           90220172
      END                                                               90220173

C                                                                       78440001
C     ------------------------------------------------------------------78440002
C                                                                       78440003
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    78440004
C                                                                       78440005
      INTEGER I,J,K,L,N,II,NM,JP1                                       78440006
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)                                  78440007
      REAL*8 F,G,H,HH,SCALE                                             78440008
      REAL*8 DSQRT,DABS,DSIGN                                           78440009
C                                                                       78440010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    78440011
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   78440012
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   78440013
C                                                                       78440014
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              78440015
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               78440016
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            78440017
C                                                                       78440018
C     ON INPUT:                                                         78440019
C                                                                       78440020
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         78440021
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          78440022
C          DIMENSION STATEMENT;                                         78440023
C                                                                       78440024
C        N IS THE ORDER OF THE MATRIX;                                  78440025
C                                                                       78440026
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          78440027
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               78440028
C                                                                       78440029
C     ON OUTPUT:                                                        78440030
C                                                                       78440031
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;    78440032
C                                                                       78440033
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         78440034
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;      78440035
C                                                                       78440036
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                78440037
C          PRODUCED IN THE REDUCTION;                                   78440038
C                                                                       78440039
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            78440040
C                                                                       78440041
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        78440042
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         78440043
C                                                                       78440044
C     ------------------------------------------------------------------78440045
C                                                                       78440046
      DO 100 I = 1, N                                                   78440047
C                                                                       78440048
         DO 100 J = 1, I                                                78440049
            Z(I,J) = A(I,J)                                             78440050
  100 CONTINUE                                                          78440051
C                                                                       78440052
      IF (N .EQ. 1) GO TO 320                                           78440053
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::               78440054
      DO 300 II = 2, N                                                  78440055
         I = N + 2 - II                                                 78440056
         L = I - 1                                                      78440057
         H = 0.0D0                                                      78440058
         SCALE = 0.0D0                                                  78440059
         IF (L .LT. 2) GO TO 130                                        78440060
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::       78440061
         DO 120 K = 1, L                                                78440062
  120    SCALE = SCALE + DABS(Z(I,K))                                   78440063
C                                                                       78440064
         IF (SCALE .NE. 0.0D0) GO TO 140                                78440065
  130    E(I) = Z(I,L)                                                  78440066
         GO TO 290                                                      78440067
C                                                                       78440068
  140    DO 150 K = 1, L                                                78440069
            Z(I,K) = Z(I,K) / SCALE                                     78440070
            H = H + Z(I,K) * Z(I,K)                                     78440071
  150    CONTINUE                                                       78440072
C                                                                       78440073
         F = Z(I,L)                                                     78440074
         G = -DSIGN(DSQRT(H),F)                                         78440075
         E(I) = SCALE * G                                               78440076
         H = H - F * G                                                  78440077
         Z(I,L) = F - G                                                 78440078
         F = 0.0D0                                                      78440079
C                                                                       78440080
         DO 240 J = 1, L                                                78440081
            Z(J,I) = Z(I,J) / H                                         78440082
            G = 0.0D0                                                   78440083
C     :::::::::: FORM ELEMENT OF A*U ::::::::::                         78440084
            DO 180 K = 1, J                                             78440085
  180       G = G + Z(J,K) * Z(I,K)                                     78440086
C                                                                       78440087
            JP1 = J + 1                                                 78440088
            IF (L .LT. JP1) GO TO 220                                   78440089
C                                                                       78440090
            DO 200 K = JP1, L                                           78440091
  200       G = G + Z(K,J) * Z(I,K)                                     78440092
C     :::::::::: FORM ELEMENT OF P ::::::::::                           78440093
  220       E(J) = G / H                                                78440094
            F = F + E(J) * Z(I,J)                                       78440095
  240    CONTINUE                                                       78440096
C                                                                       78440097
         HH = F / (H + H)                                               78440098
C     :::::::::: FORM REDUCED A ::::::::::                              78440099
         DO 260 J = 1, L                                                78440100
            F = Z(I,J)                                                  78440101
            G = E(J) - HH * F                                           78440102
            E(J) = G                                                    78440103
C                                                                       78440104
            DO 260 K = 1, J                                             78440105
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                  78440106
  260    CONTINUE                                                       78440107
C                                                                       78440108
  290    D(I) = H                                                       78440109
  300 CONTINUE                                                          78440110
C                                                                       78440111
  320 D(1) = 0.0D0                                                      78440112
      E(1) = 0.0D0                                                      78440113
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::     78440114
      DO 500 I = 1, N                                                   78440115
         L = I - 1                                                      78440116
         IF (D(I) .EQ. 0.0D0) GO TO 380                                 78440117
C                                                                       78440118
         DO 360 J = 1, L                                                78440119
            G = 0.0D0                                                   78440120
C                                                                       78440121
            DO 340 K = 1, L                                             78440122
  340       G = G + Z(I,K) * Z(K,J)                                     78440123
C                                                                       78440124
            DO 360 K = 1, L                                             78440125
               Z(K,J) = Z(K,J) - G * Z(K,I)                             78440126
  360    CONTINUE                                                       78440127
C                                                                       78440128
  380    D(I) = Z(I,I)                                                  78440129
         Z(I,I) = 1.0D0                                                 78440130
         IF (L .LT. 1) GO TO 500                                        78440131
C                                                                       78440132
         DO 400 J = 1, L                                                78440133
            Z(I,J) = 0.0D0                                              78440134
            Z(J,I) = 0.0D0                                              78440135
  400    CONTINUE                                                       78440136
C                                                                       78440137
  500 CONTINUE                                                          78440138
C                                                                       78440139
      RETURN                                                            78440140
C     :::::::::: LAST CARD OF TRED2 ::::::::::                          78440141
      END                                                               78440142

      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)                        00000010
C                                                                       00000020
C     TAKES THE SUM OF THE ABSOLUTE VALUES.                             00000030
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000040
C                                                                       00000050
      DOUBLE PRECISION DX(1),DTEMP                                      00000060
      INTEGER I,INCX,M,MP1,N,NINCX                                      00000070
C                                                                       00000080
      DASUM = 0.0D0                                                     00000090
      DTEMP = 0.0D0                                                     00000100
      IF(N.LE.0)RETURN                                                  00000110
      IF(INCX.EQ.1)GO TO 20                                             00000120
C                                                                       00000130
C        CODE FOR INCREMENT NOT EQUAL TO 1                              00000140
C                                                                       00000150
      NINCX = N*INCX                                                    00000160
      DO 10 I = 1,NINCX,INCX                                            00000170
        DTEMP = DTEMP + DABS(DX(I))                                     00000180
   10 CONTINUE                                                          00000190
      DASUM = DTEMP                                                     00000200
      RETURN                                                            00000210
C                                                                       00000220
C        CODE FOR INCREMENT EQUAL TO 1                                  00000230
C                                                                       00000240
C                                                                       00000250
C        CLEAN-UP LOOP                                                  00000260
C                                                                       00000270
   20 M = MOD(N,6)                                                      00000280
      IF( M .EQ. 0 ) GO TO 40                                           00000290
      DO 30 I = 1,M                                                     00000300
        DTEMP = DTEMP + DABS(DX(I))                                     00000310
   30 CONTINUE                                                          00000320
      IF( N .LT. 6 ) GO TO 60                                           00000330
   40 MP1 = M + 1                                                       00000340
      DO 50 I = MP1,N,6                                                 00000350
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2)) 00000360
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))           00000370
   50 CONTINUE                                                          00000380
   60 DASUM = DTEMP                                                     00000390
      RETURN                                                            00000400
      END                                                               00000410

      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)                            00000010
C                                                                       00000020
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DA                                   00000070
      INTEGER I,INCX,INCY,IXIY,M,MP1,N                                  00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF (DA .EQ. 0.0D0) RETURN                                         00000110
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000120
C                                                                       00000130
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                00000140
C          NOT EQUAL TO 1                                               00000150
C                                                                       00000160
      IX = 1                                                            00000170
      IY = 1                                                            00000180
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000190
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000200
      DO 10 I = 1,N                                                     00000210
        DY(IY) = DY(IY) + DA*DX(IX)                                     00000220
        IX = IX + INCX                                                  00000230
        IY = IY + INCY                                                  00000240
   10 CONTINUE                                                          00000250
      RETURN                                                            00000260
C                                                                       00000270
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            00000280
C                                                                       00000290
C                                                                       00000300
C        CLEAN-UP LOOP                                                  00000310
C                                                                       00000320
   20 M = MOD(N,4)                                                      00000330
      IF( M .EQ. 0 ) GO TO 40                                           00000340
      DO 30 I = 1,M                                                     00000350
        DY(I) = DY(I) + DA*DX(I)                                        00000360
   30 CONTINUE                                                          00000370
      IF( N .LT. 4 ) RETURN                                             00000380
   40 MP1 = M + 1                                                       00000390
      DO 50 I = MP1,N,4                                                 00000400
        DY(I) = DY(I) + DA*DX(I)                                        00000410
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)                            00000420
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)                            00000430
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)                            00000440
   50 CONTINUE                                                          00000450
      RETURN                                                            00000460
      END                                                               00000470

      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)                              00000010
C                                                                       00000020
C     COPIES A VECTOR, X, TO A VECTOR, Y.                               00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1)                                      00000070
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000110
C                                                                       00000120
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                00000130
C          NOT EQUAL TO 1                                               00000140
C                                                                       00000150
      IX = 1                                                            00000160
      IY = 1                                                            00000170
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000180
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000190
      DO 10 I = 1,N                                                     00000200
        DY(IY) = DX(IX)                                                 00000210
        IX = IX + INCX                                                  00000220
        IY = IY + INCY                                                  00000230
   10 CONTINUE                                                          00000240
      RETURN                                                            00000250
C                                                                       00000260
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            00000270
C                                                                       00000280
C                                                                       00000290
C        CLEAN-UP LOOP                                                  00000300
C                                                                       00000310
   20 M = MOD(N,7)                                                      00000320
      IF( M .EQ. 0 ) GO TO 40                                           00000330
      DO 30 I = 1,M                                                     00000340
        DY(I) = DX(I)                                                   00000350
   30 CONTINUE                                                          00000360
      IF( N .LT. 7 ) RETURN                                             00000370
   40 MP1 = M + 1                                                       00000380
      DO 50 I = MP1,N,7                                                 00000390
        DY(I) = DX(I)                                                   00000400
        DY(I + 1) = DX(I + 1)                                           00000410
        DY(I + 2) = DX(I + 2)                                           00000420
        DY(I + 3) = DX(I + 3)                                           00000430
        DY(I + 4) = DX(I + 4)                                           00000440
        DY(I + 5) = DX(I + 5)                                           00000450
        DY(I + 6) = DX(I + 6)                                           00000460
   50 CONTINUE                                                          00000470
      RETURN                                                            00000480
      END                                                               00000490

      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                 00000010
C                                                                       00000020
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                             00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                00000070
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 00000080
C                                                                       00000090
      DDOT = 0.0D0                                                      00000100
      DTEMP = 0.0D0                                                     00000110
      IF(N.LE.0)RETURN                                                  00000120
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000130
C                                                                       00000140
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                00000150
C          NOT EQUAL TO 1                                               00000160
C                                                                       00000170
      IX = 1                                                            00000180
      IY = 1                                                            00000190
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000200
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000210
      DO 10 I = 1,N                                                     00000220
        DTEMP = DTEMP + DX(IX)*DY(IY)                                   00000230
        IX = IX + INCX                                                  00000240
        IY = IY + INCY                                                  00000250
   10 CONTINUE                                                          00000260
      DDOT = DTEMP                                                      00000270
      RETURN                                                            00000280
C                                                                       00000290
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            00000300
C                                                                       00000310
C                                                                       00000320
C        CLEAN-UP LOOP                                                  00000330
C                                                                       00000340
   20 M = MOD(N,5)                                                      00000350
      IF( M .EQ. 0 ) GO TO 40                                           00000360
      DO 30 I = 1,M                                                     00000370
        DTEMP = DTEMP + DX(I)*DY(I)                                     00000380
   30 CONTINUE                                                          00000390
      IF( N .LT. 5 ) GO TO 60                                           00000400
   40 MP1 = M + 1                                                       00000410
      DO 50 I = MP1,N,5                                                 00000420
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +             00000430
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)00000440
   50 CONTINUE                                                          00000450
   60 DDOT = DTEMP                                                      00000460
      RETURN                                                            00000470
      END                                                               00000480

      SUBROUTINE  DSCAL(N,DA,DX,INCX)                                   00000010
C                                                                       00000020
C     SCALES A VECTOR BY A CONSTANT.                                    00000030
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DA,DX(1)                                         00000070
      INTEGER I,INCX,M,MP1,N,NINCX                                      00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF(INCX.EQ.1)GO TO 20                                             00000110
C                                                                       00000120
C        CODE FOR INCREMENT NOT EQUAL TO 1                              00000130
C                                                                       00000140
      NINCX = N*INCX                                                    00000150
      DO 10 I = 1,NINCX,INCX                                            00000160
        DX(I) = DA*DX(I)                                                00000170
   10 CONTINUE                                                          00000180
      RETURN                                                            00000190
C                                                                       00000200
C        CODE FOR INCREMENT EQUAL TO 1                                  00000210
C                                                                       00000220
C                                                                       00000230
C        CLEAN-UP LOOP                                                  00000240
C                                                                       00000250
   20 M = MOD(N,5)                                                      00000260
      IF( M .EQ. 0 ) GO TO 40                                           00000270
      DO 30 I = 1,M                                                     00000280
        DX(I) = DA*DX(I)                                                00000290
   30 CONTINUE                                                          00000300
      IF( N .LT. 5 ) RETURN                                             00000310
   40 MP1 = M + 1                                                       00000320
      DO 50 I = MP1,N,5                                                 00000330
        DX(I) = DA*DX(I)                                                00000340
        DX(I + 1) = DA*DX(I + 1)                                        00000350
        DX(I + 2) = DA*DX(I + 2)                                        00000360
        DX(I + 3) = DA*DX(I + 3)                                        00000370
        DX(I + 4) = DA*DX(I + 4)                                        00000380
   50 CONTINUE                                                          00000390
      RETURN                                                            00000400
      END                                                               00000410

      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)                               00000010
      INTEGER LDA,N,KPVT(1),INFO                                        00000020
      DOUBLE PRECISION A(LDA,1)                                         00000030
C                                                                       00000040
C     DSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION  00000050
C     WITH SYMMETRIC PIVOTING.                                          00000060
C                                                                       00000070
C     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.                        00000080
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.                 00000090
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.               00000100
C     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.                   00000110
C     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.                   00000120
C                                                                       00000130
C     ON ENTRY                                                          00000140
C                                                                       00000150
C        A       DOUBLE PRECISION(LDA,N)                                00000160
C                THE SYMMETRIC MATRIX TO BE FACTORED.                   00000170
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         00000180
C                                                                       00000190
C        LDA     INTEGER                                                00000200
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000210
C                                                                       00000220
C        N       INTEGER                                                00000230
C                THE ORDER OF THE MATRIX  A .                           00000240
C                                                                       00000250
C     ON RETURN                                                         00000260
C                                                                       00000270
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      00000280
C                WERE USED TO OBTAIN IT.                                00000290
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     00000300
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         00000310
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            00000320
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            00000330
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         00000340
C                                                                       00000350
C        KPVT    INTEGER(N)                                             00000360
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000370
C                                                                       00000380
C        INFO    INTEGER                                                00000390
C                = 0  NORMAL VALUE.                                     00000400
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS      00000410
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,       00000420
C                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY      00000430
C                     DIVIDE BY ZERO IF CALLED.                         00000440
C                                                                       00000450
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000460
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            00000470
C                                                                       00000480
C     SUBROUTINES AND FUNCTIONS                                         00000490
C                                                                       00000500
C     BLAS DAXPY,DSWAP,IDAMAX                                           00000510
C     FORTRAN DABS,DMAX1,DSQRT                                          00000520
C                                                                       00000530
C     INTERNAL VARIABLES                                                00000540
C                                                                       00000550
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T              00000560
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX                       00000570
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,IDAMAX              00000580
      LOGICAL SWAP                                                      00000590
C                                                                       00000600
C                                                                       00000610
C     INITIALIZE                                                        00000620
C                                                                       00000630
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.                       00000640
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0                             00000650
C                                                                       00000660
      INFO = 0                                                          00000670
C                                                                       00000680
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.                           00000690
C                                                                       00000700
      K = N                                                             00000710
   10 CONTINUE                                                          00000720
C                                                                       00000730
C        LEAVE THE LOOP IF K=0 OR K=1.                                  00000740
C                                                                       00000750
C     ...EXIT                                                           00000760
         IF (K .EQ. 0) GO TO 200                                        00000770
         IF (K .GT. 1) GO TO 20                                         00000780
            KPVT(1) = 1                                                 00000790
            IF (A(1,1) .EQ. 0.0D0) INFO = 1                             00000800
C     ......EXIT                                                        00000810
            GO TO 200                                                   00000820
   20    CONTINUE                                                       00000830
C                                                                       00000840
C        THIS SECTION OF CODE DETERMINES THE KIND OF                    00000850
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,            00000860
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND          00000870
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS                00000880
C        REQUIRED.                                                      00000890
C                                                                       00000900
         KM1 = K - 1                                                    00000910
         ABSAKK = DABS(A(K,K))                                          00000920
C                                                                       00000930
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN                  00000940
C        COLUMN K.                                                      00000950
C                                                                       00000960
         IMAX = IDAMAX(K-1,A(1,K),1)                                    00000970
         COLMAX = DABS(A(IMAX,K))                                       00000980
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30                         00000990
            KSTEP = 1                                                   00001000
            SWAP = .FALSE.                                              00001010
         GO TO 90                                                       00001020
   30    CONTINUE                                                       00001030
C                                                                       00001040
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN               00001050
C           ROW IMAX.                                                   00001060
C                                                                       00001070
            ROWMAX = 0.0D0                                              00001080
            IMAXP1 = IMAX + 1                                           00001090
            DO 40 J = IMAXP1, K                                         00001100
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))                   00001110
   40       CONTINUE                                                    00001120
            IF (IMAX .EQ. 1) GO TO 50                                   00001130
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)                        00001140
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))                00001150
   50       CONTINUE                                                    00001160
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60          00001170
               KSTEP = 1                                                00001180
               SWAP = .TRUE.                                            00001190
            GO TO 80                                                    00001200
   60       CONTINUE                                                    00001210
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70      00001220
               KSTEP = 1                                                00001230
               SWAP = .FALSE.                                           00001240
            GO TO 80                                                    00001250
   70       CONTINUE                                                    00001260
               KSTEP = 2                                                00001270
               SWAP = IMAX .NE. KM1                                     00001280
   80       CONTINUE                                                    00001290
   90    CONTINUE                                                       00001300
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100                 00001310
C                                                                       00001320
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.           00001330
C                                                                       00001340
            KPVT(K) = K                                                 00001350
            INFO = K                                                    00001360
         GO TO 190                                                      00001370
  100    CONTINUE                                                       00001380
         IF (KSTEP .EQ. 2) GO TO 140                                    00001390
C                                                                       00001400
C           1 X 1 PIVOT BLOCK.                                          00001410
C                                                                       00001420
            IF (.NOT.SWAP) GO TO 120                                    00001430
C                                                                       00001440
C              PERFORM AN INTERCHANGE.                                  00001450
C                                                                       00001460
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)                    00001470
               DO 110 JJ = IMAX, K                                      00001480
                  J = K + IMAX - JJ                                     00001490
                  T = A(J,K)                                            00001500
                  A(J,K) = A(IMAX,J)                                    00001510
                  A(IMAX,J) = T                                         00001520
  110          CONTINUE                                                 00001530
  120       CONTINUE                                                    00001540
C                                                                       00001550
C           PERFORM THE ELIMINATION.                                    00001560
C                                                                       00001570
            DO 130 JJ = 1, KM1                                          00001580
               J = K - JJ                                               00001590
               MULK = -A(J,K)/A(K,K)                                    00001600
               T = MULK                                                 00001610
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)                        00001620
               A(J,K) = MULK                                            00001630
  130       CONTINUE                                                    00001640
C                                                                       00001650
C           SET THE PIVOT ARRAY.                                        00001660
C                                                                       00001670
            KPVT(K) = K                                                 00001680
            IF (SWAP) KPVT(K) = IMAX                                    00001690
         GO TO 190                                                      00001700
  140    CONTINUE                                                       00001710
C                                                                       00001720
C           2 X 2 PIVOT BLOCK.                                          00001730
C                                                                       00001740
            IF (.NOT.SWAP) GO TO 160                                    00001750
C                                                                       00001760
C              PERFORM AN INTERCHANGE.                                  00001770
C                                                                       00001780
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)                  00001790
               DO 150 JJ = IMAX, KM1                                    00001800
                  J = KM1 + IMAX - JJ                                   00001810
                  T = A(J,K-1)                                          00001820
                  A(J,K-1) = A(IMAX,J)                                  00001830
                  A(IMAX,J) = T                                         00001840
  150          CONTINUE                                                 00001850
               T = A(K-1,K)                                             00001860
               A(K-1,K) = A(IMAX,K)                                     00001870
               A(IMAX,K) = T                                            00001880
  160       CONTINUE                                                    00001890
C                                                                       00001900
C           PERFORM THE ELIMINATION.                                    00001910
C                                                                       00001920
            KM2 = K - 2                                                 00001930
            IF (KM2 .EQ. 0) GO TO 180                                   00001940
               AK = A(K,K)/A(K-1,K)                                     00001950
               AKM1 = A(K-1,K-1)/A(K-1,K)                               00001960
               DENOM = 1.0D0 - AK*AKM1                                  00001970
               DO 170 JJ = 1, KM2                                       00001980
                  J = KM1 - JJ                                          00001990
                  BK = A(J,K)/A(K-1,K)                                  00002000
                  BKM1 = A(J,K-1)/A(K-1,K)                              00002010
                  MULK = (AKM1*BK - BKM1)/DENOM                         00002020
                  MULKM1 = (AK*BKM1 - BK)/DENOM                         00002030
                  T = MULK                                              00002040
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)                     00002050
                  T = MULKM1                                            00002060
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)                   00002070
                  A(J,K) = MULK                                         00002080
                  A(J,K-1) = MULKM1                                     00002090
  170          CONTINUE                                                 00002100
  180       CONTINUE                                                    00002110
C                                                                       00002120
C           SET THE PIVOT ARRAY.                                        00002130
C                                                                       00002140
            KPVT(K) = 1 - K                                             00002150
            IF (SWAP) KPVT(K) = -IMAX                                   00002160
            KPVT(K-1) = KPVT(K)                                         00002170
  190    CONTINUE                                                       00002180
         K = K - KSTEP                                                  00002190
      GO TO 10                                                          00002200
  200 CONTINUE                                                          00002210
      RETURN                                                            00002220
      END                                                               00002230

      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)                             00000010
C                                                                       00000020
C     INTERCHANGES TWO VECTORS.                                         00000030
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.                     00000040
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000050
C                                                                       00000060
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                00000070
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 00000080
C                                                                       00000090
      IF(N.LE.0)RETURN                                                  00000100
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               00000110
C                                                                       00000120
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       00000130
C         TO 1                                                          00000140
C                                                                       00000150
      IX = 1                                                            00000160
      IY = 1                                                            00000170
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 00000180
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 00000190
      DO 10 I = 1,N                                                     00000200
        DTEMP = DX(IX)                                                  00000210
        DX(IX) = DY(IY)                                                 00000220
        DY(IY) = DTEMP                                                  00000230
        IX = IX + INCX                                                  00000240
        IY = IY + INCY                                                  00000250
   10 CONTINUE                                                          00000260
      RETURN                                                            00000270
C                                                                       00000280
C       CODE FOR BOTH INCREMENTS EQUAL TO 1                             00000290
C                                                                       00000300
C                                                                       00000310
C       CLEAN-UP LOOP                                                   00000320
C                                                                       00000330
   20 M = MOD(N,3)                                                      00000340
      IF( M .EQ. 0 ) GO TO 40                                           00000350
      DO 30 I = 1,M                                                     00000360
        DTEMP = DX(I)                                                   00000370
        DX(I) = DY(I)                                                   00000380
        DY(I) = DTEMP                                                   00000390
   30 CONTINUE                                                          00000400
      IF( N .LT. 3 ) RETURN                                             00000410
   40 MP1 = M + 1                                                       00000420
      DO 50 I = MP1,N,3                                                 00000430
        DTEMP = DX(I)                                                   00000440
        DX(I) = DY(I)                                                   00000450
        DY(I) = DTEMP                                                   00000460
        DTEMP = DX(I + 1)                                               00000470
        DX(I + 1) = DY(I + 1)                                           00000480
        DY(I + 1) = DTEMP                                               00000490
        DTEMP = DX(I + 2)                                               00000500
        DX(I + 2) = DY(I + 2)                                           00000510
        DY(I + 2) = DTEMP                                               00000520
   50 CONTINUE                                                          00000530
      RETURN                                                            00000540
      END                                                               00000550

      INTEGER FUNCTION IDAMAX(N,DX,INCX)                                00000010
C                                                                       00000020
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            00000030
C     JACK DONGARRA, LINPACK, 3/11/78.                                  00000040
C                                                                       00000050
      DOUBLE PRECISION DX(1),DMAX                                       00000060
      INTEGER I,INCX,IX,N                                               00000070
C                                                                       00000080
      IDAMAX = 0                                                        00000090
      IF( N .LT. 1 ) RETURN                                             00000100
      IDAMAX = 1                                                        00000110
      IF(N.EQ.1)RETURN                                                  00000120
      IF(INCX.EQ.1)GO TO 20                                             00000130
C                                                                       00000140
C        CODE FOR INCREMENT NOT EQUAL TO 1                              00000150
C                                                                       00000160
      IX = 1                                                            00000170
      DMAX = DABS(DX(1))                                                00000180
      IX = IX + INCX                                                    00000190
      DO 10 I = 2,N                                                     00000200
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5                               00000210
         IDAMAX = I                                                     00000220
         DMAX = DABS(DX(IX))                                            00000230
    5    IX = IX + INCX                                                 00000240
   10 CONTINUE                                                          00000250
      RETURN                                                            00000260
C                                                                       00000270
C        CODE FOR INCREMENT EQUAL TO 1                                  00000280
C                                                                       00000290
   20 DMAX = DABS(DX(1))                                                00000300
      DO 30 I = 2,N                                                     00000310
         IF(DABS(DX(I)).LE.DMAX) GO TO 30                               00000320
         IDAMAX = I                                                     00000330
         DMAX = DABS(DX(I))                                             00000340
   30 CONTINUE                                                          00000350
      RETURN                                                            00000360
      END                                                               00000370

