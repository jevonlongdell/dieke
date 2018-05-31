        FUNCTION THREEJ(J1,J2,J3,M1,M2,M3)                            
C...08DEC88...MODIFIED TO CHECK FOR FACTORIAL VALUES GREATER THAN
C...LIMIT LFAC AND EXTEND LIMIT TO LFAC=66.
C...PERMIT ENTRY BY X = THRJ(J1,J2,J3,M1,M2,3) OR THREEJ() OR SIXJ()
C...G.L.GOODMAN, ARGONNE NATIONAL LABORATORY, B04451@ANLVM.BITNET
      IMPLICIT  REAL*8(A-H,O-Z)                                       
      ENTRY THRJ(J1,J2,J3,M1,M2,M3)
      S3J=0.                                                          
      CALL DELTA3(J1,J2,J3,D)                                         
      IF (M1+M2+M3.NE.0.OR.D.EQ.0) GOTO 9                             
      P=1-MOD(ABS(J1-J2-M3),4)                                        
      F1=F(J1+M1)                                                     
      F2=F(J1-M1)                                                     
      F3=F(J2+M2)                                                     
      F4=F(J2-M2)                                                     
      F5=F(J3+M3)                                                     
      F6=F(J3-M3)                                                     
      S=F1*F2*F3*F4*F5*F6                                             
      U=0.                                                              00009200
      L=MIN(J1+J2-J3,J1-M1,J2+M2)                                       00009300
      M=MAX(0,J2-J3-M1,J1-J3+M2)                                        00009400
      IF (L.LT.M) GO TO 9
      DO 1 K=M,L,2                                                      00009500
      A=1-MOD(K,4)                                                      00009600
      F1=F(K)                                                           00009700
      F2=F(J1+J2-J3-K)                                                  00009800
      F3=F(J3-J2+M1+K)                                                  00009900
      F4=F(J1-M1-K)                                                     00010000
      F5=F(J3-J1-M2+K)                                                  00010100
      F6=F(J2+M2-K)                                                     00010200
      FM=F1*F2*F3*F4*F5*F6                                              00010300
      U=U+A/FM                                                          00010400
    1 CONTINUE                                                          00010500
      S3J=P*D*SQRT(S)*U                                                 00010600
    9 THRJ = S3J
      THREEJ = S3J
      END                                                               00010700
      FUNCTION SIXJ(J1,J2,J3,L1,L2,L3)           
      IMPLICIT  REAL*8(A-H,O-Z)                                         00010900
      S6J=0.                                                            00011100
      CALL DELTA3(J1,J2,J3,D1)                                          00011200
      CALL DELTA3(J1,L2,L3,D2)                                          00011300
      CALL DELTA3(L1,J2,L3,D3)                                          00011400
      CALL DELTA3(L1,L2,J3,D4)                                          00011500
      IF(D1.LE.0..OR.D2.LE.0..OR.D3.LE.0..OR.D4.LE.0.)GO TO 9
      D=D1*D2*D3*D4                                                     00011600
      P=1-MOD(J1+J2+L1+L2,4)                                            00011800
      U=0.                                                              00011900
      L=MIN(J1+J2-J3,L1+L2-J3,J1+L2-L3,L1+J2-L3)                        00012000
      M=MAX(0,J1+L1-J3-L3,J2+L2-J3-L3)                                  00012100
      IF(L.LT.M) GO TO 9
      DO 1 K=M,L,2                                                      00012200
      A=1-MOD(K,4)                                                      00012300
      F0=F(J1+J2+L1+L2+2-K)                                             00012400
      F0=F0/F(K)                                                        00012500
      F0=F0/F(J1+J2-J3-K)                                               00012600
      F0=F0/F(L1+L2-J3-K)                                               00012700
      F0=F0/F(J1+L2-L3-K)                                               00012800
      F0=F0/F(L1+J2-L3-K)                                               00012900
      F0=F0/F(-J1-L1+J3+L3+K)                                           00013000
      F0=F0/F(-J2-L2+J3+L3+K)                                           00013100
      U=U+A*F0                                                          00013300
    1 CONTINUE                                                          00013400
      S6J=P*D*U                                                         00013500
    9 SIXJ = S6J
      END                                                               00013600
      SUBROUTINE DELTA3(J1,J2,J3,D)                                     00013700
      IMPLICIT  REAL*8(A-H,O-Z)                                         0001380
      D=0.                                                              00014000
      J0=(J1+J2+J3)/2                                                   00014100
      IF (J0.LT.J1.OR.J0.LT.J2.OR.J0.LT.J3) RETURN                      00014200
      F1=F(J1+J2-J3)                                                    00014300
      F2=F(J1-J2+J3)                                                    00014400
      F3=F(-J1+J2+J3)                                                   00014500
      F0=F(J1+J2+J3+2)                                                  00014600
      D=SQRT((F1/F0)*F2*F3)                                             00014700
      END                                                               00014800
      SUBROUTINE FACSYS                                                 00014900
      IMPLICIT  REAL*16(A-H,O-Z)                                        00015000
      COMMON/FAC/FC(66),LFAC
      LFAC = 66
      FA=1.                                                             00015200
      DO 1 I=1,LFAC                                                     00015300
      FA=FA*I                                                           00015400
      FC(I)=FA                                                          00015500
    1 CONTINUE                                                          00015600
      FC(0)=1.                                                          00015700
      END                                                               00015800
      FUNCTION F(I)
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*16 FC
      COMMON/FAC/FC(66),LFAC
      DATA IFIRST/0/
      IF (IFIRST.NE.0)GO  TO 10
      IFIRST=1
      CALL FACSYS
   10 I2 = I/2
      IF(I2.GT.LFAC.OR.I.LT.0.OR.I.NE.I2*2) GO TO 100
      F = FC(I2)
      RETURN
  100 F = FC(LFAC)
      WRITE(6,101) F,I2,LFAC
  101 FORMAT(' NB: ',E8.2,' USED FOR FACTORIAL ',I5,' LIMIT IS',I5)
      CONTINUE
      RETURN
      END 
C                                                                       23440001
C     ------------------------------------------------------------------23440002
C                                                                       23440003
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,                          23440004
     X                  IERR,RV1,RV2,RV3,RV4,RV6)                       23440005
C                                                                       23440006
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP            23440007
      REAL*8 D(N),E(N),E2(N),W(M),Z(NM,M),                              23440008
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)                         23440009
      REAL*8 U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER
      REAL*16 MACHEP           
      REAL*8 DSQRT,DABS,DFLOAT                                          23440011
      INTEGER IND(M)                                                    23440012
      REAL*16 JUNK
C                                                                       23440013
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-   23440014
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.    23440015
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).   23440016
C                                                                       23440017
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL         23440018
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,          23440019
C     USING INVERSE ITERATION.                                          23440020
C                                                                       23440021
C     ON INPUT:                                                         23440022
C                                                                       23440023
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         23440024
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          23440025
C          DIMENSION STATEMENT;                                         23440026
C                                                                       23440027
C        N IS THE ORDER OF THE MATRIX;                                  23440028
C                                                                       23440029
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          23440030
C                                                                       23440031
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        23440032
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;               23440033
C                                                                       23440034
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,    23440035
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.        23440036
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN       23440037
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM    23440038
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN    23440039
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0    23440040
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,     23440041
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,   23440042
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE;      23440043
C                                                                       23440044
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES;                      23440045
C                                                                       23440046
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER; 23440047
C                                                                       23440048
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES    23440049
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --        23440050
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM      23440051
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC. 23440052
C                                                                       23440053
C     ON OUTPUT:                                                        23440054
C                                                                       23440055
C        ALL INPUT ARRAYS ARE UNALTERED;                                23440056
C                                                                       23440057
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.     23440058
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO;           23440059
C                                                                       23440060
C        IERR IS SET TO                                                 23440061
C          ZERO       FOR NORMAL RETURN,                                23440062
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH      23440063
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS;     23440064
C                                                                       23440065
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.      23440066
C                                                                       23440067
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        23440068
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         23440069
C                                                                       23440070
C     ------------------------------------------------------------------23440071
C                                                                       23440072
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     23440073
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   23440074
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        23440075
C                ON S360 ::::::::::                                     23440076
C      DATA MACHEP/Z3410000000000000/                                    23440077
C       DATA MACHEP/1.D-17/
C
      JUNK =16.0D0**(-13)
      MACHEP = JUNK
C      PRINT *,'MACHEP =',MACHEP
C                                                                       23440078
      IERR = 0                                                          23440079
      IF (M .EQ. 0) GO TO 1001                                          23440080
      TAG = 0                                                           23440081
      ORDER = 1.0D0 - E2(1)                                             23440082
      Q = 0                                                             23440083
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX ::::::::::        23440084
  100 P = Q + 1                                                         23440085
C                                                                       23440086
      DO 120 Q = P, N                                                   23440087
         IF (Q .EQ. N) GO TO 140                                        23440088
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140                              23440089
  120 CONTINUE                                                          23440090
C     :::::::::: FIND VECTORS BY INVERSE ITERATION ::::::::::           23440091
  140 TAG = TAG + 1                                                     23440092
      S = 0                                                             23440093
C                                                                       23440094
      DO 920 R = 1, M                                                   23440095
         IF (IND(R) .NE. TAG) GO TO 920                                 23440096
         ITS = 1                                                        23440097
         X1 = W(R)                                                      23440098
         IF (S .NE. 0) GO TO 510                                        23440099
C     :::::::::: CHECK FOR ISOLATED ROOT ::::::::::                     23440100
         XU = 1.0D0                                                     23440101
         IF (P .NE. Q) GO TO 490                                        23440102
         RV6(P) = 1.0D0                                                 23440103
         GO TO 870                                                      23440104
  490    NORM = DABS(D(P))                                              23440105
         IP = P + 1                                                     23440106
C                                                                       23440107
         DO 500 I = IP, Q                                               23440108
  500    NORM = NORM + DABS(D(I)) + DABS(E(I))                          23440109
C     :::::::::: EPS2 IS THE CRITERION FOR GROUPING,                    23440110
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL                    23440111
C                ROOTS ARE MODIFIED BY EPS3,                            23440112
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ::::::::::  23440113
         EPS2 = 1.0D-3 * NORM                                           23440114
         EPS3 = MACHEP * NORM                                           23440115
         UK = DFLOAT(Q-P+1)                                             23440116
         EPS4 = UK * EPS3                                               23440117
         UK = EPS4 / DSQRT(UK)                                          23440118
         S = P                                                          23440119
  505    GROUP = 0                                                      23440120
         GO TO 520                                                      23440121
C     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS ::::::::::          23440122
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505                           23440123
         GROUP = GROUP + 1                                              23440124
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3       23440125
C     :::::::::: ELIMINATION WITH INTERCHANGES AND                      23440126
C                INITIALIZATION OF VECTOR ::::::::::                    23440127
  520    V = 0.0D0                                                      23440128
C                                                                       23440129
         DO 580 I = P, Q                                                23440130
            RV6(I) = UK                                                 23440131
            IF (I .EQ. P) GO TO 560                                     23440132
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540                      23440133
C     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF            23440134
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ::::::::::   23440135
            XU = U / E(I)                                               23440136
            RV4(I) = XU                                                 23440137
            RV1(I-1) = E(I)                                             23440138
            RV2(I-1) = D(I) - X1                                        23440139
            RV3(I-1) = 0.0D0                                            23440140
            IF (I .NE. Q) RV3(I-1) = E(I+1)                             23440141
            U = V - XU * RV2(I-1)                                       23440142
            V = -XU * RV3(I-1)                                          23440143
            GO TO 580                                                   23440144
  540       XU = E(I) / U                                               23440145
            RV4(I) = XU                                                 23440146
            RV1(I-1) = U                                                23440147
            RV2(I-1) = V                                                23440148
            RV3(I-1) = 0.0D0                                            23440149
  560       U = D(I) - X1 - XU * V                                      23440150
            IF (I .NE. Q) V = E(I+1)                                    23440151
  580    CONTINUE                                                       23440152
C                                                                       23440153
         IF (U .EQ. 0.0D0) U = EPS3                                     23440154
         RV1(Q) = U                                                     23440155
         RV2(Q) = 0.0D0                                                 23440156
         RV3(Q) = 0.0D0                                                 23440157
C     :::::::::: BACK SUBSTITUTION                                      23440158
C                FOR I=Q STEP -1 UNTIL P DO -- ::::::::::               23440159
  600    DO 620 II = P, Q                                               23440160
            I = P + Q - II                                              23440161
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)        23440162
            V = U                                                       23440163
            U = RV6(I)                                                  23440164
  620    CONTINUE                                                       23440165
C     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS                 23440166
C                MEMBERS OF GROUP ::::::::::                            23440167
         IF (GROUP .EQ. 0) GO TO 700                                    23440168
         J = R                                                          23440169
C                                                                       23440170
         DO 680 JJ = 1, GROUP                                           23440171
  630       J = J - 1                                                   23440172
            IF (IND(J) .NE. TAG) GO TO 630                              23440173
            XU = 0.0D0                                                  23440174
C                                                                       23440175
            DO 640 I = P, Q                                             23440176
  640       XU = XU + RV6(I) * Z(I,J)                                   23440177
C                                                                       23440178
            DO 660 I = P, Q                                             23440179
  660       RV6(I) = RV6(I) - XU * Z(I,J)                               23440180
C                                                                       23440181
  680    CONTINUE                                                       23440182
C                                                                       23440183
  700    NORM = 0.0D0                                                   23440184
C                                                                       23440185
         DO 720 I = P, Q                                                23440186
  720    NORM = NORM + DABS(RV6(I))                                     23440187
C                                                                       23440188
         IF (NORM .GE. 1.0D0) GO TO 840                                 23440189
C     :::::::::: FORWARD SUBSTITUTION ::::::::::                        23440190
         IF (ITS .EQ. 5) GO TO 830                                      23440191
         IF (NORM .NE. 0.0D0) GO TO 740                                 23440192
         RV6(S) = EPS4                                                  23440193
         S = S + 1                                                      23440194
         IF (S .GT. Q) S = P                                            23440195
         GO TO 780                                                      23440196
  740    XU = EPS4 / NORM                                               23440197
C                                                                       23440198
         DO 760 I = P, Q                                                23440199
  760    RV6(I) = RV6(I) * XU                                           23440200
C     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR                  23440201
C                ITERATE ::::::::::                                     23440202
  780    DO 820 I = IP, Q                                               23440203
            U = RV6(I)                                                  23440204
C     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE               23440205
C                WAS PERFORMED EARLIER IN THE                           23440206
C                TRIANGULARIZATION PROCESS ::::::::::                   23440207
            IF (RV1(I-1) .NE. E(I)) GO TO 800                           23440208
            U = RV6(I-1)                                                23440209
            RV6(I-1) = RV6(I)                                           23440210
  800       RV6(I) = U - RV4(I) * RV6(I-1)                              23440211
  820    CONTINUE                                                       23440212
C                                                                       23440213
         ITS = ITS + 1                                                  23440214
         GO TO 600                                                      23440215
C     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR ::::::::::      23440216
  830    IERR = -R                                                      23440217
         XU = 0.0D0                                                     23440218
         GO TO 870                                                      23440219
C     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS                    23440220
C                1 AND EXPAND TO FULL ORDER ::::::::::                  23440221
  840    U = 0.0D0                                                      23440222
C                                                                       23440223
         DO 860 I = P, Q                                                23440224
  860    U = U + RV6(I)**2                                              23440225
C                                                                       23440226
         XU = 1.0D0 / DSQRT(U)                                          23440227
C                                                                       23440228
  870    DO 880 I = 1, N                                                23440229
  880    Z(I,R) = 0.0D0                                                 23440230
C                                                                       23440231
         DO 900 I = P, Q                                                23440232
  900    Z(I,R) = RV6(I) * XU                                           23440233
C                                                                       23440234
         X0 = X1                                                        23440235
  920 CONTINUE                                                          23440236
C                                                                       23440237
      IF (Q .LT. N) GO TO 100                                           23440238
 1001 RETURN                                                            23440239
C     :::::::::: LAST CARD OF TINVIT ::::::::::                         23440240
      END                                                               23440241

C                                                                       79440001
C     ------------------------------------------------------------------79440002
C                                                                       79440003
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)                                   79440004
C                                                                       79440005
      INTEGER I,J,K,L,M,N,NM                                            79440006
      REAL*8 A(NM,N),E(N),Z(NM,M)                                       79440007
      REAL*8 S                                                          79440008
C                                                                       79440009
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,   79440010
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   79440011
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   79440012
C                                                                       79440013
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC        79440014
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING            79440015
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1.                79440016
C                                                                       79440017
C     ON INPUT:                                                         79440018
C                                                                       79440019
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         79440020
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          79440021
C          DIMENSION STATEMENT;                                         79440022
C                                                                       79440023
C        N IS THE ORDER OF THE MATRIX;                                  79440024
C                                                                       79440025
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             79440026
C          FORMATIONS USED IN THE REDUCTION BY  TRED1                   79440027
C          IN ITS STRICT LOWER TRIANGLE;                                79440028
C                                                                       79440029
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         79440030
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;        79440031
C                                                                       79440032
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED;        79440033
C                                                                       79440034
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED             79440035
C          IN ITS FIRST M COLUMNS.                                      79440036
C                                                                       79440037
C     ON OUTPUT:                                                        79440038
C                                                                       79440039
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS                        79440040
C          IN ITS FIRST M COLUMNS.                                      79440041
C                                                                       79440042
C     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.                79440043
C                                                                       79440044
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        79440045
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         79440046
C                                                                       79440047
C     ------------------------------------------------------------------79440048
C                                                                       79440049
      IF (M .EQ. 0) GO TO 200                                           79440050
      IF (N .EQ. 1) GO TO 200                                           79440051
C                                                                       79440052
      DO 140 I = 2, N                                                   79440053
         L = I - 1                                                      79440054
         IF (E(I) .EQ. 0.0D0) GO TO 140                                 79440055
C                                                                       79440056
         DO 130 J = 1, M                                                79440057
            S = 0.0D0                                                   79440058
C                                                                       79440059
            DO 110 K = 1, L                                             79440060
  110       S = S + A(I,K) * Z(K,J)                                     79440061
C     :::::::::: DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.        79440062
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ::::::::::   79440063
            S = (S / A(I,L)) / E(I)                                     79440064
C                                                                       79440065
            DO 120 K = 1, L                                             79440066
  120       Z(K,J) = Z(K,J) + S * A(I,K)                                79440067
C                                                                       79440068
  130    CONTINUE                                                       79440069
C                                                                       79440070
  140 CONTINUE                                                          79440071
C                                                                       79440072
  200 RETURN                                                            79440073
C     :::::::::: LAST CARD OF TRBAK1 ::::::::::                         79440074
      END                                                               79440075

      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   00000010
C                                                                       00000020
      INTEGER I,J,K,L,N,II,NM,JP1                                       00000030
C      REAL*8  A(NM,N),D(N),E(N),E2(N)                                  00000040
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)                          00000040
      DOUBLE PRECISION F,G,H,SCALE                                      00000050
C      REAL*8  F,G,H,SCALE                                              00000050
C                                                                       00000060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    00000070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   00000080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   00000090
C                                                                       00000100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   00000110
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           00000120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            00000130
C                                                                       00000140
C     ON INPUT                                                          00000150
C                                                                       00000160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         00000170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          00000180
C          DIMENSION STATEMENT.                                         00000190
C                                                                       00000200
C        N IS THE ORDER OF THE MATRIX.                                  00000210
C                                                                       00000220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          00000230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               00000240
C                                                                       00000250
C     ON OUTPUT                                                         00000260
C                                                                       00000270
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             00000280
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         00000290
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.        00000300
C                                                                       00000310
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    00000320
C                                                                       00000330
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         00000340
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      00000350
C                                                                       00000360
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    00000370
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        00000380
C                                                                       00000390
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    00000400
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 00000410
C                                                                       00000420
C     THIS VERSION DATED AUGUST 1983.                                   00000430
C                                                                       00000440
C     ------------------------------------------------------------------00000450
C                                                                       00000460
      DO 100 I = 1, N                                                   00000470
         D(I) = A(N,I)                                                  00000480
         A(N,I) = A(I,I)                                                00000490
  100 CONTINUE                                                          00000500
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........               00000510
      DO 300 II = 1, N                                                  00000520
         I = N + 1 - II                                                 00000530
         L = I - 1                                                      00000540
         H = 0.0D0                                                      00000550
         SCALE = 0.0D0                                                  00000560
         IF (L .LT. 1) GO TO 130                                        00000570
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       00000580
         DO 120 K = 1, L                                                00000590
  120    SCALE = SCALE + DABS(D(K))                                     00000600
C                                                                       00000610
         IF (SCALE .NE. 0.0D0) GO TO 140                                00000620
C                                                                       00000630
         DO 125 J = 1, L                                                00000640
            D(J) = A(L,J)                                               00000650
            A(L,J) = A(I,J)                                             00000660
            A(I,J) = 0.0D0                                              00000670
  125    CONTINUE                                                       00000680
C                                                                       00000690
  130    E(I) = 0.0D0                                                   00000700
         E2(I) = 0.0D0                                                  00000710
         GO TO 300                                                      00000720
C                                                                       00000730
  140    DO 150 K = 1, L                                                00000740
            D(K) = D(K) / SCALE                                         00000750
            H = H + D(K) * D(K)                                         00000760
  150    CONTINUE                                                       00000770
C                                                                       00000780
         E2(I) = SCALE * SCALE * H                                      00000790
         F = D(L)                                                       00000800
         G = -DSIGN(DSQRT(H),F)                                         00000810
         E(I) = SCALE * G                                               00000820
         H = H - F * G                                                  00000830
         D(L) = F - G                                                   00000840
         IF (L .EQ. 1) GO TO 285                                        00000850
C     .......... FORM A*U ..........                                    00000860
         DO 170 J = 1, L                                                00000870
  170    E(J) = 0.0D0                                                   00000880
C                                                                       00000890
         DO 240 J = 1, L                                                00000900
            F = D(J)                                                    00000910
            G = E(J) + A(J,J) * F                                       00000920
            JP1 = J + 1                                                 00000930
            IF (L .LT. JP1) GO TO 220                                   00000940
C                                                                       00000950
            DO 200 K = JP1, L                                           00000960
               G = G + A(K,J) * D(K)                                    00000970
               E(K) = E(K) + A(K,J) * F                                 00000980
  200       CONTINUE                                                    00000990
C                                                                       00001000
  220       E(J) = G                                                    00001010
  240    CONTINUE                                                       00001020
C     .......... FORM P ..........                                      00001030
         F = 0.0D0                                                      00001040
C                                                                       00001050
         DO 245 J = 1, L                                                00001060
            E(J) = E(J) / H                                             00001070
            F = F + E(J) * D(J)                                         00001080
  245    CONTINUE                                                       00001090
C                                                                       00001100
         H = F / (H + H)                                                00001110
C     .......... FORM Q ..........                                      00001120
         DO 250 J = 1, L                                                00001130
  250    E(J) = E(J) - H * D(J)                                         00001140
C     .......... FORM REDUCED A ..........                              00001150
         DO 280 J = 1, L                                                00001160
            F = D(J)                                                    00001170
            G = E(J)                                                    00001180
C                                                                       00001190
            DO 260 K = J, L                                             00001200
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)                       00001210
C                                                                       00001220
  280    CONTINUE                                                       00001230
C                                                                       00001240
  285    DO 290 J = 1, L                                                00001250
            F = D(J)                                                    00001260
            D(J) = A(L,J)                                               00001270
            A(L,J) = A(I,J)                                             00001280
            A(I,J) = F * SCALE                                          00001290
  290    CONTINUE                                                       00001300
C                                                                       00001310
  300 CONTINUE                                                          00001320
C                                                                       00001330
      RETURN                                                            00001340
      END                                                               00001350
*** FROM NETLIB, WED APR 27 21:20:12 CDT 1988 ***                       00001360
C     ------------------------------------------------------------------37720002
C                                                                       37720003
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)   37720004
C                                                                       37720005
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM      37720006
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)               37720007
C      REAL*8 D(N),E(N),E2(N),W(M),RV4(N),RV5(N)                        37720007
      REAL*8 U,V,LB,T1,T2,UB,XU,X0,X1,EPS1
      REAL*16 MACHEP                       
      REAL*8 DABS,DMAX1,DMIN1,DFLOAT                                    37720009
      INTEGER IND(M)                                                    37720010
      REAL*16 JUNK
C                                                                       37720011
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,   37720012
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.      37720013
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).   37720014
C                                                                       37720015
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL          37720016
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,              37720017
C     USING BISECTION.                                                  37720018
C                                                                       37720019
C     ON INPUT:                                                         37720020
C                                                                       37720021
C        N IS THE ORDER OF THE MATRIX;                                  37720022
C                                                                       37720023
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED           37720024
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,             37720025
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,           37720026
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE            37720027
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX;                   37720028
C                                                                       37720029
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          37720030
C                                                                       37720031
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        37720032
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;               37720033
C                                                                       37720034
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    37720035
C          E2(1) IS ARBITRARY;                                          37720036
C                                                                       37720037
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED         37720038
C          EIGENVALUES;                                                 37720039
C                                                                       37720040
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER      37720041
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.          37720042
C                                                                       37720043
C     ON OUTPUT:                                                        37720044
C                                                                       37720045
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS              37720046
C          (LAST) DEFAULT VALUE;                                        37720047
C                                                                       37720048
C        D AND E ARE UNALTERED;                                         37720049
C                                                                       37720050
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED        37720051
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE        37720052
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.            37720053
C          E2(1) IS ALSO SET TO ZERO;                                   37720054
C                                                                       37720055
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED    37720056
C          EIGENVALUES;                                                 37720057
C                                                                       37720058
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES          37720059
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER;              37720060
C                                                                       37720061
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES    37720062
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --        37720063
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM      37720064
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.;37720065
C                                                                       37720066
C        IERR IS SET TO                                                 37720067
C          ZERO       FOR NORMAL RETURN,                                37720068
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE         37720069
C                     UNIQUE SELECTION IMPOSSIBLE,                      37720070
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE         37720071
C                     UNIQUE SELECTION IMPOSSIBLE;                      37720072
C                                                                       37720073
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.                      37720074
C                                                                       37720075
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER  37720076
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.        37720077
C                                                                       37720078
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        37720079
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         37720080
C                                                                       37720081
C     ------------------------------------------------------------------37720082
C                                                                       37720083
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     37720084
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   37720085
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        37720086
C                ON S360 ::::::::::                                     37720087
C      DATA MACHEP/Z3410000000000000/                                    37720088
C       DATA MACHEP/1.D-17/
      JUNK = 16.0D0**(-13)
      MACHEP = JUNK
C      PRINT *,'MACHEP =',MACHEP
C                                                                       37720089
C      PRINT *,'IN THE CALL 2'
      IERR = 0                                                          37720090
      TAG = 0                                                           37720091
      XU = D(1)                                                         37720092
      X0 = D(1)                                                         37720093
      U = 0.0D0                                                         37720094
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN   37720095
C                INTERVAL CONTAINING ALL THE EIGENVALUES ::::::::::     37720096
C      PRINT *, N
      DO 40 I = 1, N                                                    37720097
C      PRINT *, 'IN LOOP 1'
         X1 = U                                                         37720098
         U = 0.0D0                                                      37720099
         IF (I .NE. N) U = DABS(E(I+1))                                 37720100
         XU = DMIN1(D(I)-(X1+U),XU)                                     37720101
         X0 = DMAX1(D(I)+(X1+U),X0)                                     37720102
         IF (I .EQ. 1) GO TO 20                                         37720103
         IF (DABS(E(I)) .GT. MACHEP * (DABS(D(I)) + DABS(D(I-1))))      37720104
     X      GO TO 40                                                    37720105
   20    E2(I) = 0.0D0                                                  37720106
   40 CONTINUE                                                          37720107
C                                                                       37720108
      X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP * DFLOAT(N)                37720109
      XU = XU - X1                                                      37720110
      T1 = XU                                                           37720111
      X0 = X0 + X1                                                      37720112
      T2 = X0                                                           37720113
C     :::::::::: DETERMINE AN INTERVAL CONTAINING EXACTLY               37720114
C                THE DESIRED EIGENVALUES ::::::::::                     37720115
      P = 1                                                             37720116
      Q = N                                                             37720117
      M1 = M11 - 1                                                      37720118
      IF (M1 .EQ. 0) GO TO 75                                           37720119
      ISTURM = 1                                                        37720120
   50 V = X1                                                            37720121
      X1 = XU + (X0 - XU) * 0.5D0                                       37720122
      IF (X1 .EQ. V) GO TO 980                                          37720123
      GO TO 320                                                         37720124
   60 IF (S - M1) 65, 73, 70                                            37720125
   65 XU = X1                                                           37720126
      GO TO 50                                                          37720127
   70 X0 = X1                                                           37720128
      GO TO 50                                                          37720129
   73 XU = X1                                                           37720130
      T1 = X1                                                           37720131
   75 M22 = M1 + M                                                      37720132
      IF (M22 .EQ. N) GO TO 90                                          37720133
      X0 = T2                                                           37720134
      ISTURM = 2                                                        37720135
      GO TO 50                                                          37720136
   80 IF (S - M22) 65, 85, 70                                           37720137
   85 T2 = X1                                                           37720138
   90 Q = 0                                                             37720139
      R = 0                                                             37720140
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING         37720141
C                INTERVAL BY THE GERSCHGORIN BOUNDS ::::::::::          37720142
  100 IF (R .EQ. M) GO TO 1001                                          37720143
      TAG = TAG + 1                                                     37720144
      P = Q + 1                                                         37720145
      XU = D(P)                                                         37720146
      X0 = D(P)                                                         37720147
      U = 0.0D0                                                         37720148
C                                                                       37720149
      DO 120 Q = P, N                                                   37720150
C      PRINT *, 'IN LOOP 2'
         X1 = U                                                         37720151
         U = 0.0D0                                                      37720152
         V = 0.0D0                                                      37720153
         IF (Q .EQ. N) GO TO 110                                        37720154
         U = DABS(E(Q+1))                                               37720155
         V = E2(Q+1)                                                    37720156
  110    XU = DMIN1(D(Q)-(X1+U),XU)                                     37720157
         X0 = DMAX1(D(Q)+(X1+U),X0)                                     37720158
         IF (V .EQ. 0.0D0) GO TO 140                                    37720159
  120 CONTINUE                                                          37720160
C                                                                       37720161
  140 X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP                            37720162
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1                                   37720163
      IF (P .NE. Q) GO TO 180                                           37720164
C     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL ::::::::::     37720165
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940                     37720166
      M1 = P                                                            37720167
      M2 = P                                                            37720168
      RV5(P) = D(P)                                                     37720169
      GO TO 900                                                         37720170
  180 X1 = X1 * DFLOAT(Q-P+1)                                           37720171
      LB = DMAX1(T1,XU-X1)                                              37720172
      UB = DMIN1(T2,X0+X1)                                              37720173
      X1 = LB                                                           37720174
      ISTURM = 3                                                        37720175
      GO TO 320                                                         37720176
  200 M1 = S + 1                                                        37720177
      X1 = UB                                                           37720178
      ISTURM = 4                                                        37720179
      GO TO 320                                                         37720180
  220 M2 = S                                                            37720181
      IF (M1 .GT. M2) GO TO 940                                         37720182
C     :::::::::: FIND ROOTS BY BISECTION ::::::::::                     37720183
      X0 = UB                                                           37720184
      ISTURM = 5                                                        37720185
C                                                                       37720186
      DO 240 I = M1, M2                                                 37720187
C      PRINT *, 'IN LOOP 3'
         RV5(I) = UB                                                    37720188
         RV4(I) = LB                                                    37720189
  240 CONTINUE                                                          37720190
C     :::::::::: LOOP FOR K-TH EIGENVALUE                               37720191
C                FOR K=M2 STEP -1 UNTIL M1 DO --                        37720192
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ::::::::::37720193
      K = M2                                                            37720194
  250    XU = LB                                                        37720195
C     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- ::::::::::              37720196
         DO 260 II = M1, K                                              37720197
C      PRINT *, 'IN LOOP 4'
            I = M1 + K - II                                             37720198
            IF (XU .GE. RV4(I)) GO TO 260                               37720199
            XU = RV4(I)                                                 37720200
            GO TO 280                                                   37720201
  260    CONTINUE                                                       37720202
C                                                                       37720203
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)                                37720204
C     :::::::::: NEXT BISECTION STEP ::::::::::                         37720205
  300    X1 = (XU + X0) * 0.5D0                                         37720206
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP *                           37720207
     X      (DABS(XU) + DABS(X0)) + DABS(EPS1))) GO TO 420              37720208
C     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE ::::::::::        37720209
  320    S = P - 1                                                      37720210
         U = 1.0D0                                                      37720211
C                                                                       37720212
C      PRINT *, P,Q
         DO 340 I = P, Q                                                37720213
C      PRINT *, 'IN LOOP 5'
C      PRINT *, I
            IF (U .NE. 0.0D0) GO TO 325                                 37720214
            V = DABS(E(I)) / MACHEP                                     37720215
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0                             37720216
            GO TO 330                                                   37720217
  325       V = E2(I) / U                                               37720218
  330       U = D(I) - X1 - V                                           37720219
            IF (U .LT. 0.0D0) S = S + 1                                 37720220
  340    CONTINUE                                                       37720221
C                                                                       37720222
C      PRINT *,'ISTURM =', ISTURM
         GO TO (60,80,200,220,360), ISTURM                              37720223
C     :::::::::: REFINE INTERVALS ::::::::::                            37720224
  360    IF (S .GE. K) GO TO 400                                        37720225
         XU = X1                                                        37720226
         IF (S .GE. M1) GO TO 380                                       37720227
         RV4(M1) = X1                                                   37720228
         GO TO 300                                                      37720229
  380    RV4(S+1) = X1                                                  37720230
         IF (RV5(S) .GT. X1) RV5(S) = X1                                37720231
         GO TO 300                                                      37720232
  400    X0 = X1                                                        37720233
         GO TO 300                                                      37720234
C     :::::::::: K-TH EIGENVALUE FOUND ::::::::::                       37720235
  420    RV5(K) = X1                                                    37720236
      K = K - 1                                                         37720237
      IF (K .GE. M1) GO TO 250                                          37720238
C     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR                    37720239
C                SUBMATRIX ASSOCIATIONS ::::::::::                      37720240
  900 S = R                                                             37720241
      R = R + M2 - M1 + 1                                               37720242
      J = 1                                                             37720243
      K = M1                                                            37720244
C                                                                       37720245
      DO 920 L = 1, R                                                   37720246
C      PRINT *, 'IN LOOP 6'
         IF (J .GT. S) GO TO 910                                        37720247
         IF (K .GT. M2) GO TO 940                                       37720248
         IF (RV5(K) .GE. W(L)) GO TO 915                                37720249
C                                                                       37720250
         DO 905 II = J, S                                               37720251
C      PRINT *, 'IN LOOP 7'
            I = L + S - II                                              37720252
            W(I+1) = W(I)                                               37720253
            IND(I+1) = IND(I)                                           37720254
  905    CONTINUE                                                       37720255
C                                                                       37720256
  910    W(L) = RV5(K)                                                  37720257
         IND(L) = TAG                                                   37720258
         K = K + 1                                                      37720259
         GO TO 920                                                      37720260
  915    J = J + 1                                                      37720261
  920 CONTINUE                                                          37720262
C                                                                       37720263
  940 IF (Q .LT. N) GO TO 100                                           37720264
      GO TO 1001                                                        37720265
C     :::::::::: SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING       37720266
C                EXACTLY THE DESIRED EIGENVALUES ::::::::::             37720267
  980 IERR = 3 * N + ISTURM                                             37720268
 1001 LB = T1                                                           37720269
      UB = T2                                                           37720270
      RETURN                                                            37720271
C     :::::::::: LAST CARD OF TRIDIB ::::::::::                         37720272
      END                                                               37720273

