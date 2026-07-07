      SUBROUTINE   EIGVEC(NM,N,A,WR,WI,Z)                               
      INTEGER      NM,N,LOW,IGH,IERR                                    
      REAL         A(NM,N),Z(NM,N), WR(N),WI(N)                         
      REAL         SCALE(300)                                           
COMMENT*XXXX;                                                           
C**PURPOSE     FINDS ALL N EIGENVALUES AND NORMALIZED EIGENVECTORS OF   
C   THE GENERAL REAL MATRIX A.                                          
COMMENT*YYYY;                                                           
C INPUT       NM - ROW DIMENSION OF THE TW0-DIM. ARRAYS A AND Z AS      
C                  DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT  
C             N  - ORDER OF MATRIX                                      
C             A  - CONTAINS THE GENERAL REAL INPUT MATRIX.              
C                                                                       
C OUTPUT      A  - THE INPUT MATRIX IS DESTROYED BY THE COMPUTATION     
C             WR - THE REAL PARTS OF THE EIGENVALUES ARE STORED IN THE  
C                  ARRAY WR                                             
C             WI - THE IMAGINARY PARTS OF THE EIGENVALUES ARE STORED    
C                  IN THE ARRAY WI.  COMPLEX CONJUGATE EIGENVALUES      
C                  ARE STORED CONSECUTIVELY.                            
C             Z  - TW0-DIMENSIONAL ARRAY CONTAINING THE NORMALIZED      
C                  EIGENVECTORS OF A IN ITS FIRST N COLUMNS: IF THE     
C                  K-TH EIGENVALUE IS REAL, COLUMN K CONTAINS ITS       
C                  EIGENVECTOR; IF THE K-TH EIGENVALUE IS COMPLEX WITH  
C                  POSITIVE IMAGINARY PART, COLUMN K CONTAINS THE REAL  
C                  PART AND COLUMN K+1 THE IMAGINARY PART OF THE KTH    
C                  EIGENVECTOR. THE COMPLEX CONJUGATE EIGENVECTOR IS    
C                  NOT STORED.                                          
C                                                                       
C SUBROUTINES REQUIRED:  BALANC, ORTHES, ORTRAN, HQR2, BALBAK, NRMLZ    
C                                                                       
C RESTRICTION  THE ORDER N OF THE MATRIX MUST NOT EXCEED 100.  FOR      
C   LARGER MATRICES, THE DIMENSION OF THE WORKING STORAGE ARRAY         
C   SCALE MUST BE INCREASED.  THE ARRAYS A AND Z  MUST HAVE THE SAME    
C   ACTUAL ROW DIMENSION AND MUST OCCUPY SEPARATE STORAGE.              
C                                                                       
C METHOD      A IS FIRST PREPARED BY SUBROUTINE BALANC, WHICH DETECTS   
C   ANY ISOLATED EIGENVALUES AND REDUCES THE NORM OF A WITHOUT INTRO-   
C   DUCING ANY ROUNDOFF ERRORS.  ORTHES IS THEN CALLED TO REDUCE A TO   
C   UPPER HESSENBERG FORM USING ELEMENTARY ORTHOGONAL (HOUSEHOLDER)     
C   TRANSFORMATIONS.  THE MATRIX DEFINING THE SIMILARITY REDUCTION OF   
C   A TO UPPER HESSENBERG FORM IS THEN FORMED AND STORED IN THE ARRAY   
C   U BY SUBROUTINE ORTRANS USING THE INFORMATION LEFT BELOW THE SUB-   
C   DIAGONAL OF A BY ORTHES, COMPLETING THE PRELIMINARY REDUCTIONS AND  
C   TRANSFORMATIONS.  SUBROUTINE HQR2 IS EMPLOYED TO OBTAIN THE EIGEN-  
C   VALUES OF A AND THE EIGENVECTORS OF THE BALANCED MATRIX.  HQR2      
C   UTILIZES THE FRANCIS DOUBLE-STEP QR ALGORITHM TO COMPUTE THE        
C   EIGENVALUES AND A METHOD OF BACK SUBSTITUTION FOR THE EIGENVECTORS. 
C   THE EIGENVECTORS OF THE BALANCED MATRIX ARE NEXT TRANSFORMED INTO   
C   EIGENVECTORS OF THE ORIGINAL MATRIX BY SUBROUTINE BALBAK.  FINALLY, 
C   THE EIGENVECTORS OF A ARE NORMALIZED SO THAT A REAL EIGENVECTOR     
C   HAS ITS COMPONENT OF MAXIMUM MAGNITUDE EQUAL TO 1.0, AND A COMPLEX  
C   EIGENVECTOR HAS ITS MAXIMUM MODULUS COMPONENT EQUAL TO 1.0 + 0*I.   
C                                                                       
C REMARKS     BALANCING IS STRONGLY RECOMMENDED SINCE IT GENERALLY      
C   IMPROVES ACCURACY AT A SMALL ADDITIONAL COST IN TIME.  HOWEVER, IT  
C   IS NOT NECESSARY TO BALANCE ORTHOGONAL OR NORMAL MATRICES.  IF      
C   BALANC IS NOT USED, LOW MUST BE SET TO 1 AND UPP TO N.              
C             IF THE ORIGINAL MATRIX IS ALREADY OF UPPER HESSENBERG     
C   FORM (ZEROES BELOW THE FIRST SUB-DIAGONAL), THE CALLS TO ORTHES     
C   AND ORTRANS MAY BE REMOVED (SUBROUTINE BALANC WILL PRESERVE THE     
C   HESSENBERG FORM).  IN THIS CASE U MUST BE SET TO THE IDENTITY       
C   MATRIX BEFORE HQR2 IS ENTERED.                                      
C             IF UNORMALIZED EIGENVECTORS ARE SUFFICIENT, THE CALL TO   
C   SUBROUTINE NORMALIZE MAY BE REMOVED.                                
C             IF THE K-TH EIGENVALUE HAS NOT BEEN DETERMINED AFTER      
C   30 ITERATIONS, THE VALUE OF K WILL BE PRINTED AND THE SUBROUTINE    
C   EXITED NORMALLY WITH THE EIGENVALUES OF INDICES K+1, K+2, ..., N    
C   STORED CORRECTLY.  FAILURES WILL BE EXCEEDINGLY RARE.               
C                                                                       
C REFERENCE   FOR MORE INFORMATION REFER TO G. PETERS AND J. H.         
C   WILKINSON, "EIGENVECTORS OF REAL AND COMPLEX MATRICES BY LR AND     
C   QR TRIANGULARIZATIONS," NUMERISCHE MATHEMATIK 16, 1970, PP 181-204. 
C                                                                       
C ACKNOWLEDGEMNT  THESE FORTRAN TRANSLATIONS OF BALANC, ORTHES,         
C   ORTRAN, HQR2 AND BALBAK ARE PART OF THE EISPACK SERIES OF           
C   SUBROUTINES DEVELOPED AND MADE AVAILABLE BY  NATS PROJECT OF THE    
C   THE APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORIES,    
C   DIVISION, ARGONNE NATIONAL LABORATORIES, ARGONNE, ILLINOIS.         
C                                                                       
COMMENT*ZZZZ;                                                           
C     print *,'n = ',n
C     print *,'nm = ',nm
C     print *,'A:'
C     do i=1,n
C      do j=1,n
C       print '(2i,e)',i,j,a(i,j)
C      enddo
C     enddo
      CALL BALANC(NM,N,A,LOW,IGH,SCALE)                                 
      CALL ORTHES(NM,N,LOW,IGH,A,WR)                                    
C     ARRAY WR IS USED TEMPORARILY TO PASS INFO FROM ORTHES TO ORTRAN   
      CALL ORTRAN(NM,N,LOW,IGH,A,WR,Z)                                  
      CALL HQR2(NM,N,LOW,IGH,A,WR,WI,Z,IERR)                            
      CALL BALBAK(NM,N,LOW,IGH,SCALE,N,Z)                               
      CALL NRMLZ(NM,N,WI,Z)                                             
      IF (IERR .NE. 0) PRINT 1,IERR                                     
    1 FORMAT(I3,'-TH EIGENVALUE WAS NOT FOUND AFTER 3O ITERATIONS')     
C     do i=1,n
C      do j=1,n
C       print '(2i,e)',i,j,z(i,j)
C      enddo
C     enddo
      RETURN                                                            
      END                                                               
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)                           
                                                                        
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC                            
      REAL A(NM,N),SCALE(N)                                             
      REAL C,F,G,R,S,B2,RADIX                                           
      REAL ABS                                                          
      LOGICAL NOCONV                                                    
COMMENT*XXXX;                                                           
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE BALANCE,             
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.              
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).   
C                                                                       
C     BALANC   BALANCES A REAL MATRIX AND ISOLATES                      
C     EIGENVALUES WHENEVER POSSIBLE.                                    
COMMENT*YYYY;                                                           
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.                    
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        A CONTAINS THE BALANCED MATRIX,                                
C                                                                       
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)                  
C          IS EQUAL TO ZERO IF                                          
C           (1) I IS GREATER THAN J AND                                 
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,                         
C                                                                       
C        SCALE CONTAINS INFORMATION DETERMINING THE                     
C           PERMUTATIONS AND SCALING FACTORS USED.                      
C                                                                       
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH      
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED       
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS         
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN          
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1                        
C                 = D(J,J),      J = LOW,...,IGH                        
C                 = P(J)         J = IGH+1,...,N.                       
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,       
C     THEN 1 TO LOW-1.                                                  
C                                                                       
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.          
C                                                                       
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN           
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS       
C     K,L HAVE BEEN REVERSED.)                                          
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING      
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION. 
C                                                                       
C                **********                                             
COMMENT*ZZZZ;                                                           
      RADIX = 8.                                                        
C                                                                       
      B2 = RADIX * RADIX                                                
      K = 1                                                             
      L = N                                                             
      GO TO 100                                                         
C     ********** IN-LINE PROCEDURE FOR ROW AND                          
C                COLUMN EXCHANGE **********                             
   20 SCALE(M) = J                                                      
      IF (J .EQ. M) GO TO 50                                            
C                                                                       
      DO 30 I = 1, L                                                    
         F = A(I,J)                                                     
         A(I,J) = A(I,M)                                                
         A(I,M) = F                                                     
   30 CONTINUE                                                          
C                                                                       
      DO 40 I = K, N                                                    
         F = A(J,I)                                                     
         A(J,I) = A(M,I)                                                
         A(M,I) = F                                                     
   40 CONTINUE                                                          
C                                                                       
   50 GO TO (80,130), IEXC                                              
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE                
C                AND PUSH THEM DOWN **********                          
   80 IF (L .EQ. 1) GO TO 280                                           
      L = L - 1                                                         
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********               
  100 DO 120 JJ = 1, L                                                  
         J = L + 1 - JJ                                                 
C                                                                       
         DO 110 I = 1, L                                                
            IF (I .EQ. J) GO TO 110                                     
            IF (A(J,I) .NE. 0.0) GO TO 120                              
  110    CONTINUE                                                       
C                                                                       
         M = L                                                          
         IEXC = 1                                                       
         GO TO 20                                                       
  120 CONTINUE                                                          
C                                                                       
      GO TO 140                                                         
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE             
C                AND PUSH THEM LEFT **********                          
  130 K = K + 1                                                         
C                                                                       
  140 DO 170 J = K, L                                                   
C                                                                       
         DO 150 I = K, L                                                
            IF (I .EQ. J) GO TO 150                                     
            IF (A(I,J) .NE. 0.0) GO TO 170                              
  150    CONTINUE                                                       
C                                                                       
         M = K                                                          
         IEXC = 2                                                       
         GO TO 20                                                       
  170 CONTINUE                                                          
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********    
      DO 180 I = K, L                                                   
  180 SCALE(I) = 1.0                                                    
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********           
  190 NOCONV = .FALSE.                                                  
C                                                                       
      DO 270 I = K, L                                                   
         C = 0.0                                                        
         R = 0.0                                                        
C                                                                       
         DO 200 J = K, L                                                
            IF (J .EQ. I) GO TO 200                                     
            C = C + ABS(A(J,I))                                         
            R = R + ABS(A(I,J))                                         
  200    CONTINUE                                                       
C                                                                       
         G = R / RADIX                                                  
         F = 1.0                                                        
         S = C + R                                                      
  210    IF (C .GE. G) GO TO 220                                        
         F = F * RADIX                                                  
         C = C * B2                                                     
         GO TO 210                                                      
  220    G = R * RADIX                                                  
  230    IF (C .LT. G) GO TO 240                                        
         F = F / RADIX                                                  
         C = C / B2                                                     
         GO TO 230                                                      
C     ********** NOW BALANCE **********                                 
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270                       
         G = 1.0 / F                                                    
         SCALE(I) = SCALE(I) * F                                        
         NOCONV = .TRUE.                                                
C                                                                       
         DO 250 J = K, N                                                
  250    A(I,J) = A(I,J) * G                                            
C                                                                       
         DO 260 J = 1, L                                                
  260    A(J,I) = A(J,I) * F                                            
C                                                                       
  270 CONTINUE                                                          
C                                                                       
      IF (NOCONV) GO TO 190                                             
C                                                                       
  280 LOW = K                                                           
      IGH = L                                                           
      RETURN                                                            
C     ********** LAST CARD OF BALANC **********                         
      END                                                               
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)                             
                                                                        
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW                        
      REAL A(NM,N),ORT(IGH)                                             
      REAL F,G,H,SCALE                                                  
      REAL SQRT,ABS,SIGN                                                
COMMENT*XXXX;                                                           
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE ORTHES,              
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.             
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).   
C                                                                       
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE                      
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS                  
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY                       
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            
COMMENT*YYYY;                                                           
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING           
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,          
C          SET LOW=1, IGH=N,                                            
C                                                                       
C        A CONTAINS THE INPUT MATRIX.                                   
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT           
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION         
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE                
C          HESSENBERG MATRIX,                                           
C                                                                       
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.    
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.                      
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
COMMENT*ZZZZ;                                                           
      LA = IGH - 1                                                      
      KP1 = LOW + 1                                                     
      IF (LA .LT. KP1) GO TO 200                                        
C                                                                       
      DO 180 M = KP1, LA                                                
         H = 0.0                                                        
         ORT(M) = 0.0                                                   
         SCALE = 0.0                                                    
C     ********** SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) **********    
         DO 90 I = M, IGH                                               
   90    SCALE = SCALE + ABS(A(I,M-1))                                  
C                                                                       
         IF (SCALE .EQ. 0.0) GO TO 180                                  
         MP = M + IGH                                                   
C     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********             
         DO 100 II = M, IGH                                             
            I = MP - II                                                 
            ORT(I) = A(I,M-1) / SCALE                                   
            H = H + ORT(I) * ORT(I)                                     
  100    CONTINUE                                                       
C                                                                       
         G = -SIGN(SQRT(H),ORT(M))                                      
         H = H - ORT(M) * G                                             
         ORT(M) = ORT(M) - G                                            
C     ********** FORM (I-(U*UT)/H) * A **********                       
         DO 130 J = M, N                                                
            F = 0.0                                                     
C     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********             
            DO 110 II = M, IGH                                          
               I = MP - II                                              
               F = F + ORT(I) * A(I,J)                                  
  110       CONTINUE                                                    
C                                                                       
            F = F / H                                                   
C                                                                       
            DO 120 I = M, IGH                                           
  120       A(I,J) = A(I,J) - F * ORT(I)                                
C                                                                       
  130    CONTINUE                                                       
C     ********** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) **********            
         DO 160 I = 1, IGH                                              
            F = 0.0                                                     
C     ********** FOR J=IGH STEP -1 UNTIL M DO -- **********             
            DO 140 JJ = M, IGH                                          
               J = MP - JJ                                              
               F = F + ORT(J) * A(I,J)                                  
  140       CONTINUE                                                    
C                                                                       
            F = F / H                                                   
C                                                                       
            DO 150 J = M, IGH                                           
  150       A(I,J) = A(I,J) - F * ORT(J)                                
C                                                                       
  160    CONTINUE                                                       
C                                                                       
         ORT(M) = SCALE * ORT(M)                                        
         A(M,M-1) = SCALE * G                                           
  180 CONTINUE                                                          
C                                                                       
  200 RETURN                                                            
C     ********** LAST CARD OF ORTHES **********                         
      END                                                               
      SUBROUTINE ORTRAN(NM,N,LOW,IGH,A,ORT,Z)                           
                                                                        
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1                             
      REAL A(NM,IGH),ORT(IGH),Z(NM,N)                                   
      REAL G,H                                                          
COMMENT*XXXX;                                                           
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,             
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.             
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).   
C                                                                       
C     ORTRAN   ACCUMULATES THE ORTHOGONAL SIMILARITY                    
C     TRANSFORMATIONS USED IN THE REDUCTION OF A REAL GENERAL           
C     MATRIX TO UPPER HESSENBERG FORM BY  ORTHES.                       
COMMENT*YYYY;                                                           
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING           
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,          
C          SET LOW=1, IGH=N,                                            
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES                  
C          IN ITS STRICT LOWER TRIANGLE,                                
C                                                                       
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-              
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES.                 
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.                      
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C          REDUCTION BY  ORTHES,                                        
C                                                                       
C        ORT HAS BEEN ALTERED.                                          
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
COMMENT*ZZZZ;                                                           
C     ********** INITIALIZE Z TO IDENTITY MATRIX **********             
      DO 80 I = 1, N                                                    
C                                                                       
         DO 60 J = 1, N                                                 
   60    Z(I,J) = 0.0                                                   
C                                                                       
         Z(I,I) = 1.0                                                   
   80 CONTINUE                                                          
C                                                                       
      KL = IGH - LOW - 1                                                
      IF (KL .LT. 1) GO TO 200                                          
C     ********** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********      
      DO 140 MM = 1, KL                                                 
         MP = IGH - MM                                                  
C     ********** H BELOW IS NEGATIVE OF H FORMED IN ORTHES **********   
         H = A(MP,MP-1) * ORT(MP)                                       
         IF (H .EQ. 0.0) GO TO 140                                      
         MP1 = MP + 1                                                   
C                                                                       
         DO 100 I = MP1, IGH                                            
  100    ORT(I) = A(I,MP-1)                                             
C                                                                       
         DO 130 J = MP, IGH                                             
            G = 0.0                                                     
C                                                                       
            DO 110 I = MP, IGH                                          
  110       G = G + ORT(I) * Z(I,J)                                     
C                                                                       
            G = G / H                                                   
C                                                                       
            DO 120 I = MP, IGH                                          
  120       Z(I,J) = Z(I,J) + G * ORT(I)                                
C                                                                       
  130    CONTINUE                                                       
C                                                                       
  140 CONTINUE                                                          
C                                                                       
  200 RETURN                                                            
C     ********** LAST CARD OF ORTRAN **********                         
      END                                                               
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)                      
                                                                        
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,                      
     *IGH,ITS,LOW,MP2,ENM2,IERR                                 ORZ15720
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)                                  
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP                   
      REAL SQRT,ABS,SIGN                                                
      INTEGER MIN0                                                      
      LOGICAL NOTLAS                                                    
      COMPLEX Z3                                                        
      COMPLEX CMPLX                                                     
      REAL T3(2)                                                        
      EQUIVALENCE (Z3,T3(1))                                            
COMMENT*XXXX;                                                           
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE HQR2,                
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.             
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).   
C                                                                       
C     HQR2   FINDS THE EIGENVALUES AND EIGENVECTORS                     
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE          
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND           
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE            
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM        
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.                 
COMMENT*YYYY;                                                           
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING           
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,          
C          SET LOW=1, IGH=N,                                            
C                                                                       
C        H CONTAINS THE UPPER HESSENBERG MATRIX,                        
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN       
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE     
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS     
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE     
C          IDENTITY MATRIX.                                             
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        H HAS BEEN DESTROYED,                                          
C                                                                       
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,                
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES           
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS            
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE           
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN             
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT        
C          FOR INDICES IERR+1,...,N,                                    
C                                                                       
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.   
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z         
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX 
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH          
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS     
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN      
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND, 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL        
C     PROCEDURE CDIV BY COMPLEX DIVISION.                               
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
COMMENT*ZZZZ;                                                           
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = .5*8.0**(-12)                                            
C                                                                       
      IERR = 0                                                          
C     ********** STORE ROOTS ISOLATED BY BALANC **********              
      DO 50 I = 1, N                                                    
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50                      
         WR(I) = H(I,I)                                                 
         WI(I) = 0.0                                                    
   50 CONTINUE                                                          
C                                                                       
      EN = IGH                                                          
      T = 0.0                                                           
C     ********** SEARCH FOR NEXT EIGENVALUES **********                 
   60 IF (EN .LT. LOW) GO TO 340                                        
      ITS = 0                                                           
      NA = EN - 1                                                       
      ENM2 = NA - 1                                                     
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT             
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********            
   70 DO 80 LL = LOW, EN                                                
         L = EN + LOW - LL                                              
         IF (L .EQ. LOW) GO TO 100                                      
         IF (ABS(H(L,L-1)) .LE. MACHEP * (ABS(H(L-1,L-1))               
     X      + ABS(H(L,L)))) GO TO 100                                   
   80 CONTINUE                                                          
C     ********** FORM SHIFT **********                                  
  100 X = H(EN,EN)                                                      
      IF (L .EQ. EN) GO TO 270                                          
      Y = H(NA,NA)                                                      
      W = H(EN,NA) * H(NA,EN)                                           
      IF (L .EQ. NA) GO TO 280                                          
      IF (ITS .EQ. 30) GO TO 1000                                       
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130                      
C     ********** FORM EXCEPTIONAL SHIFT **********                      
      T = T + X                                                         
C                                                                       
      DO 120 I = LOW, EN                                                
  120 H(I,I) = H(I,I) - X                                               
C                                                                       
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))                               
      X = 0.75 * S                                                      
      Y = X                                                             
      W = -0.4375 * S * S                                               
  130 ITS = ITS + 1                                                     
C     ********** LOOK FOR TWO CONSECUTIVE SMALL                         
C                SUB-DIAGONAL ELEMENTS.                                 
C                FOR M=EN-2 STEP -1 UNTIL L DO -- **********            
      DO 140 MM = L, ENM2                                               
         M = ENM2 + L - MM                                              
         ZZ = H(M,M)                                                    
         R = X - ZZ                                                     
         S = Y - ZZ                                                     
         P = (R * S - W) / H(M+1,M) + H(M,M+1)                          
         Q = H(M+1,M+1) - ZZ - R - S                                    
         R = H(M+2,M+1)                                                 
         S = ABS(P) + ABS(Q) + ABS(R)                                   
         P = P / S                                                      
         Q = Q / S                                                      
         R = R / S                                                      
         IF (M .EQ. L) GO TO 150                                        
         IF (ABS(H(M,M-1)) * (ABS(Q) + ABS(R)) .LE. MACHEP * ABS(P)     
     X    * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))) GO TO 150    
  140 CONTINUE                                                          
C                                                                       
  150 MP2 = M + 2                                                       
C                                                                       
      DO 160 I = MP2, EN                                                
         H(I,I-2) = 0.0                                                 
         IF (I .EQ. MP2) GO TO 160                                      
         H(I,I-3) = 0.0                                                 
  160 CONTINUE                                                          
C     ********** DOUBLE QR STEP INVOLVING ROWS L TO EN AND              
C                COLUMNS M TO EN **********                             
      DO 260 K = M, NA                                                  
         NOTLAS = K .NE. NA                                             
         IF (K .EQ. M) GO TO 170                                        
         P = H(K,K-1)                                                   
         Q = H(K+1,K-1)                                                 
         R = 0.0                                                        
         IF (NOTLAS) R = H(K+2,K-1)                                     
         X = ABS(P) + ABS(Q) + ABS(R)                                   
         IF (X .EQ. 0.0) GO TO 260                                      
         P = P / X                                                      
         Q = Q / X                                                      
         R = R / X                                                      
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)                                  
         IF (K .EQ. M) GO TO 180                                        
         H(K,K-1) = -S * X                                              
         GO TO 190                                                      
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)                             
  190    P = P + S                                                      
         X = P / S                                                      
         Y = Q / S                                                      
         ZZ = R / S                                                     
         Q = Q / P                                                      
         R = R / P                                                      
C     ********** ROW MODIFICATION **********                            
         DO 210 J = K, N                                                
            P = H(K,J) + Q * H(K+1,J)                                   
            IF (.NOT. NOTLAS) GO TO 200                                 
            P = P + R * H(K+2,J)                                        
            H(K+2,J) = H(K+2,J) - P * ZZ                                
  200       H(K+1,J) = H(K+1,J) - P * Y                                 
            H(K,J) = H(K,J) - P * X                                     
  210    CONTINUE                                                       
C                                                                       
         J = MIN0(EN,K+3)                                               
C     ********** COLUMN MODIFICATION **********                         
         DO 230 I = 1, J                                                
            P = X * H(I,K) + Y * H(I,K+1)                               
            IF (.NOT. NOTLAS) GO TO 220                                 
            P = P + ZZ * H(I,K+2)                                       
            H(I,K+2) = H(I,K+2) - P * R                                 
  220       H(I,K+1) = H(I,K+1) - P * Q                                 
            H(I,K) = H(I,K) - P                                         
  230    CONTINUE                                                       
C     ********** ACCUMULATE TRANSFORMATIONS **********                  
         DO 250 I = LOW, IGH                                            
            P = X * Z(I,K) + Y * Z(I,K+1)                               
            IF (.NOT. NOTLAS) GO TO 240                                 
            P = P + ZZ * Z(I,K+2)                                       
            Z(I,K+2) = Z(I,K+2) - P * R                                 
  240       Z(I,K+1) = Z(I,K+1) - P * Q                                 
            Z(I,K) = Z(I,K) - P                                         
  250    CONTINUE                                                       
C                                                                       
  260 CONTINUE                                                          
C                                                                       
      GO TO 70                                                          
C     ********** ONE ROOT FOUND **********                              
  270 H(EN,EN) = X  +  T                                                
      WR(EN) = H(EN,EN)                                                 
      WI(EN) = 0.0                                                      
      EN = NA                                                           
      GO TO 60                                                          
C     ********** TWO ROOTS FOUND **********                             
  280 P = (Y - X) / 2.0                                                 
      Q = P * P + W                                                     
      ZZ = SQRT(ABS(Q))                                                 
      H(EN,EN) = X + T                                                  
      X = H(EN,EN)                                                      
      H(NA,NA) = Y + T                                                  
      IF (Q .LT. 0.0) GO TO 320                                         
C     ********** REAL PAIR **********                                   
      ZZ = P + SIGN(ZZ,P)                                               
      WR(NA) = X + ZZ                                                   
      WR(EN) = WR(NA)                                                   
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ                              
      WI(NA) = 0.0                                                      
      WI(EN) = 0.0                                                      
      X = H(EN,NA)                                                      
      R = SQRT(X*X+ZZ*ZZ)                                               
      P = X / R                                                         
      Q = ZZ / R                                                        
C     ********** ROW MODIFICATION **********                            
      DO 290 J = NA, N                                                  
         ZZ = H(NA,J)                                                   
         H(NA,J) = Q * ZZ + P * H(EN,J)                                 
         H(EN,J) = Q * H(EN,J) - P * ZZ                                 
  290 CONTINUE                                                          
C     ********** COLUMN MODIFICATION **********                         
      DO 300 I = 1, EN                                                  
         ZZ = H(I,NA)                                                   
         H(I,NA) = Q * ZZ + P * H(I,EN)                                 
         H(I,EN) = Q * H(I,EN) - P * ZZ                                 
  300 CONTINUE                                                          
C     ********** ACCUMULATE TRANSFORMATIONS **********                  
      DO 310 I = LOW, IGH                                               
         ZZ = Z(I,NA)                                                   
         Z(I,NA) = Q * ZZ + P * Z(I,EN)                                 
         Z(I,EN) = Q * Z(I,EN) - P * ZZ                                 
  310 CONTINUE                                                          
C                                                                       
      GO TO 330                                                         
C     ********** COMPLEX PAIR **********                                
  320 WR(NA) = X + P                                                    
      WR(EN) = X + P                                                    
      WI(NA) = ZZ                                                       
      WI(EN) = -ZZ                                                      
  330 EN = ENM2                                                         
      GO TO 60                                                          
C     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND               
C                VECTORS OF UPPER TRIANGULAR FORM **********            
  340 NORM = 0.0                                                        
      K = 1                                                             
C                                                                       
      DO 360 I = 1, N                                                   
C                                                                       
         DO 350 J = K, N                                                
  350    NORM = NORM + ABS(H(I,J))                                      
C                                                                       
         K = I                                                          
  360 CONTINUE                                                          
C                                                                       
      IF (NORM .EQ. 0.0) GO TO 1001                                     
C     ********** FOR EN=N STEP -1 UNTIL 1 DO -- **********              
      DO 800 NN = 1, N                                                  
         EN = N + 1 - NN                                                
         P = WR(EN)                                                     
         Q = WI(EN)                                                     
         NA = EN - 1                                                    
         IF (Q) 710, 600, 800                                           
C     ********** REAL VECTOR **********                                 
  600    M = EN                                                         
         H(EN,EN) = 1.0                                                 
         IF (NA .EQ. 0) GO TO 800                                       
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********            
         DO 700 II = 1, NA                                              
            I = EN - II                                                 
            W = H(I,I) - P                                              
            R = H(I,EN)                                                 
             IF (M .GT. NA) GO TO 620                                   
C                                                                       
            DO 610 J = M, NA                                            
  610       R = R + H(I,J) * H(J,EN)                                    
C                                                                       
  620       IF (WI(I) .GE. 0.0) GO TO 630                               
            ZZ = W                                                      
            S = R                                                       
            GO TO 700                                                   
  630       M = I                                                       
            IF (WI(I) .NE. 0.0) GO TO 640                               
            T = W                                                       
            IF (W .EQ. 0.0) T = MACHEP * NORM                           
            H(I,EN) = -R / T                                            
            GO TO 700                                                   
C     ********** SOLVE REAL EQUATIONS **********                        
  640       X = H(I,I+1)                                                
            Y = H(I+1,I)                                                
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)               
            T = (X * S - ZZ * R) / Q                                    
            H(I,EN) = T                                                 
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650                          
            H(I+1,EN) = (-R - W * T) / X                                
            GO TO 700                                                   
  650       H(I+1,EN) = (-S - Y * T) / ZZ                               
  700 CONTINUE                                                          
C     ********** END REAL VECTOR **********                             
         GO TO 800                                                      
C     ********** COMPLEX VECTOR **********                              
  710    M = NA                                                         
C     ********** LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT         
C                EIGENVECTOR MATRIX IS TRIANGULAR **********            
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720                
         H(NA,NA) = Q / H(EN,NA)                                        
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)                          
         GO TO 730                                                      
  720    Z3 = CMPLX(0.0,-H(NA,EN)) / CMPLX(H(NA,NA)-P,Q)                
         H(NA,NA) = T3(1)                                               
         H(NA,EN) = T3(2)                                               
  730    H(EN,NA) = 0.0                                                 
         H(EN,EN) = 1.0                                                 
         ENM2 = NA - 1                                                  
         IF (ENM2 .EQ. 0) GO TO 800                                     
C                                                                       
         DO 790 II = 1, ENM2                                            
            I = NA - II                                                 
            W = H(I,I) - P                                              
            RA = 0.0                                                    
            SA = H(I,EN)                                                
C                                                                       
            DO 760 J = M, NA                                            
               RA = RA + H(I,J) * H(J,NA)                               
               SA = SA + H(I,J) * H(J,EN)                               
  760       CONTINUE                                                    
C                                                                       
            IF (WI(I) .GE. 0.0) GO TO 770                               
            ZZ = W                                                      
            R = RA                                                      
            S = SA                                                      
            GO TO 790                                                   
  770       M = I                                                       
            IF (WI(I) .NE. 0.0) GO TO 780                               
            Z3 = CMPLX(-RA,-SA) / CMPLX(W,Q)                            
            H(I,NA) = T3(1)                                             
            H(I,EN) = T3(2)                                             
            GO TO 790                                                   
C     ********** SOLVE COMPLEX EQUATIONS **********                     
  780       X = H(I,I+1)                                                
            Y = H(I+1,I)                                                
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q      
            VI = (WR(I) - P) * 2.0 * Q                                  
            IF (VR .EQ. 0.0 .AND. VI .EQ. 0.0) VR = MACHEP * NORM       
     X       * (ABS(W) + ABS(Q) + ABS(X) + ABS(Y) + ABS(ZZ))            
            Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA) / CMPLX(VR,VI)    
            H(I,NA) = T3(1)                                             
            H(I,EN) = T3(2)                                             
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785                 
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X           
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X           
            GO TO 790                                                   
  785       Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN)) / CMPLX(ZZ,Q)         
            H(I+1,NA) = T3(1)                                           
            H(I+1,EN) = T3(2)                                           
  790    CONTINUE                                                       
C     ********** END COMPLEX VECTOR **********                          
  800 CONTINUE                                                          
C     ********** END BACK SUBSTITUTION.                                 
C                VECTORS OF ISOLATED ROOTS **********                   
      DO 840 I = 1, N                                                   
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840                     
C                                                                       
         DO 820 J = I, N                                                
  820    Z(I,J) = H(I,J)                                                
C                                                                       
  840 CONTINUE                                                          
C     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE              
C                VECTORS OF ORIGINAL FULL MATRIX.                       
C                FOR J=N STEP -1 UNTIL LOW DO -- **********             
      DO 880 JJ = LOW, N                                                
         J = N + LOW - JJ                                               
         M = MIN0(J,IGH)                                                
C                                                                       
         DO 880 I = LOW, IGH                                            
            ZZ = 0.0                                                    
C                                                                       
            DO 860 K = LOW, M                                           
  860       ZZ = ZZ + Z(I,K) * H(K,J)                                   
C                                                                       
            Z(I,J) = ZZ                                                 
  880 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000 IERR = EN                                                         
 1001 RETURN                                                            
C     ********** LAST CARD OF HQR2 **********                           
      END                                                               
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)                         
                                                                        
      INTEGER I,J,K,M,N,II,NM,IGH,LOW                                   
      REAL SCALE(N),Z(NM,M)                                             
      REAL S                                                            
COMMENT*XXXX;                                                           
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE BALBAK,              
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.              
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).   
C                                                                       
C     BALBAK   FORMS THE EIGENVECTORS OF A REAL GENERAL                 
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING            
C     BALANCED MATRIX DETERMINED BY  BALANC.                            
COMMENT*YYYY;                                                           
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC,                
C                                                                       
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS        
C          AND SCALING FACTORS USED BY  BALANC,                         
C                                                                       
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED,        
C                                                                       
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-          
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.       
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE                 
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.             
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
COMMENT*ZZZZ;                                                           
      IF (IGH .EQ. LOW) GO TO 120                                       
C                                                                       
      DO 110 I = LOW, IGH                                               
         S = SCALE(I)                                                   
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED            
C                IF THE FOREGOING STATEMENT IS REPLACED BY              
C                S=1.0/SCALE(I). **********                             
         DO 100 J = 1, M                                                
  100    Z(I,J) = Z(I,J) * S                                            
C                                                                       
  110 CONTINUE                                                          
C     ********- FOR I=LOW-1 STEP -1 UNTIL 1,                            
C               IGH+1 STEP 1 UNTIL N DO -- **********                   
  120 DO 140 II = 1, N                                                  
         I = II                                                         
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140                     
         IF (I .LT. LOW) I = LOW - II                                   
         K = SCALE(I)                                                   
         IF (K .EQ. I) GO TO 140                                        
C                                                                       
         DO 130 J = 1, M                                                
            S = Z(I,J)                                                  
            Z(I,J) = Z(K,J)                                             
            Z(K,J) = S                                                  
  130    CONTINUE                                                       
C                                                                       
  140 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C     ********** LAST CARD OF BALBAK **********                         
      END                                                               
      SUBROUTINE   NRMLZ(NM,N,WI,Z)                                     
      INTEGER      I,J,ISAVE,N,NM                                       
      REAL         S,S1,RECIPR,RECIPI                                   
      REAL         Z(NM,N), WI(N)                                       
COMMENT*XXXX                                                            
C     NORMALIZES A REAL EIGENVECTOR SO THAT ITS COMPONENT OF MAXIMUM    
C     MAGNITUDE IS 1.0 AND A COMPLEX EIGENVECTOR SO THAT ITS COMPONENT  
C     OF MAXIMUM MAGNITUDE IS 1.0 + 0.0*I                               
COMMENT*YYYY;                                                           
COMMENT*ZZZZ;                                                           
C     DO 500 J=1,N  QUESTA RIGA SOSTITUITA CO LE TRE SEGUENTI           
      J=0                                                               
  555 J=J+1                                                             
      IF (J.GT.N) GO TO 556                                             
         S = 0.0                                                        
         IF (WI(J) .NE. 0) GO TO 250                                    
C        REAL EIGENVECTOR                                               
         DO 100 I=1,N                                                   
            S1 = ABS(Z(I,J))                                            
            IF(S1 .LE. S) GO TO 100                                     
            S = S1                                                      
            ISAVE = I                                                   
 100     CONTINUE                                                       
         IF (S .NE. 0) S = 1./Z(ISAVE,J)                                
         DO 200 I=1,N                                                   
 200     Z(I,J) = Z(I,J)*S                                              
         GO TO 500                                                      
C        COMPLEX EIGENVECTOR STORED IN COLUMN J AND J+1                 
 250     DO 300 I=1,N                                                   
            S1 = Z(I,J)**2 + Z(I,J+1)**2                                
            IF (S1 .LE. S) GO TO 300                                    
            S = S1                                                      
            ISAVE = I                                                   
 300     CONTINUE                                                       
         IF (S .EQ. 0) GO TO 500                                        
         RECIPR = Z(ISAVE,J)/S                                          
         RECIPI =-Z(ISAVE,J+1)/S                                        
         DO 400 I=1,N                                                   
            S1       = Z(I,J)*RECIPR - Z(I,J+1)*RECIPI                  
            Z(I,J+1) = Z(I,J)*RECIPI + Z(I,J+1)*RECIPR                  
            Z(I,J)   = S1                                               
 400     CONTINUE                                                       
         Z(ISAVE,J) = 1.0                                               
         Z(ISAVE,J+1) = 0.0                                             
         J = J+1                                                        
C 500  CONTINUE  RIGA SOSTITUITA CON LE DUE SEGUENTI                    
 500  GO TO 555                                                         
 556  CONTINUE                                                          
      RETURN                                                            
      END                                                               
