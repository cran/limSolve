
C*********************************************************************
C LEAST DISTANCE SUBROUTINE
C*********************************************************************

      SUBROUTINE ldp(G,H,NUnknowns,NConstraints,NW,X,XNorm,W,xIndex,           &
     &               Mode,verbose, IsError)


      INTEGER           :: NUnknowns,NConstraints,NW
      DOUBLE PRECISION  :: G(NConstraints,NUnknowns) 
      DOUBLE PRECISION  :: H(NConstraints)

      DOUBLE PRECISION  :: X(NUnknowns)
      DOUBLE PRECISION  :: XNorm

      LOGICAL           :: IsError, verbose

      DOUBLE PRECISION  :: W(NW)
      INTEGER           :: xIndex(NConstraints)
      INTEGER           :: Mode
      INTEGER :: xLDPSucces,xLDPNoUnknownsOrEquations,                         &
     &       xLDPToomanyIterations,xLDPIncompatibleConstraints,                &
     &       xLDPUnsolvable


      PARAMETER (xLDPSucces = 1,xLDPNoUnknownsOrEquations     = 2,             &
     &       xLDPToomanyIterations= 3,xLDPIncompatibleConstraints   = 4,       &
     &       xLDPUnsolvable  = -1)


      CALL xLDP(G,NConstraints,NConstraints,NUnknowns,H,X,Xnorm,W,             &
     &          xINdex,Mode)

      IsError=.TRUE.
      if (mode == xLDPSucces) IsError = .FALSE.

      IF (verbose) THEN 
      SELECT CASE (mode)
      CASE (xLDPNoUnknownsOrEquations)
       CALL xMESSAGE ("No unknowns or equations")
      CASE (xLDPToomanyIterations)
       CALL xMESSAGE ("Too many iterations")
      CASE (xLDPIncompatibleConstraints)
       CALL xMESSAGE ("Incompatible constraints ")
      CASE (xLDPUnsolvable       )
       CALL xMESSAGE ("LDP problem unsolvable")
      END SELECT
      ENDIF

      RETURN
      END SUBROUTINE LDP



C****************************************************************************
C SOLVING LINEAR LEAST SQUARES WITH LINEAR INEQUALITIES and EQUALITIES
C****************************************************************************
  
C----------------------------------------------------------------------------
C  Solves the least squares model with equality and inequality constraints
C
C  Mathematically: minimise SUM(A*x-B)^2 
C  where:          E*x=F
C                  G*x>=h   x> 0
C  x(NumUnknowns), A(NumApproximate,NumUnknowns), B(NumApproximate)
C               E(NumEquations,NumUnknowns)     , F(NumEquations)
C               G(NumInequalities,NumUnknowns)  , H(NumInequalities)
C----------------------------------------------------------------------------

      SUBROUTINE lsei (NUnknowns,NEquations,NConstraints,NApproximate,         &
     &          A,B,E,F,G,H,X,mIP,mdW,mWS,IP,W,WS,lpr,ProgOpt,                 &
     &          verbose,IsError)
 
      IMPLICIT NONE

C The arrays and dimensions for the mass balance inversion

      INTEGER          :: NUnknowns,NEquations,NConstraints,NApproximate
      INTEGER          :: mIP,mdW,mWS,lpr
      LOGICAL          :: IsError, verbose
      DOUBLE PRECISION  ::  A (NApproximate,NUnknowns),                        &
     &                      B (NApproximate)          ,                        &
     &                      E (NEquations,NUnknowns)  ,                        &
     &                      F (NEquations)            ,                        &
     &                      G (NConstraints,NUnknowns),                        &
     &                      H (NConstraints)          ,                        &
     &                      X (NUnknowns)                     
C work arrays
      DOUBLE PRECISION :: W(MDW,NUnknowns+1),WS(mWS)
      INTEGER          :: IP(mip)

      INTEGER          :: I,J,K,MOde,ME, MA, MG,N
      DOUBLE PRECISION :: RNORME, RNORML,ProgOpt(lpr)

C---------------------------------------------------------------------
C a linear Least Squares solver with linear Equality and Inequality constraints: (LSEI)
C
C Find MIN ||A*X=B||
C
C and where
C    E*x=F, G*X>H
C    x> 0
C---------------------------------------------------------------------

      N  = NUnknowns
      ME = NEquations
      MA = NApproximate
      MG = NConstraints

C W is Working array with E,A,G   ,F,B,H 

      RNormE = 0.D0
      RNORML = 0.D0
      MODE   = 0

       DO I = 1, ME
         DO J = 1, N
           W(I,J) = E(I,J)
         ENDDO
         W(I,N+1) = F(I)
       ENDDO        

      K = ME
       DO I = 1, MA
         DO J = 1, N
           W(I+K,J) = A(I,J)
         ENDDO
         W(I+K,N+1) = B(I)
       ENDDO        
      K = ME+MA
      DO I = 1,NConstraints
         DO J = 1, N
           W(I+K,J) = G(I,J)
         ENDDO
         W(I+K,N+1) = H(I)
       ENDDO        
      K = K + NConstraints

c      ProgOpt(1) = 1.D0


    
C CALLING SOLVER!

        CALL xdLSEI(W,                                                   &     ! Linear equations array
     &              MDW,                                                 &     ! Dimensions
     &              ME,                                                  & 
     &              MA,                                                  &
     &              MG,                                                  &
     &              N,                                                   &
     &              ProgOpt,                                             &     ! Program options
     &              X,                                                   &     ! Return vector
     &              RNORME,                                              &     ! Residual of F-EX (should be 0)
     &              RNORML,                                              &     ! Residual of B-AX
     &              MODE,                                                &     ! 0 = success; 1: equality contradictory,2:inequality contradict
     &              WS,                                                  &
     &              IP)


      IF (verbose) THEN
       SELECT CASE (Mode) 
       CASE(1)
           CALL XMESSAGE ("LSEI error: equalities contradictory")

       CASE(2)
           CALL XMESSAGE ("LSEI error: inequalities contradictory")

       CASE(3)
           CALL XMESSAGE                                                  &
     &    ("LSEI error: equalities + inequalities contradictory")

       CASE(4)
           CALL XMESSAGE("LSEI error: wrong input")       
       END SELECT
      ENDIF
      IsError = .FALSE.
      IF (mode > 0) Iserror = .TRUE.
      RETURN

      END SUBROUTINE LSEI

C************************************************************************
C TRIDIAGONAL MATRIX SOLVERS
C************************************************************************
 

      SUBROUTINE tridia(Nmx,au,bu,cu,du,value)

      IMPLICIT NONE

C-------------------------------------------------------------------------
C Solves a tridiagonal system of linear equations A*X=DU by backsubstitution
C where the non-zero diagonal elements are stored in au (below diagonal), 
C bu (on the diagonal), cu (above diagonal)
C and the equations are from first to last
C-------------------------------------------------------------------------

      INTEGER            Nmx                   
      DOUBLE PRECISION   au(Nmx)               
      DOUBLE PRECISION   bu(Nmx)               
      DOUBLE PRECISION   cu(Nmx)               
      DOUBLE PRECISION   du(Nmx)               

      DOUBLE PRECISION   value(Nmx)    
C
      DOUBLE PRECISION   ru(Nmx),qu(Nmx)
      INTEGER            I
C-------------------------------------------------------------------------
C Backsubstitution
      ru(Nmx)=au(Nmx)/bu(Nmx)
      qu(Nmx)=du(Nmx)/bu(Nmx)

      DO 10 I=Nmx-1,2,-1
         ru(I)=au(I)/(bu(I)-cu(I)*ru(I+1))
         qu(I)=(du(I)-cu(I)*qu(I+1))/(bu(I)-cu(I)*ru(I+1))
10    CONTINUE
 
      qu(1)=(du(1)-cu(1)*qu(2))/(bu(1)-cu(1)*ru(2))

C Forward substitution
      value(1)=qu(1)

      DO 20 I=2,Nmx
         value(I)=qu(I)-ru(I)*value(I-1)
20    CONTINUE

      RETURN
      END SUBROUTINE TRIDIA

C************************************************************************
C BANDED MATRIX SOLVER
C************************************************************************


      SUBROUTINE banded(abd,beta,NumABD,NumSvar,BandDown,BandUp,                 &
     &                  indx,info,IsError)
      INTEGER           NumABD,NumSvar,BandDown,BandUp
      INTEGER           indx(NumSvar), info
      DOUBLE PRECISION  ABD(NumABD,NumSVar),beta(NumSvar)
      LOGICAL           isError

        IsError = .FALSE.

        CALL dgbfa(abd,NumABD,NumSvar,BandDown,BandUp,indx,info)

        IF (Info .LT. 0) THEN
C xdgbsl will divide by 0 (singular matrix)
          IsError = .TRUE.
          RETURN
        ENDIF 
        CALL dgbsl(abd,NumABD,NumSvar,BandDown,BandUp,indx,beta,info)
      RETURN
      END SUBROUTINE BANDED   




C######################################################################
C Programs and subroutines relating to reading/writing from and to ascii files
C
C######################################################################


C                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
C                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
C                !            ERROR HANDLING          !
C                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
C                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!


      SUBROUTINE XMESSAGE (String)

      CHARACTER (LEN=*)   :: String

C Check whether it is safe to write
        CALL rwarn(String)

      END SUBROUTINE XMESSAGE



C######################################################################
C 
C Routines for solving quadratic and linear programming problem
C
C SUBROUTINES FROM LINPACK LIBRARY
C
C MAIN ROUTINES ARE xLDP
C xLSEI and  xdbocls
C 
C######################################################################



C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C    Quadratic programming solver    C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C


C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 15, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C  made compatible with fortran 95 by karline Soetaert
C  added x-prefix
C  captured writing to screen -> XMESSAGE

C*************************************************************************C
C LEAST DISTANCE SUBROUTINE
C*************************************************************************C


      SUBROUTINE xLDP (G,MDG,M,N,H,X,XNORM,W,xINDEX,MODE)     
C
C  Algorithm LDP: LEAST DISTANCE PROGRAMMING
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1974 MAR 1, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER  :: xLDPSucces ,xLDPNoUnknownsOrEquations 
      INTEGER  :: xLDPToomanyIterations,xLDPIncompatibleConstraints
      INTEGER  :: xLDPUnsolvable 
      PARAMETER(xLDPSucces= 1,xLDPNoUnknownsOrEquations     = 2,                 &
     &                      xLDPToomanyIterations         = 3,                   &
     &                      xLDPIncompatibleConstraints   = 4,                   &
     &                      xLDPUnsolvable                = -1)
C Number of unknowns
      INTEGER  :: M, MDG,N         
C Succes or failure
      INTEGER  :: MODE             
      INTEGER  :: xINDEX(*)  
C Constraints G*X>H ; W=workarray; xnorm=residual norm
      DOUBLE PRECISION :: G(MDG,*), H(*), X(*)
      DOUBLE PRECISION :: W(*)                
      DOUBLE PRECISION :: XNORM               

      INTEGER          :: I, IW, IWDUAL, IY, IZ, J, JF, NP1
      DOUBLE PRECISION :: xDIFF, FAC, RNORM
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------

      MODE=xLDPSucces

      IF (N.LE.0) THEN
        MODE=xLDPNoUnknownsOrEquations
        RETURN
      ENDIF

      DO J=1,N   
        X(J)=ZERO     
      ENDDO
      XNORM=ZERO

      IF (M.LE.0) THEN
        MODE=xLDPNoUnknownsOrEquations           
        RETURN
      ENDIF
C   
C     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.   
C   
C      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.    
C     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.   
C     COPY H**T INTO ROW N+1 OF E.  
C   
      IW=0  
      DO J=1,M   
          DO I=1,N   
               IW=IW+1   
               W(IW)=G(J,I)  
          ENDDO
          IW=IW+1   
          W(IW)=H(J) 
      ENDDO   
      JF=IW+1   
C                                STORE N ZEROS FOLLOWED BY A ONE INTO F.
      DO I=1,N   
          IW=IW+1   
          W(IW)=ZERO    
      ENDDO
      W(IW+1)=ONE   
C   
      NP1=N+1   
      IZ=IW+2   
      IY=IZ+NP1 
      IWDUAL=IY+M   
C   
      CALL xNNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),         &
     &  xINDEX,MODE)  
C                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
      IF (MODE.NE.xLDPSucces) RETURN 

C      IF (RNORM) 130,130,50     KARLINE changed
C  50  
      IF (RNORM .LE. 0) THEN
C karline:remove ????
        MODE=xLDPUnsolvable
        RETURN
      ENDIF

      FAC=ONE   
      IW=IY-1   
      DO I=1,M   
         IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
         FAC=FAC-H(I)*W(IW)
      ENDDO
C   
C      IF (DIFF(ONE+FAC,ONE)) 130,130,70 
C   70 FAC=ONE/FAC   

      IF (xDIFF(ONE+FAC,ONE) .LE. 0) THEN
        MODE=xLDPIncompatibleConstraints
        RETURN
      ENDIF
      FAC=ONE/FAC   

      DO J=1,N   
         IW=IY-1   
         DO I=1,M   
            IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
            X(J)=X(J)+G(I,J)*W(IW)
         ENDDO 
         X(J)=X(J)*FAC 
      ENDDO

      DO J=1,N  
         XNORM=XNORM+X(J)**2   
      ENDDO

      XNORM=sqrt(XNORM) 
      RETURN
      END SUBROUTINE xLDP



C*************************************************************************C

C     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C   
C  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
C   
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 15, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM   
C   
C                      A * X = B  SUBJECT TO X .GE. 0   
C     ------------------------------------------------------------------
C                     Subroutine Arguments
C
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
C                     MATRIX, A.           ON EXIT A() CONTAINS 
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
C                     THIS SUBROUTINE.  
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
C             CONTAIN THE SOLUTION VECTOR.  
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
C             RESIDUAL VECTOR.  
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
C     ZZ()     AN M-ARRAY OF WORKING SPACE.     
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
C                 P AND Z AS FOLLOWS..  
C   
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N   
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
C             MEANINGS. 
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. 
C   
C     ------------------------------------------------------------------
      SUBROUTINE xnnls (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) 
C     ------------------------------------------------------------------
      integer I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ, L
      integer M, MDA, MODE,N, NPP1, NSETP, RTNKEY
C     integer INDEX(N)  
C     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M)   
      integer INDEX(*)  
      double precision A(MDA,*), B(*), W(*), X(*), ZZ(*)   
      double precision ALPHA, ASAVE, CC, xDIFF, DUMMY, FACTOR, RNORM
      double precision SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision ZERO, ZTEST
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=3*N 
C   
C                    INITIALIZE THE ARRAYS INDEX() AND X(). 
C   
      DO I=1,N   
         X(I)=ZERO     
         INDEX(I)=I    
      ENDDO
C   
      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
C                             &*****  MAIN LOOP BEGINS HERE  ******     
   30 CONTINUE  
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
C   
      IF (IZ1 .GT.IZ2.OR.NSETP.GE.M) GO TO 350   
C   
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C   
      DO IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         SM=ZERO   
         DO L=NPP1,M
           SM=SM+A(L,J)*B(L)   
         ENDDO  
         W(J)=SM   
      ENDDO
C                                   FIND LARGEST POSITIVE W(J). 
   60 continue
      WMAX=ZERO 
      DO IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)     
            IZMAX=IZ  
         endif
      ENDDO
C   
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C   
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX  
      J=INDEX(IZ)   
C   
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
C     NEAR LINEAR DEPENDENCE.   
C   
      ASAVE=A(NPP1,J)   
      CALL xH12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO L=1,NSETP   
            UNORM=UNORM+A(L,J)**2     
          ENDDO
      endif
      UNORM=sqrt(UNORM) 
      IF (xDIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
C   
C        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
C        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
C   
        DO L=1,M  
          ZZ(L)=B(L)    
        ENDDO
        CALL xH12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
        ZTEST=ZZ(NPP1)/A(NPP1,J)  
C   
C                                     SEE IF ZTEST IS POSITIVE  
C   
         IF (ZTEST .gt. ZERO) go to 140
      endif
C   
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.     
C   
      A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
C   
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
C     COL J,  SET W(J)=0.   
C   
  140 continue
      DO L=1,M  
        B(L)=ZZ(L)    
      ENDDO
C   
      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   
C   
      IF (IZ1 .le. IZ2) then
         DO JZ=IZ1,IZ2 
            JJ=INDEX(JZ)  
            CALL xH12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
         ENDDO
      endif
C   
      IF (NSETP .ne. M) then
         DO L=NPP1,M   
           A(L,J)=ZERO   
         ENDDO
      endif
C   
      W(J)=ZERO 
C                                SOLVE THE TRIANGULAR SYSTEM.   
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      RTNKEY = 1
      GO TO 400 
  200 CONTINUE  
C   
C                       &*****  SECONDARY LOOP BEGINS HERE ******   
C   
C                          ITERATION COUNTER.   
C 
  210 continue  
      ITER=ITER+1   
      IF (ITER .gt. ITMAX) then
         MODE=3
C         write (*,'(/a)') ' NNLS quitting on iteration count.'
      CALL XMESSAGE ('error in LDP - NNLS quitting on iteration count.')
         GO TO 350 
      endif
C   
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
C                                  IF NOT COMPUTE ALPHA.    
C   
      ALPHA=TWO 
      DO IP=1,NSETP 
         L=INDEX(IP)   
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))     
            IF (ALPHA .gt. T) then
               ALPHA=T   
               JJ=IP 
            endif
         endif
      ENDDO 
C   
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
C   
      IF (ALPHA.EQ.TWO) GO TO 330   
C   
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
C   
      DO IP=1,NSETP 
         L=INDEX(IP)   
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
      ENDDO
C   
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
C        FROM SET P TO SET Z.   
C   
      I=INDEX(JJ)   
  260 continue
      X(I)=ZERO 
C   
      IF (JJ .ne. NSETP) then
         JJ=JJ+1   
         DO 280 J=JJ,NSETP 
            II=INDEX(J)   
            INDEX(J-1)=II 
            CALL xG1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
            A(J,II)=ZERO  
            DO L=1,N  
               IF (L.NE.II) then
C
C                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))  
C
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
            ENDDO   
C
C                 Apply procedure G2 (CC,SS,B(J-1),B(J))   
C
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)    
            B(J)   =-SS*TEMP + CC*B(J)    
  280    continue
      endif
C
      NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
C   
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
C        AND MOVED FROM SET P TO SET Z. 
C   
      DO JJ=1,NSETP 
         I=INDEX(JJ)   
         IF (X(I) .le. ZERO) go to 260
      ENDDO 
C   
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C   
      DO I=1,M  
        ZZ(I)=B(I)    
      ENDDO
      RTNKEY = 2
      GO TO 400 
  320 CONTINUE  
      GO TO 210 
C                      &*****  END OF SECONDARY LOOP  ******
C   
  330 continue
      DO IP=1,NSETP 
        I=INDEX(IP)   
        X(I)=ZZ(IP)   
      ENDDO
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
      GO TO 30  
C   
C                        &*****  END OF MAIN LOOP  ******   
C   
C                        COME TO HERE FOR TERMINATION.  
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
C 
  350 continue  
      SM=ZERO   
      IF (NPP1 .le. M) then
         DO I=NPP1,M   
           SM=SM+B(I)**2 
         ENDDO
      else
         DO J=1,N  
           W(J)=ZERO     
         ENDDO
      endif
      RNORM=sqrt(SM)    
      RETURN
C   
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
C   
  400 continue
      DO L=1,NSETP  
         IP=NSETP+1-L  
         IF (L .ne. 1) then
            DO II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
            ENDDO
         endif
         JJ=INDEX(IP)  
         ZZ(IP)=ZZ(IP)/A(IP,JJ)    
      ENDDO
C
      go to (200, 320), RTNKEY
      END SUBROUTINE xNNLS


     

C*************************************************************************C

C     SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)    
C
C  QR ALGORITHM FOR SINGULAR VALUES OF A BIDIAGONAL MATRIX.
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 12, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------   
C     THE BIDIAGONAL MATRIX 
C   
C                       (Q1,E2,0...    )
C                       (   Q2,E3,0... )
C                D=     (       .      )
C                       (         .   0)
C                       (           .EN)
C                       (          0,QN)
C   
C                 IS PRE AND POST MULTIPLIED BY 
C                 ELEMENTARY ROTATION MATRICES  
C                 RI AND PI SO THAT 
C   
C                 RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SN) 
C   
C                 TO WITHIN WORKING ACCURACY.   
C   
C  1. EI AND QI OCCUPY E(I) AND Q(I) AS INPUT.  
C   
C  2. RM...R1*C REPLACES 'C' IN STORAGE AS OUTPUT.  
C   
C  3. V*P1**(T)...PM**(T) REPLACES 'V' IN STORAGE AS OUTPUT.
C   
C  4. SI OCCUPIES Q(I) AS OUTPUT.   
C   
C  5. THE SI'S ARE NONINCREASING AND NONNEGATIVE.   
C   
C     THIS CODE IS BASED ON THE PAPER AND 'ALGOL' CODE..    
C REF..     
C  1. REINSCH,C.H. AND GOLUB,G.H. 'SINGULAR VALUE DECOMPOSITION 
C     AND LEAST SQUARES SOLUTIONS' (NUMER. MATH.), VOL. 14,(1970).  
C   
C     ------------------------------------------------------------------   
      SUBROUTINE xQRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)    
C     ------------------------------------------------------------------   
      integer MDC, MDV, NCC, NN, NRV
C     double precision C(MDC,NCC), E(NN), Q(NN),V(MDV,NN)
      double precision C(MDC,*  ), E(* ), Q(* ),V(MDV,* )
      integer I, II, IPASS, J, K, KK, L, LL, LP1, N, N10, NQRS
      double precision CS, xDIFF, DNORM, F, G, H, SMALL
      double precision ONE, SN, T, TEMP, TWO, X, Y, Z, ZERO
      
      logical WNTV ,HAVERS,FAIL     
      parameter(ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------   
      N=NN  
      IPASS=1   
      IF (N.LE.0) RETURN
      N10=10*N  
      WNTV=NRV.GT.0     
      HAVERS=NCC.GT.0   
      FAIL=.FALSE.  
      NQRS=0
      E(1)=ZERO 
      DNORM=ZERO
      DO J=1,N  
        DNORM=max(abs(Q(J))+abs(E(J)),DNORM)   
      ENDDO
      DO 200 KK=1,N
           K=N+1-KK     
C   
C     TEST FOR SPLITTING OR RANK DEFICIENCIES.. 
C         FIRST MAKE TEST FOR LAST DIAGONAL TERM, Q(K), BEING SMALL.    
   20       IF(K.EQ.1) GO TO 50     
            IF(xDIFF(DNORM+Q(K),DNORM) .ne. ZERO) go to 50
C   
C     SINCE Q(K) IS SMALL WE WILL MAKE A SPECIAL PASS TO    
C     TRANSFORM E(K) TO ZERO.   
C   
           CS=ZERO  
           SN=-ONE  
                DO 40 II=2,K
                I=K+1-II
                F=-SN*E(I+1)
                E(I+1)=CS*E(I+1)    
                CALL xG1 (Q(I),F,CS,SN,Q(I))     
C         TRANSFORMATION CONSTRUCTED TO ZERO POSITION (I,K).
C   
                IF (.NOT.WNTV) GO TO 40 
                     DO 30 J=1,NRV  
C
C                          Apply procedure G2 (CS,SN,V(J,I),V(J,K))  
C
                        TEMP = V(J,I)
                        V(J,I) = CS*TEMP + SN*V(J,K)
                        V(J,K) =-SN*TEMP + CS*V(J,K)
   30                continue
C              ACCUMULATE RT. TRANSFORMATIONS IN V. 
C   
   40           CONTINUE
C   
C         THE MATRIX IS NOW BIDIAGONAL, AND OF LOWER ORDER  
C         SINCE E(K) .EQ. ZERO..    
C   
   50           DO LL=1,K
                  L=K+1-LL
                  IF(xDIFF(DNORM+E(L),DNORM) .eq. ZERO) go to 100
                  IF(xDIFF(DNORM+Q(L-1),DNORM) .eq. ZERO) go to 70
                ENDDO
C     THIS LOOP CAN'T COMPLETE SINCE E(1) = ZERO.   
C   
           GO TO 100    
C   
C         CANCELLATION OF E(L), L.GT.1. 
   70      CS=ZERO  
           SN=-ONE  
              DO I=L,K 
                F=-SN*E(I)  
                E(I)=CS*E(I)
                IF(xDIFF(DNORM+F,DNORM) .eq. ZERO) go to 100
                CALL xG1 (Q(I),F,CS,SN,Q(I))     
                IF (HAVERS) then
                     DO J=1,NCC  
C
C                          Apply procedure G2 ( CS, SN, C(I,J), C(L-1,J)
C
                        TEMP = C(I,J)
                        C(I,J)   = CS*TEMP + SN*C(L-1,J)
                        C(L-1,J) =-SN*TEMP + CS*C(L-1,J)
                     ENDDO
                endif
              ENDDO  
C   
C         TEST FOR CONVERGENCE..    
  100      Z=Q(K)   
           IF (L.EQ.K) GO TO 170    
C   
C         SHIFT FROM BOTTOM 2 BY 2 MINOR OF B**(T)*B.   
           X=Q(L)   
           Y=Q(K-1)     
           G=E(K-1)     
           H=E(K)   
           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
           G=sqrt(ONE+F**2) 
           IF (F .ge. ZERO) then
              T=F+G
           else
              T=F-G
           endif
           F=((X-Z)*(X+Z)+H*(Y/T-H))/X  
C   
C         NEXT QR SWEEP..   
           CS=ONE   
           SN=ONE   
           LP1=L+1  
                DO 160 I=LP1,K  
                G=E(I)  
                Y=Q(I)  
                H=SN*G  
                G=CS*G  
                CALL xG1 (F,H,CS,SN,E(I-1))  
                F=X*CS+G*SN 
                G=-X*SN+G*CS
                H=Y*SN  
                Y=Y*CS  
                IF (WNTV) then
C   
C              ACCUMULATE ROTATIONS (FROM THE RIGHT) IN 'V' 
C
                     DO J=1,NRV 
C                    
C                          Apply procedure G2 (CS,SN,V(J,I-1),V(J,I))
C
                        TEMP = V(J,I-1)
                        V(J,I-1) = CS*TEMP + SN*V(J,I)
                        V(J,I)   =-SN*TEMP + CS*V(J,I)
                     ENDDO
                endif
                CALL xG1 (F,H,CS,SN,Q(I-1))  
                F=CS*G+SN*Y 
                X=-SN*G+CS*Y
                IF (HAVERS) then
                     DO J=1,NCC 
C
C                          Apply procedure G2 (CS,SN,C(I-1,J),C(I,J))
C
                        TEMP = C(I-1,J)
                        C(I-1,J) = CS*TEMP + SN*C(I,J)
                        C(I,J)   =-SN*TEMP + CS*C(I,J)
                     ENDDO
                endif
C
C              APPLY ROTATIONS FROM THE LEFT TO 
C              RIGHT HAND SIDES IN 'C'..
C   
  160           CONTINUE
           E(L)=ZERO    
           E(K)=F   
           Q(K)=X   
           NQRS=NQRS+1  
           IF (NQRS.LE.N10) GO TO 20
C          RETURN TO 'TEST FOR SPLITTING'.  
C 
           SMALL=ABS(E(K))
           I=K     
C          IF FAILURE TO CONVERGE SET SMALLEST MAGNITUDE
C          TERM IN OFF-DIAGONAL TO ZERO.  CONTINUE ON.
C      ..           
                DO 165 J=L,K 
                TEMP=ABS(E(J))
                IF(TEMP .EQ. ZERO) GO TO 165
                IF(TEMP .LT. SMALL) THEN
                     SMALL=TEMP
                     I=J
                end if  
  165           CONTINUE
           E(I)=ZERO
           NQRS=0
           FAIL=.TRUE.  
           GO TO 20
C     ..    
C     CUTOFF FOR CONVERGENCE FAILURE. 'NQRS' WILL BE 2*N USUALLY.   
  170      IF (Z.GE.ZERO) GO TO 190 
           Q(K)=-Z  
           IF (WNTV) then
           DO J=1,NRV  
             V(J,K)=-V(J,K)  
           ENDDO 
           endif
  190      CONTINUE     
C         CONVERGENCE. Q(K) IS MADE NONNEGATIVE..   
C   
  200      CONTINUE     
      IF (N.EQ.1) RETURN
           DO 210 I=2,N 
           IF (Q(I).GT.Q(I-1)) GO TO 220
  210      CONTINUE     
      IF (FAIL) IPASS=2 
      RETURN
C     ..    
C     EVERY SINGULAR VALUE IS IN ORDER..
  220      DO 270 I=2,N 
           T=Q(I-1)     
           K=I-1
                DO 230 J=I,N
                IF (T.GE.Q(J)) GO TO 230
                T=Q(J)  
                K=J     
  230           CONTINUE
           IF (K.EQ.I-1) GO TO 270  
           Q(K)=Q(I-1)  
           Q(I-1)=T     
           IF (HAVERS) then
           DO J=1,NCC  
             T=C(I-1,J)  
             C(I-1,J)=C(K,J)     
             C(K,J)=T
           ENDDO
           endif

           IF (WNTV) then
             DO J=1,NRV  
               T=V(J,I-1)  
               V(J,I-1)=V(J,K)     
               V(J,K)=T
             ENDDO
           endif
  270      CONTINUE     
C         END OF ORDERING ALGORITHM.
C   
      IF (FAIL) IPASS=2 
      RETURN
      END SUBROUTINE xQRBD


C*************************************************************************C

C     SUBROUTINE SVDRS (A, MDA, M1, N1, B, MDB, NB, S, WORK) 
C
C  SINGULAR VALUE DECOMPOSITION ALSO TREATING RIGHT SIDE VECTOR.
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1974 SEP 25, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------   
C  This 1995 version differs from the original 1974 version by adding
C  the argument WORK().
C  WORK() provides 2*N1 locations of work space.  Originally S() was
C  required to have 3*N1 elements, of which the last 2*N1 were used for
C  work space.  Now S() only needs N1 elements.
C     ------------------------------------------------------------------
C  This subroutine computes the singular value decomposition of the
C  given M1 x N1 matrix, A, and optionally applys the transformations
C  from the left to the NB column vectors of the M1 x NB matrix B.
C  Either M1 .ge. N1  or  M1 .lt. N1 is permitted.
C
C       The singular value decomposition of A is of the form
C
C                  A  =  U * S * V**t
C
C  where U is M1 x M1 orthogonal, S is M1 x N1 diagonal with the
C  diagonal terms nonnegative and ordered from large to small, and
C  V is N1 x N1 orthogonal.  Note that these matrices also satisfy
C
C                  S = (U**t) * A * V
C
C       The matrix V is returned in the leading N1 rows and
C  columns of the array A(,).
C
C       The singular values, i.e. the diagonal terms of the matrix S,
C  are returned in the array S().  If M1 .lt. N1, positions M1+1
C  through N1 of S() will be set to zero.
C
C       The product matrix  G = U**t * B replaces the given matrix B
C  in the array B(,).
C
C       If the user wishes to obtain a minimum length least squares
C  solution of the linear system
C
C                        A * X ~=~ B
C
C  the solution X can be constructed, following use of this subroutine,
C  by computing the sum for i = 1, ..., R of the outer products
C
C          (Col i of V) * (1/S(i)) * (Row i of G)
C
C  Here R denotes the pseudorank of A which the user may choose
C  in the range 0 through Min(M1, N1) based on the sizes of the
C  singular values.
C     ------------------------------------------------------------------
C                    Subroutine Arguments
C
C  A(,)     (In/Out)  On input contains the M1 x N1 matrix A.
C           On output contains the N1 x N1 matrix V.
C
C  LDA      (In)  First dimensioning parameter for A(,).
C           Require LDA .ge. Max(M1, N1).
C
C  M1       (In)  No. of rows of matrices A, B, and G.
C           Require M1 > 0.
C
C  N1       (In)  No. of cols of matrix A, No. of rows and cols of
C           matrix V.  Permit M1 .ge. N1  or  M1 .lt. N1.
C           Require N1 > 0.
C
C  B(,)     (In/Out)  If NB .gt. 0 this array must contain an
C           M1 x NB matrix on input and will contain the
C           M1 x NB product matrix, G = (U**t) * B on output.
C
C  LDB      (In)  First dimensioning parameter for B(,).
C           Require LDB .ge. M1.
C
C  NB       (In)  No. of cols in the matrices B and G.
C           Require NB .ge. 0.
C
C  S()      (Out)  Must be dimensioned at least N1.  On return will
C           contain the singular values of A, with the ordering
C                S(1) .ge. S(2) .ge. ... .ge. S(N1) .ge. 0.
C           If M1 .lt. N1 the singular values indexed from M1+1
C           through N1 will be zero.
C           If the given integer arguments are not consistent, this
C           subroutine will return immediately, setting S(1) = -1.0.
C
C  WORK()  (Scratch)  Work space of total size at least 2*N1.
C           Locations 1 thru N1 will hold the off-diagonal terms of
C           the bidiagonal matrix for subroutine QRBD.  Locations N1+1
C           thru 2*N1 will save info from one call to the next of
C           H12.
C     ------------------------------------------------------------------
C  This code gives special treatment to rows and columns that are
C  entirely zero.  This causes certain zero sing. vals. to appear as
C  exact zeros rather than as about MACHEPS times the largest sing. val.
C  It similarly cleans up the associated columns of U and V.  
C
C  METHOD..  
C   1. EXCHANGE COLS OF A TO PACK NONZERO COLS TO THE LEFT.  
C      SET N = NO. OF NONZERO COLS.  
C      USE LOCATIONS A(1,N1),A(1,N1-1),...,A(1,N+1) TO RECORD THE    
C      COL PERMUTATIONS. 
C   2. EXCHANGE ROWS OF A TO PACK NONZERO ROWS TO THE TOP.   
C      QUIT PACKING IF FIND N NONZERO ROWS.  MAKE SAME ROW EXCHANGES 
C      IN B.  SET M SO THAT ALL NONZERO ROWS OF THE PERMUTED A   
C      ARE IN FIRST M ROWS.  IF M .LE. N THEN ALL M ROWS ARE 
C      NONZERO.  IF M .GT. N THEN      THE FIRST N ROWS ARE KNOWN    
C      TO BE NONZERO,AND ROWS N+1 THRU M MAY BE ZERO OR NONZERO.     
C   3. APPLY ORIGINAL ALGORITHM TO THE M BY N PROBLEM.   
C   4. MOVE PERMUTATION RECORD FROM A(,) TO S(I),I=N+1,...,N1.   
C   5. BUILD V UP FROM  N BY N  TO  N1 BY N1  BY PLACING ONES ON     
C      THE DIAGONAL AND ZEROS ELSEWHERE.  THIS IS ONLY PARTLY DONE   
C      EXPLICITLY.  IT IS COMPLETED DURING STEP 6.   
C   6. EXCHANGE ROWS OF V TO COMPENSATE FOR COL EXCHANGES OF STEP 2. 
C   7. PLACE ZEROS IN  S(I),I=N+1,N1  TO REPRESENT ZERO SING VALS.   
C     ------------------------------------------------------------------
      SUBROUTINE xSVDRS (A, MDA, M1, N1, B, MDB, NB, S, WORK) 
      integer I, IPASS, J, K, L, M, MDA, MDB, M1
      integer N, NB, N1, NP1, NS, NSP1
C     double precision A(MDA,N1),B(MDB,NB), S(N1)
      double precision A(MDA, *),B(MDB, *), S( *)
      double precision ONE, T, WORK(N1,2), ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                             BEGIN.. SPECIAL FOR ZERO ROWS AND COLS.   
C   
C                             PACK THE NONZERO COLS TO THE LEFT 
C   
      N=N1  
      IF (N.LE.0.OR.M1.LE.0) RETURN 
      J=N   

   10 CONTINUE  

      DO I=1,M1  
         IF (A(I,J) .ne. ZERO) go to 50
      ENDDO
C   
C           COL J  IS ZERO. EXCHANGE IT WITH COL N.   
C   
      IF (J .ne. N) THEN
         DO I=1,M1  
           A(I,J)=A(I,N) 
         ENDDO
      ENDIF

      A(1,N)=J  
      N=N-1 

   50    CONTINUE  
      J=J-1 
      IF (J.GE.1) GO TO 10  

C                             IF N=0 THEN A IS ENTIRELY ZERO AND SVD    
C                             COMPUTATION CAN BE SKIPPED    
      NS=0  
      IF (N.EQ.0) GO TO 240 
C                             PACK NONZERO ROWS TO THE TOP  
C                             QUIT PACKING IF FIND N NONZERO ROWS   
      I=1   
      M=M1  
   60 IF (I.GT.N.OR.I.GE.M) GO TO 150   

C KARLINE - CHANGED
C   IF (A(I,I)) 90,70,90  
C   70     DO 80 J=1,N   
C          IF (A(I,J)) 90,80,90  
C   80     CONTINUE  
C   90 I=I+1 
C KARLINE - CHANGED - tricky


      IF (A(I,I) .NE. 0) GOTO 90
      DO J=1,N   
          IF (A(I,J) .NE. 0) GOTO 90
      ENDDO

      GO TO 100 
   90 I=I+1 
      GO TO 60  
C                             ROW I IS ZERO     
C                             EXCHANGE ROWS I AND M 
  100 IF(NB.LE.0) GO TO 115 
      DO J=1,NB 
         T=B(I,J)  
         B(I,J)=B(M,J) 
         B(M,J)=T  
      ENDDO
  115 DO J=1,N  
         A(I,J)=A(M,J) 
      ENDDO
      IF (M.GT.N) GO TO 140 
      DO J=1,N  
         A(M,J)=ZERO   
      ENDDO
  140 CONTINUE  
C                             EXCHANGE IS FINISHED  
      M=M-1 
      GO TO 60  
C   
  150 CONTINUE  
C                             END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
C                             BEGIN.. SVD ALGORITHM 
C     METHOD..  
C     (1)     REDUCE THE MATRIX TO UPPER BIDIAGONAL FORM WITH   
C     HOUSEHOLDER TRANSFORMATIONS.  
C          H(N)...H(1)AQ(1)...Q(N-2) = (D**T,0)**T  
C     WHERE D IS UPPER BIDIAGONAL.  
C   
C     (2)     APPLY H(N)...H(1) TO B.  HERE H(N)...H(1)*B REPLACES B    
C     IN STORAGE.   
C   
C     (3)     THE MATRIX PRODUCT W= Q(1)...Q(N-2) OVERWRITES THE FIRST  
C     N ROWS OF A IN STORAGE.   
C   
C     (4)     AN SVD FOR D IS COMPUTED.  HERE K ROTATIONS RI AND PI ARE 
C     COMPUTED SO THAT  
C          RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SM)    
C     TO WORKING ACCURACY.  THE SI ARE NONNEGATIVE AND NONINCREASING.   
C     HERE RK...R1*B OVERWRITES B IN STORAGE WHILE  
C     A*P1**(T)...PK**(T)  OVERWRITES A IN STORAGE. 
C   
C     (5)     IT FOLLOWS THAT,WITH THE PROPER DEFINITIONS,  
C     U**(T)*B OVERWRITES B, WHILE V OVERWRITES THE FIRST N ROW AND     
C     COLUMNS OF A.     
C   
      L=min(M,N)   
C             THE FOLLOWING LOOP REDUCES A TO UPPER BIDIAGONAL AND  
C             ALSO APPLIES THE PREMULTIPLYING TRANSFORMATIONS TO B.     
C   CKarline: change from M to M-1 ???
      DO J=1,L  
       IF (J.GE.M) GO TO 160       
       CALL xH12 (1,J,J+1,M,A(1,J),1,T,A(1,J+1),1,MDA,N-J)
       CALL xH12 (2,J,J+1,M,A(1,J),1,T,B,1,MDB,NB)
C  160     IF (J.GE.N-1) GO TO 170   
  160  IF (J.GE.N-1) CYCLE
       CALL xH12 (1,J+1,J+2,N,A(J,1),MDA,work(J,2),A(J+1,1),MDA,1,M-J)
      ENDDO
C   
C     COPY THE BIDIAGONAL MATRIX INTO S() and WORK() FOR QRBD.   
C 1986 Jan 8. C. L. Lawson. Changed N to L in following 2 statements.
      IF (L.EQ.1) GO TO 190 
      DO J=2,L  
          S(J)=A(J,J) 
          WORK(J,1)=A(J-1,J)   
      ENDDO  
  190 S(1)=A(1,1)     
C   
      NS=N  
      IF (M.GE.N) GO TO 200 
      NS=M+1
      S(NS)=ZERO  
      WORK(NS,1)=A(M,M+1)  
  200 CONTINUE  
C   
C     CONSTRUCT THE EXPLICIT N BY N PRODUCT MATRIX, W=Q1*Q2*...*QL*I    
C     IN THE ARRAY A(). 
C   
       DO K=1,N  
         I=N+1-K   
         IF (I .GT. min(M,N-2)) GO TO 210     
        CALL xH12 (2,I+1,I+2,N,A(I,1),MDA,WORK(I,2),A(1,I+1),1,MDA,N-I)
  210    DO J=1,N  
             A(I,J)=ZERO   
         ENDDO
         A(I,I)=ONE    
       ENDDO
C   
C          COMPUTE THE SVD OF THE BIDIAGONAL MATRIX 
C   
      CALL xQRBD (IPASS,S(1),WORK(1,1),NS,A,MDA,N,B,MDB,NB)   
C   
      if(IPASS .eq. 2) then
C         write (*,'(/a)')                                     &
C            ' FULL ACCURACY NOT ATTAINED IN BIDIAGONAL SVD'
          CALL XMESSAGE ('FULL ACCURACY NOT ATTAINED IN BIDIAGONAL SVD')

      endif

  240 CONTINUE  
      IF (NS.GE.N) GO TO 260
      NSP1=NS+1 
          DO J=NSP1,N   
          S(J)=ZERO   
          ENDDO
  260 CONTINUE  
      IF (N.EQ.N1) RETURN   
      NP1=N+1   
C                                  MOVE RECORD OF PERMUTATIONS  
C                                  AND STORE ZEROS  
      DO J=NP1,N1   
          S(J)=A(1,J) 
          DO I=1,N  
             A(I,J)=ZERO   
          ENDDO 
      ENDDO
C                             PERMUTE ROWS AND SET ZERO SINGULAR VALUES.
      DO K=NP1,N1   
          I=S(K)  
          S(K)=ZERO   
          DO J=1,N1 
             A(K,J)=A(I,J) 
             A(I,J)=ZERO   
          ENDDO
          A(I,K)=ONE    
      ENDDO
C                             END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
      RETURN
      END SUBROUTINE xSVDRS

C*************************************************************************C

      double precision FUNCTION xDIFF(X,Y)
C
C  Function used in tests that depend on machine precision.
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 7, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
      double precision X, Y
      xDIFF=X-Y  
      RETURN
      END FUNCTION xDIFF  



C*************************************************************************C

      SUBROUTINE xG1 (A,B,CTERM,STERM,SIG)   
C
C     COMPUTE ORTHOGONAL ROTATION MATRIX..  
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 12, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C   
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))   
C                        (-S,C)         (-S,C)(B)   (   0          )    
C     COMPUTE SIG = SQRT(A**2+B**2) 
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT 
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C     ------------------------------------------------------------------
      double precision A, B, CTERM, ONE, SIG, STERM, XR, YR, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      if (abs(A) .gt. abs(B)) then
         XR=B/A
         YR=sqrt(ONE+XR**2)
         CTERM=sign(ONE/YR,A)
         STERM=CTERM*XR
         SIG=abs(A)*YR     
         RETURN
      endif

      if (B .ne. ZERO) then
         XR=A/B
         YR=sqrt(ONE+XR**2)
         STERM=sign(ONE/YR,B)
         CTERM=STERM*XR
         SIG=abs(B)*YR     
         RETURN
      endif

      SIG=ZERO  
      CTERM=ZERO  
      STERM=ONE   
      RETURN
      END SUBROUTINE xG1

C*************************************************************************C

C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C   
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE   
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B   
C   
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 12, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
C                     Subroutine Arguments
C
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
C            Householder transformation, or Algorithm H2 to apply a
C            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
C            vector.  IUE is the storage increment between elements.  
C            On exit when MODE = 1, U() and UP contain quantities
C            defining the vector U of the Householder transformation.
C            on entry with MODE = 2, U() and UP should contain
C            quantities previously computed with MODE = 1.  These will
C            not be modified during the entry with MODE = 2.   
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
C            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
C            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
C            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
C            NO OPERATIONS WILL BE DONE ON C(). 
C     ------------------------------------------------------------------
      SUBROUTINE xH12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C     ------------------------------------------------------------------
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
C     double precision U(IUE,M)
      double precision U(IUE,*)
      double precision UP
      parameter(ONE = 1.0d0)
C     ------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=abs(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
C                            &***** CONSTRUCT THE TRANSFORMATION. ******
          DO J=L1,M  
           CL=MAX(abs(U(1,J)),CL)  
          ENDDO
C KARLINE      IF (CL) 130,130,20
      IF (CL .LE. 0) GOTO 130
      CLINV=ONE/CL  
      SM=(U(1,LPIVOT)*CLINV)**2   
      DO J=L1,M  
        SM=SM+(U(1,J)*CLINV)**2 
      ENDDO
      CL=CL*SQRT(SM)   
C KARLINE     IF (U(1,LPIVOT)) 50,50,40     
      IF (U(1,LPIVOT) .LE. 0) GOTO 50
      CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70  
C            &***** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C   
CKARLINE   60 IF (CL) 130,130,70
   60 IF (CL .LE. 0) GOTO 130
   70 IF (NCV.LE.0) RETURN  
      B= UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C   
CKARLINE      IF (B) 80,130,130 
      IF (B .GE. 0) GOTO 130 
      B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*UP
           DO I=L1,M  
             SM=SM+C(I3)*U(1,I)
             I3=I3+ICE 
           ENDDO
CKARLINE          IF (SM) 100,120,100   
          IF (SM .EQ. 0) GOTO 120
          SM=SM*B   
          C(I2)=C(I2)+SM*UP
            DO I=L1,M 
              C(I4)=C(I4)+SM*U(1,I)
              I4=I4+ICE 
            ENDDO  
  120     CONTINUE  
  130 RETURN
      END SUBROUTINE xH12



C##########################################################################
C  SUBROUTNE SVDCMP performs a Singular Value Decomposition on matrix A  ##
C                                                                        ##
C  A = U * W * V                                                         ##
C                                                                        ##
C##########################################################################

C linpack routines

      SUBROUTINE XDSVDC (X, LDX, N, P, S, E, U, LDU, V, LDV, WORK, JOB,   &
     &                   INFO)
C***BEGIN PROLOGUE  DSVDC
C***PURPOSE  Perform the singular value decomposition of a rectangular
C            matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D6
C***TYPE      DOUBLE PRECISION (SSVDC-S, DSVDC-D, CSVDC-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX,
C             SINGULAR VALUE DECOMPOSITION
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     DSVDC is a subroutine to reduce a double precision NxP matrix X
C     by orthogonal transformations U and V to diagonal form.  The
C     diagonal elements S(I) are the singular values of X.  The
C     columns of U are the corresponding left singular vectors,
C     and the columns of V the right singular vectors.
C
C     On Entry
C
C         X         DOUBLE PRECISION(LDX,P), where LDX .GE. N.
C                   X contains the matrix whose singular value
C                   decomposition is to be computed.  X is
C                   destroyed by DSVDC.
C
C         LDX       INTEGER.
C                   LDX is the leading dimension of the array X.
C
C         N         INTEGER.
C                   N is the number of rows of the matrix X.
C
C         P         INTEGER.
C                   P is the number of columns of the matrix X.
C
C         LDU       INTEGER.
C                   LDU is the leading dimension of the array U.
C                   (See below).
C
C         LDV       INTEGER.
C                   LDV is the leading dimension of the array V.
C                   (See below).
C
C         WORK      DOUBLE PRECISION(N).
C                   WORK is a scratch array.
C
C         JOB       INTEGER.
C                   JOB controls the computation of the singular
C                   vectors.  It has the decimal expansion AB
C                   with the following meaning
C
C                        A .EQ. 0    do not compute the left singular
C                                  vectors.
C                        A .EQ. 1    return the N left singular vectors
C                                  in U.
C                        A .GE. 2    return the first MIN(N,P) singular
C                                  vectors in U.
C                        B .EQ. 0    do not compute the right singular
C                                  vectors.
C                        B .EQ. 1    return the right singular vectors
C                                  in V.
C
C     On Return
C
C         S         DOUBLE PRECISION(MM), where MM=MIN(N+1,P).
C                   The first MIN(N,P) entries of S contain the
C                   singular values of X arranged in descending
C                   order of magnitude.
C
C         E         DOUBLE PRECISION(P).
C                   E ordinarily contains zeros.  However see the
C                   discussion of INFO for exceptions.
C
C         U         DOUBLE PRECISION(LDU,K), where LDU .GE. N.
C                   If JOBA .EQ. 1, then K .EQ. N.
C                   If JOBA .GE. 2, then K .EQ. MIN(N,P).
C                   U contains the matrix of right singular vectors.
C                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
C                   or if JOBA .EQ. 2, then U may be identified with X
C                   in the subroutine call.
C
C         V         DOUBLE PRECISION(LDV,P), where LDV .GE. P.
C                   V contains the matrix of right singular vectors.
C                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
C                   then V may be identified with X in the
C                   subroutine call.
C
C         INFO      INTEGER.
C                   The singular values (and their corresponding
C                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
C                   are correct (here M=MIN(N,P)).  Thus if
C                   INFO .EQ. 0, all the singular values and their
C                   vectors are correct.  In any event, the matrix
C                   B = TRANS(U)*X*V is the bidiagonal matrix
C                   with the elements of S on its diagonal and the
C                   elements of E on its super-diagonal (TRANS(U)
C                   is the transpose of U).  Thus the singular
C                   values of X and B are the same.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  XDAXPY, XDDOT,XDNRM2, XDROT, XDROTG, XDSCAL, XDSWAP
C***REVISION HISTORY  (YYMMDD)
C   790319  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSVDC
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      DOUBLE PRECISION X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
C
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,    &
     &        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      DOUBLE PRECISION xDDOT,T
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,XDNRM2,SCALE,SHIFT,SL,SM,SN,   &
     &                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
C***FIRST EXECUTABLE STATEMENT  DSVDC
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
      MAXIT = 30
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN(N-1,P)
      NRT = MAX(0,MIN(P-2,N))
      LU = MAX(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = XDNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = SIGN(S(L),X(L,L))
               CALL XDSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -XDDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL XDAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = XDNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = SIGN(E(L),E(LP1))
               CALL XDSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL XDAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL XDAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -XDDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL XDAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL XDSCAL(N-L+1,-1.0D0,U(L,L),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
                  T = -XDDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL XDAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
            IF (L .EQ. 0) GO TO 400
            TEST = ABS(S(L)) + ABS(S(L+1))
            ZTEST = TEST + ABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + ABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + ABS(E(LS-1))
               ZTEST = TEST + ABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C

         GO TO (490,520,540,570), KASE

C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL XDROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL XDROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL XDROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL XDROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = MAX(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),                &
     &                    ABS(S(L)),ABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) - SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL XDROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL XDROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL XDROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)                                   &
     &            CALL XDROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL XDSCAL(P,-1.0D0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)                                  &
     &            CALL XDSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)                                  &
     &            CALL XDSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END



C*************************************************************************C

      DOUBLE PRECISION FUNCTION xDNRM2 (N, DX, INCX)
C***BEGIN PROLOGUE  xDNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, xDNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    xDNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDNRM2
      INTEGER NEXT, N,NN,INCX,I,J
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,      &
     &                 ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  xDNRM2
      IF (N .GT. 0) GO TO 10
         XDNRM2  = ZERO
         GO TO 300
C
C   10 ASSIGN 30 TO NEXT               CHANGED INTO:
   10 NEXT = 30      
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
C  20    GO TO NEXT,(30, 50, 70, 110)  CHANGED INTO:
   20 SELECT CASE (NEXT)

        CASE (30)
          GOTO 30
        CASE (50)
          GOTO 50
        CASE (70)
          GOTO 70
        CASE (110)
          GOTO 110
                   
      END SELECT
      
      
        
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C     ASSIGN 50 TO NEXT               CHANGED INTO:
      NEXT = 50

      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
C      ASSIGN 70 TO NEXT              CHANGED INTO:
      NEXT = 70
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
C      ASSIGN 110 TO NEXT              CHANGED INTO:
      NEXT = 110
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
C KARLINE

C      DO 95 J = I,NN,INCX
C      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
C   95    SUM = SUM + DX(J)**2

      DO J = I,NN,INCX
         IF (ABS(DX(J)) .GE. HITEST) GO TO 100
         SUM = SUM + DX(J)**2
      ENDDO 

      XDNRM2 = SQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      XDNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END

C*************************************************************************C

CDECK DROT
      SUBROUTINE XDROT (N, DX, INCX, DY, INCY, DC, DS)
C***BEGIN PROLOGUE  DROT
C***PURPOSE  Apply a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      DOUBLE PRECISION (SROT-S, DROT-D, CSROT-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C       DC  D.P. element of rotation matrix
C       DS  D.P. element of rotation matrix
C
C     --Output--
C       DX  rotated vector DX (unchanged if N .LE. 0)
C       DY  rotated vector DY (unchanged if N .LE. 0)
C
C     Multiply the 2 x 2 matrix  ( DC DS) times the 2 x N matrix (DX**T)
C                                (-DS DC)                        (DY**T)
C     where **T indicates transpose.  The elements of DX are in
C     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DROT
      DOUBLE PRECISION DX, DY, DC, DS, ZERO, ONE, W, Z
      DIMENSION DX(*), DY(*)
      INTEGER N,INCX,INCY,NSTEPS,I,KX,KY
      SAVE ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C***FIRST EXECUTABLE STATEMENT  DROT
      IF (N .LE. 0 .OR. (DS .EQ. ZERO .AND. DC .EQ. ONE)) GO TO 40
      IF (.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
C          Code for equal and positive increments.
C
           NSTEPS=INCX*N
           DO 10 I = 1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=DC*W+DS*Z
                DY(I)=-DS*W+DC*Z
   10           CONTINUE
           GO TO 40
C
C     Code for unequal or nonpositive increments.
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF (INCX .LT. 0) KX = 1-(N-1)*INCX
           IF (INCY .LT. 0) KY = 1-(N-1)*INCY
C
           DO 30 I = 1,N
                W=DX(KX)
                Z=DY(KY)
                DX(KX)=DC*W+DS*Z
                DY(KY)=-DS*W+DC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END

C*************************************************************************C

      SUBROUTINE xDROTG (DA, DB, DC, DS)
C***BEGIN PROLOGUE  DROTG
C***PURPOSE  Construct a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      DOUBLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C       DA  double precision scalar
C       DB  double precision scalar
C
C     --Output--
C       DA  double precision result R
C       DB  double precision result Z
C       DC  double precision result
C       DS  double precision result
C
C     Construct the Givens transformation
C
C         ( DC  DS )
C     G = (        ) ,    DC**2 + DS**2 = 1 ,
C         (-DS  DC )
C
C     which zeros the second entry of the 2-vector  (DA,DB)**T .
C
C     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
C     storage.  The value of DB is overwritten by a value Z which
C     allows DC and DS to be recovered by the following algorithm.
C
C           If Z=1  set  DC=0.0  and  DS=1.0
C           If ABS(Z) .LT. 1  set  DC=SQRT(1-Z**2)  and  DS=Z
C           If ABS(Z) .GT. 1  set  DC=1/Z  and  DS=SQRT(1-DC**2)
C
C     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will
C     next be called to apply the transformation to a 2 by N matrix.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DROTG
      DOUBLE PRECISION  DA, DB, DC, DS, U, V, R
C***FIRST EXECUTABLE STATEMENT  DROTG
      IF (ABS(DA) .LE. ABS(DB)) GO TO 10
C
C *** HERE ABS(DA) .GT. ABS(DB) ***
C
      U = DA + DA
      V = DB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DA
C
      R = SQRT(0.25D0 + V**2) * U
C
C     NOTE THAT DC IS POSITIVE
C
      DC = DA / R
      DS = V * (DC + DC)
      DB = DS
      DA = R
      RETURN
C
C *** HERE ABS(DA) .LE. ABS(DB) ***
C
   10 IF (DB .EQ. 0.0D0) GO TO 20
      U = DB + DB
      V = DA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DB
C     (R IS IMMEDIATELY STORED IN DA)
C
      DA = SQRT(0.25D0 + V**2) * U
C
C     NOTE THAT DS IS POSITIVE
C
      DS = DB / DA
      DC = V * (DS + DS)
      IF (DC .EQ. 0.0D0) GO TO 15
      DB = 1.0D0 / DC
      RETURN
   15 DB = 1.0D0
      RETURN
C
C *** HERE DA = DB = 0.0 ***
C
   20 DC = 1.0D0
      DS = 0.0D0
      RETURN
C
      END


C########################################################################
C  SUBROUTNE BackSubstitution calculates X:                            ##
C   A * X = B                                                          ##
C   A = U * W * VT                                                     ##
C   X = V * 1/W *UT * B                                                ##
C                                                                      ##
C   and X is the desired simplest unconstrained solution               ##
C                                                                      ##
C########################################################################

      SUBROUTINE xBackSubstitution (UT,W,V,B,m,n,MX,NX,X,NMX)
      IMPLICIT NONE
      INTEGER            :: m,n,MX,NX,NMX 
      DOUBLE PRECISION   :: B(MX),UT(MX,MX),V(NX,NX),W(NMX),X(NX),S(NMX)
      INTEGER I,J

      DO I = 1,nx
      X(I) = 0.d0
      ENDDO
      DO I = 1,NMX
      S(I) = 0.d0
      ENDDO
C Calculate (UT*B)/w, only if W(I) <> 0 (else 'infinity')

      DO I=1,M
       IF (W(I) <= 0.D0) EXIT
        DO J=1,M
          S(I) = S(I) + UT(I,J)*B(J)
        ENDDO
      ENDDO

      DO I=1,NMX
      IF (W(I) <= 0.D0) EXIT
       S(I) = S(I)/W(I)
      ENDDO

      DO I=1,N
       DO J=1,N
       X(I) = X(I) + V(I,J)*S(J)
       ENDDO
      ENDDO

      RETURN
      END







C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C     SOLVING SYSTEMS OF EQUATIONS   C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C



C********************************************************************

      subroutine  xdscal(n,da,dx,incx)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     scales a vector by a constant.
C     uses unrolled loops for increment equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C***BEGIN PROLOGUE  xDSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, xDSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDSCAL

      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
C
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
C
C        code for increment not equal to 1
C
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end                                        
C of xDscal




C********************************************************************
      subroutine xdaxpy(n,da,dx,incx,dy,incy)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C***BEGIN PROLOGUE  xDAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, xDAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDAXPY

C
C     constant times a vector plus a vector.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
C
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C        code for unequal increments or equal increments
C          not equal to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end                                        
C of xdaxpy

C********************************************************************

      double precision function xddot(n,dx,incx,dy,incy)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     forms the dot product of two vectors.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C
      xddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C        code for unequal increments or equal increments
C          not equal to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      xddot = dtemp
      return
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +                &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 xddot = dtemp
      return
      end                                        
C of xDDOT



      subroutine xdswap (n,dx,incx,dy,incy)
C
C     interchanges two vectors.
C     uses unrolled loops for increments equal one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C       code for unequal increments or equal increments not equal
C         to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
C       code for both increments equal to 1
C
C
C       clean-up loop
C
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end


C***********************************************************************
      SUBROUTINE XDCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END


C********************************************************************
      integer function xidamax(n,dx,incx)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     finds the index of element having max. absolute value.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dmax
      integer i,incx,ix,n
C
      xidamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      xidamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
C
C        code for increment not equal to 1
C
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         xidamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
C
C        code for increment equal to 1
C
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         xidamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end                                        
C of xidamax


C***********************************************************************

      DOUBLE PRECISION FUNCTION D1MACH (IDUM)
      INTEGER  :: IDUM
C-----------------------------------------------------------------------
C This routine computes:
C the unit roundoff of the machine when IDUM = 4
C the largest value                when IDUM = 2
C Unit roundoffis defined as the smallest positive machine number
C u such that  1.0 + u .ne. 1.0
C Largest value is simply imposed...
C Subroutines/functions called by D1MACH.. None
C-----------------------------------------------------------------------
      DOUBLE PRECISION U, COMP
      DOUBLE PRECISION :: Prec(4)  
      LOGICAL          :: First(4) 
      SAVE Prec, FIRST
      DATA FIRST /.TRUE.,.TRUE.,.TRUE.,.TRUE./
      DATA Prec /1.D-8,1.D-8,1.D-8,1.D-8/


      IF (Idum > 4 .OR. Idum < 0) THEN
C         Write (*,*) "Error in function D1MACH"
C         Write (*,*) "NOT DEFINED FOR IDUM = ", Idum
       CALL XMESSAGE("Error in function D1MACH-NOT DEFINED FOR IDUM  ") 
      ENDIF

      IF (First(Idum)) THEN 

       First(Idum) = .FALSE.

       SELECT CASE (IDUM)

        CASE (2)
C Very large number
         D1MACH = 1.D300

        CASE (4)
C Unit roundoff
         U = 1.0D0
 10      U = U*0.5D0
         COMP = 1.0D0 + U
         IF (COMP .NE. 1.0D0) GO TO 10
         D1MACH = U*2.0D0
        CASE Default
C         Write (*,*) "Error in function D1MACH"
C         Write (*,*) "NOT DEFINED FOR IDUM = ", Idum
         CALL XMESSAGE("Error in function D1MACH-NOT DEFINED FOR IDUM ")
       END SELECT

       PREC (Idum) = D1MACH

      ELSE

       D1mach = Prec(IDUM)

      ENDIF

      RETURN

      END FUNCTION d1Mach





C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C    LEAST SQUARES with constraints  C
C                DLSEI               C                 
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C


C   LINPACK routine





                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
                !    LEAST SQUARES with constraints  !
                !                DLSEI               !                 
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!


!   LINPACK routine

      SUBROUTINE xDLSEI (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME,       &
     &   RNORML, MODE, WS, IP)
!***BEGIN PROLOGUE  xDLSEI
!***PURPOSE  Solve a linearly constrained least squares problem with
!            equality and inequality constraints, and optionally compute
!            a covariance matrix.
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A, D9
!***TYPE      DOUBLE PRECISION (LSEI-S, xDLSEI-D)
!***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
!             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
!             QUADRATIC PROGRAMMING
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Abstract
!
!     This subprogram solves a linearly constrained least squares
!     problem with both equality and inequality constraints, and, if the
!     user requests, obtains a covariance matrix of the solution
!     parameters.
!
!     Suppose there are given matrices E, A and G of respective
!     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
!     respective lengths ME, MA and MG.  This subroutine solves the
!     linearly constrained least squares problem
!
!                   EX = F, (E ME by N) (equations to be exactly
!                                       satisfied)
!                   AX = B, (A MA by N) (equations to be
!                                       approximately satisfied,
!                                       least squares sense)
!                   GX .GE. H,(G MG by N) (inequality constraints)
!
!     The inequalities GX .GE. H mean that every component of the
!     product GX must be .GE. the corresponding component of H.
!
!     In case the equality constraints cannot be satisfied, a
!     generalized inverse solution residual vector length is obtained
!     for F-EX.  This is the minimal length possible for F-EX.
!
!     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
!     rank of the matrix E is estimated during the computation.  We call
!     this value KRANKE.  It is an output parameter in IP(1) defined
!     below.  Using a generalized inverse solution of EX=F, a reduced
!     least squares problem with inequality constraints is obtained.
!     The tolerances used in these tests for determining the rank
!     of E and the rank of the reduced least squares problem are
!     given in Sandia Tech. Rept. SAND-78-1290.  They can be
!     modified by the user if new values are provided in
!     the option list of the array PRGOPT(*).
!
!     The user must dimension all arrays appearing in the call list..
!     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
!     where K=MAX(MA+MG,N).  This allows for a solution of a range of
!     problems in the given working space.  The dimension of WS(*)
!     given is a necessary overestimate.  Once a particular problem
!     has been run, the output parameter IP(3) gives the actual
!     dimension required for that problem.
!
!     The parameters for xDLSEI( ) are
!
!     Input.. All TYPE REAL variables are DOUBLE PRECISION
!
!     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
!     ME,MA,MG,N    first dimensioning parameter equal to MDW.
!                   For this discussion let us call M = ME+MA+MG.  Then
!                   MDW must satisfy MDW .GE. M.  The condition
!                   MDW .LT. M is an error.
!
!                   The array W(*,*) contains the matrices and vectors
!
!                                  (E  F)
!                                  (A  B)
!                                  (G  H)
!
!                   in rows and columns 1,...,M and 1,...,N+1
!                   respectively.
!
!                   The integers ME, MA, and MG are the
!                   respective matrix row dimensions
!                   of E, A and G.  Each matrix has N columns.
!
!     PRGOPT(*)    This real-valued array is the option vector.
!                  If the user is satisfied with the nominal
!                  subprogram features set
!
!                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!                  Otherwise PRGOPT(*) is a linked list consisting of
!                  groups of data of the following form
!
!                  LINK
!                  KEY
!                  DATA SET
!
!                  The parameters LINK and KEY are each one word.
!                  The DATA SET can be comprised of several words.
!                  The number of items depends on the value of KEY.
!                  The value of LINK points to the first
!                  entry of the next group of data within
!                  PRGOPT(*).  The exception is when there are
!                  no more options to change.  In that
!                  case, LINK=1 and the values KEY and DATA SET
!                  are not referenced.  The general layout of
!                  PRGOPT(*) is as follows.
!
!               ...PRGOPT(1) = LINK1 (link to first entry of next group)
!               .  PRGOPT(2) = KEY1 (key to the option change)
!               .  PRGOPT(3) = data value (data value for this change)
!               .       .
!               .       .
!               .       .
!               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
!               .                       next group)
!               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
!               .  PRGOPT(LINK1+2) = data value
!               ...     .
!               .       .
!               .       .
!               ...PRGOPT(LINK) = 1 (no more options to change)
!
!                  Values of LINK that are nonpositive are errors.
!                  A value of LINK .GT. NLINK=100000 is also an error.
!                  This helps prevent using invalid but positive
!                  values of LINK that will probably extend
!                  beyond the program limits of PRGOPT(*).
!                  Unrecognized values of KEY are ignored.  The
!                  order of the options is arbitrary and any number
!                  of options can be changed with the following
!                  restriction.  To prevent cycling in the
!                  processing of the option array, a count of the
!                  number of options changed is maintained.
!                  Whenever this count exceeds NOPT=1000, an error
!                  message is printed and the subprogram returns.
!
!                  Options..
!
!                  KEY=1
!                         Compute in W(*,*) the N by N
!                  covariance matrix of the solution variables
!                  as an output parameter.  Nominally the
!                  covariance matrix will not be computed.
!                  (This requires no user input.)
!                  The data set for this option is a single value.
!                  It must be nonzero when the covariance matrix
!                  is desired.  If it is zero, the covariance
!                  matrix is not computed.  When the covariance matrix
!                  is computed, the first dimensioning parameter
!                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
!
!                  KEY=10
!                         Suppress scaling of the inverse of the
!                  normal matrix by the scale factor RNORM**2/
!                  MAX(1, no. of degrees of freedom).  This option
!                  only applies when the option for computing the
!                  covariance matrix (KEY=1) is used.  With KEY=1 and
!                  KEY=10 used as options the unscaled inverse of the
!                  normal matrix is returned in W(*,*).
!                  The data set for this option is a single value.
!                  When it is nonzero no scaling is done.  When it is
!                  zero scaling is done.  The nominal case is to do
!                  scaling so if option (KEY=1) is used alone, the
!                  matrix will be scaled on output.
!
!                  KEY=2
!                         Scale the nonzero columns of the
!                         entire data matrix.
!                  (E)
!                  (A)
!                  (G)
!
!                  to have length one.  The data set for this
!                  option is a single value.  It must be
!                  nonzero if unit length column scaling
!                  is desired.
!
!                  KEY=3
!                         Scale columns of the entire data matrix
!                  (E)
!                  (A)
!                  (G)
!
!                  with a user-provided diagonal matrix.
!                  The data set for this option consists
!                  of the N diagonal scaling factors, one for
!                  each matrix column.
!
!                  KEY=4
!                         Change the rank determination tolerance for
!                  the equality constraint equations from
!                  the nominal value of SQRT(DRELPR).  This quantity can
!                  be no smaller than DRELPR, the arithmetic-
!                  storage precision.  The quantity DRELPR is the
!                  largest positive number such that T=1.+DRELPR
!                  satisfies T .EQ. 1.  The quantity used
!                  here is internally restricted to be at
!                  least DRELPR.  The data set for this option
!                  is the new tolerance.
!
!                  KEY=5
!                         Change the rank determination tolerance for
!                  the reduced least squares equations from
!                  the nominal value of SQRT(DRELPR).  This quantity can
!                  be no smaller than DRELPR, the arithmetic-
!                  storage precision.  The quantity used
!                  here is internally restricted to be at
!                  least DRELPR.  The data set for this option
!                  is the new tolerance.
!
!                  For example, suppose we want to change
!                  the tolerance for the reduced least squares
!                  problem, compute the covariance matrix of
!                  the solution parameters, and provide
!                  column scaling for the data matrix.  For
!                  these options the dimension of PRGOPT(*)
!                  must be at least N+9.  The Fortran statements
!                  defining these options would be as follows:
!
!                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
!                  PRGOPT(2)=1 (covariance matrix key)
!                  PRGOPT(3)=1 (covariance matrix wanted)
!
!                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
!                  PRGOPT(5)=5 (least squares equas.  tolerance key)
!                  PRGOPT(6)=... (new value of the tolerance)
!
!                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
!                  PRGOPT(8)=3 (user-provided column scaling key)
!
!                  CALL xDCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
!                    scaling factors from the user array D(*)
!                    to PRGOPT(9)-PRGOPT(N+8))
!
!                  PRGOPT(N+9)=1 (no more options to change)
!
!                  The contents of PRGOPT(*) are not modified
!                  by the subprogram.
!                  The options for WNNLS( ) can also be included
!                  in this array.  The values of KEY recognized
!                  by WNNLS( ) are 6, 7 and 8.  Their functions
!                  are documented in the usage instructions for
!                  subroutine WNNLS( ).  Normally these options
!                  do not need to be modified when using xDLSEI( ).
!
!     IP(1),       The amounts of working storage actually
!     IP(2)        allocated for the working arrays WS(*) and
!                  IP(*), respectively.  These quantities are
!                  compared with the actual amounts of storage
!                  needed by xDLSEI( ).  Insufficient storage
!                  allocated for either WS(*) or IP(*) is an
!                  error.  This feature was included in xDLSEI( )
!                  because miscalculating the storage formulas
!                  for WS(*) and IP(*) might very well lead to
!                  subtle and hard-to-find execution errors.
!
!                  The length of WS(*) must be at least
!
!                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
!
!                  where K = max(MA+MG,N)
!                  This test will not be made if IP(1).LE.0.
!
!                  The length of IP(*) must be at least
!
!                  LIP = MG+2*N+2
!                  This test will not be made if IP(2).LE.0.
!
!     Output.. All TYPE REAL variables are DOUBLE PRECISION
!
!     X(*),RNORME,  The array X(*) contains the solution parameters
!     RNORML        if the integer output flag MODE = 0 or 1.
!                   The definition of MODE is given directly below.
!                   When MODE = 0 or 1, RNORME and RNORML
!                   respectively contain the residual vector
!                   Euclidean lengths of F - EX and B - AX.  When
!                   MODE=1 the equality constraint equations EX=F
!                   are contradictory, so RNORME .NE. 0.  The residual
!                   vector F-EX has minimal Euclidean length.  For
!                   MODE .GE. 2, none of these parameters is defined.
!
!     MODE          Integer flag that indicates the subprogram
!                   status after completion.  If MODE .GE. 2, no
!                   solution has been computed.
!
!                   MODE =
!
!                   0  Both equality and inequality constraints
!                      are compatible and have been satisfied.
!
!                   1  Equality constraints are contradictory.
!                      A generalized inverse solution of EX=F was used
!                      to minimize the residual vector length F-EX.
!                      In this sense, the solution is still meaningful.
!
!                   2  Inequality constraints are contradictory.
!
!                   3  Both equality and inequality constraints
!                      are contradictory.
!
!                   The following interpretation of
!                   MODE=1,2 or 3 must be made.  The
!                   sets consisting of all solutions
!                   of the equality constraints EX=F
!                   and all vectors satisfying GX .GE. H
!                   have no points in common.  (In
!                   particular this does not say that
!                   each individual set has no points
!                   at all, although this could be the
!                   case.)
!
!                   4  Usage error occurred.  The value
!                      of MDW is .LT. ME+MA+MG, MDW is
!                      .LT. N and a covariance matrix is
!                      requested, or the option vector
!                      PRGOPT(*) is not properly defined,
!                      or the lengths of the working arrays
!                      WS(*) and IP(*), when specified in
!                      IP(1) and IP(2) respectively, are not
!                      long enough.
!
!     W(*,*)        The array W(*,*) contains the N by N symmetric
!                   covariance matrix of the solution parameters,
!                   provided this was requested on input with
!                   the option vector PRGOPT(*) and the output
!                   flag is returned with MODE = 0 or 1.
!
!     IP(*)         The integer working array has three entries
!                   that provide rank and working array length
!                   information after completion.
!
!                      IP(1) = rank of equality constraint
!                              matrix.  Define this quantity
!                              as KRANKE.
!
!                      IP(2) = rank of reduced least squares
!                              problem.
!
!                      IP(3) = the amount of storage in the
!                              working array WS(*) that was
!                              actually used by the subprogram.
!                              The formula given above for the length
!                              of WS(*) is a necessary overestimate.
!                              If exactly the same problem matrices
!                              are used in subsequent executions,
!                              the declared dimension of WS(*) can
!                              be reduced to this output value.
!     User Designated
!     Working Arrays..
!
!     WS(*),IP(*)              These are respectively type real
!                              and type integer working arrays.
!                              Their required minimal lengths are
!                              given above.
!
!***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Report SAND77-0552, Sandia
!                 Laboratories, June 1978.
!               K. H. Haskell and R. J. Hanson, Selected algorithms for
!                 the linearly constrained least squares problem - a
!                 users guide, Report SAND78-1290, Sandia Laboratories,
!                 August 1979.
!               K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Mathematical Programming
!                 21 (1981), pp. 98-118.
!               R. J. Hanson and K. H. Haskell, Two algorithms for the
!                 linearly constrained least squares problem, ACM
!                 Transactions on Mathematical Software, September 1982.
!***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDDOT, xDH12, DLSI,
!                    xDNRM2, xDSCAL, xDSWAP, xXERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
!   900510  Convert XERRWV calls to xXERMSG calls.  (RWC)
!   900604  DP version created from SP version.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  xDLSEI
      INTEGER IP(3), MA, MDW, ME, MG, MODE, N
      DOUBLE PRECISION PRGOPT(*), RNORME, RNORML, W(MDW,*), WS(*), X(*)
!
      EXTERNAL D1MACH, xDASUM, xDAXPY, xDCOPY,xDDOT,xDH12,DLSI,xDNRM2,               &
     &   xDSCAL, xDSWAP, xXERMSG
      DOUBLE PRECISION D1MACH, xDASUM, xDDOT, xDNRM2
!
      DOUBLE PRECISION DRELPR, ENORM, FNORM, GAM, RB, RN, RNMAX, SIZE,             &
     &   SN, SNMAX, T, TAU, UJ, UP, VJ, XNORM, XNRME
      INTEGER I, IMAX, J, JP1, K, KEY, KRANKE, LAST, LCHK, LINK, M,                &
     &   MAPKE1, MDEQC, MEND, MEP1, N1, N2, NEXT, NLINK, NOPT, NP1,                &
     &   NTIMES
      LOGICAL COV, FIRST
      CHARACTER*8 XERN1, XERN2, XERN3, XERN4
      SAVE FIRST, DRELPR
!
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  xDLSEI
!
!     Set the nominal tolerance used in the code for the equality
!     constraint equations.
!
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
      TAU = SQRT(DRELPR)
!
!     Check that enough storage was allocated in WS(*) and IP(*).
!
      MODE = 4
      IF (MIN(N,ME,MA,MG) .LT. 0) THEN
         WRITE (XERN1, '(I8)') N
         WRITE (XERN2, '(I8)') ME
         WRITE (XERN3, '(I8)') MA
         WRITE (XERN4, '(I8)') MG
         CALL xXERMSG ('SLATEC', 'LSEI', 'ALL OF THE VARIABLES N, ME,'//            &
     &      ' MA, MG MUST BE .GE. 0 ENTERED ROUTINE WITH' //                        &
     &      ' N  = ' // XERN1 //                                                    &
     &      ' ME = ' // XERN2 //                                                    &
     &      ' MA = ' // XERN3 //                                                    &
     &      ' MG = ' // XERN4, 2, 1)
         RETURN
      ENDIF
!
      IF (IP(1).GT.0) THEN
         LCHK = 2*(ME+N) + MAX(MA+MG,N) + (MG+2)*(N+7)
         IF (IP(1).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL xXERMSG ('SLATEC', 'xDLSEI', 'INSUFFICIENT STORAGE ' //             &
     &         'ALLOCATED FOR WS(*), NEED LW = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
!
      IF (IP(2).GT.0) THEN
         LCHK = MG + 2*N + 2
         IF (IP(2).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL xXERMSG ('SLATEC', 'xDLSEI', 'INSUFFICIENT STORAGE ' //             &
     &         'ALLOCATED FOR IP(*), NEED LIP = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
!
!     Compute number of possible right multiplying Householder
!     transformations.
!
      M = ME + MA + MG
      IF (N.LE.0 .OR. M.LE.0) THEN
         MODE = 0
         RNORME = 0
         RNORML = 0
         RETURN
      ENDIF
!
      IF (MDW.LT.M) THEN
        CALL xXERMSG ('SLATEC', 'xDLSEI', 'MDW.LT.ME+MA+MG IS AN ERROR',            &
     &      2, 1)
         RETURN
      ENDIF
!
      NP1 = N + 1
      KRANKE = MIN(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
!
!     Set nominal values.
!
!     The nominal column scaling used in the code is
!     the identity scaling.
!
      CALL xDCOPY (N, 1.D0, 0, WS(N1), 1)
!
!     No covariance matrix is nominally computed.
!
      COV = .FALSE.
!
!     Process option vector.
!     Define bound for number of options to change.
!
      NOPT = 1000
      NTIMES = 0
!
!     Define bound for positive values of LINK.
!
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.EQ.0 .OR. LINK.GT.NLINK) THEN
         CALL xXERMSG('SLATEC','xDLSEI','THE OPTION VECTOR IS UNDEFINED'            &
     &   ,2,1)
         RETURN
      ENDIF
!
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
            CALL xXERMSG ('SLATEC','xDLSEI',                                         &
     &         'THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 2, 1)
            RETURN
         ENDIF
!
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) THEN
            COV = PRGOPT(LAST+2) .NE. 0.D0
         ELSEIF (KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = xDNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               WS(J+N1-1) = T
  110       CONTINUE
         ELSEIF (KEY.EQ.3) THEN
            CALL xDCOPY (N, PRGOPT(LAST+2), 1, WS(N1), 1)
         ELSEIF (KEY.EQ.4) THEN
            TAU = MAX(DRELPR,PRGOPT(LAST+2))
         ENDIF
!
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
         CALL xXERMSG ('SLATEC', 'xDLSEI',                                           &
     &      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
            RETURN
         ENDIF
!
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
!
      DO 120 J = 1,N
         CALL xDSCAL (M, WS(N1+J-1), W(1,J), 1)
  120 CONTINUE
!
      IF (COV .AND. MDW.LT.N) THEN
         CALL xXERMSG ('SLATEC', 'xDLSEI',                                           &
     &      'MDW .LT. N WHEN COV MATRIX NEEDED, IS AN ERROR', 2, 1)
         RETURN
      ENDIF
!
!     Problem definition and option vector OK.
!
      MODE = 0
!
!     Compute norm of equality constraint matrix and right side.
!
      ENORM = 0.D0
      DO 130 J = 1,N
         ENORM = MAX(ENORM,xDASUM(ME,W(1,J),1))
  130 CONTINUE
!
      FNORM = xDASUM(ME,W(1,NP1),1)
      SNMAX = 0.D0
      RNMAX = 0.D0
      DO 150 I = 1,KRANKE
!
!        Compute maximum ratio of vector lengths. Partition is at
!        column I.
!
         DO 140 K = I,ME
            SN = xDDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
            RN = xDDOT(I-1,W(K,1),MDW,W(K,1),MDW)
            IF (RN.EQ.0.D0 .AND. SN.GT.SNMAX) THEN
               SNMAX = SN
               IMAX = K
            ELSEIF (K.EQ.I .OR. SN*RNMAX.GT.RN*SNMAX) THEN
               SNMAX = SN
               RNMAX = RN
               IMAX = K
            ENDIF
  140    CONTINUE
!
!        Interchange rows if necessary.
!
         IF (I.NE.IMAX) CALL xDSWAP (NP1, W(I,1), MDW, W(IMAX,1), MDW)
         IF (SNMAX.GT.RNMAX*TAU**2) THEN
!
!        Eliminate elements I+1,...,N in row I.
!
            CALL xDH12 (1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW,             &    
     &                1, M-I)
         ELSE
            KRANKE = I - 1
            GO TO 160
         ENDIF
  150 CONTINUE
!
!     Save diagonal terms of lower trapezoidal matrix.
!
  160 CALL xDCOPY (KRANKE, W, MDW+1, WS(KRANKE+1), 1)
!
!     Use Householder transformation from left to achieve
!     KRANKE by KRANKE upper triangular form.
!
      IF (KRANKE.LT.ME) THEN
         DO 170 K = KRANKE,1,-1
!
!           Apply transformation to matrix cols. 1,...,K-1.
!
         CALL xDH12 (1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW,K-1)
!
!           Apply to rt side vector.
!
        CALL xDH12 (2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1),1,1, 1)
  170    CONTINUE
      ENDIF
!
!     Solve for variables 1,...,KRANKE in new coordinates.
!
      CALL xDCOPY (KRANKE, W(1, NP1), 1, X, 1)
      DO 180 I = 1,KRANKE
         X(I) = (X(I)-xDDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  180 CONTINUE
!
!     Compute residuals for reduced problem.
!
      MEP1 = ME + 1
      RNORML = 0.D0
      DO 190 I = MEP1,M
         W(I,NP1) = W(I,NP1) - xDDOT(KRANKE,W(I,1),MDW,X,1)
         SN = xDDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
         RN = xDDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
         IF (RN.LE.SN*TAU**2 .AND. KRANKE.LT.N)                                     &
     &      CALL xDCOPY (N-KRANKE, 0.D0, 0, W(I,KRANKE+1), MDW)
  190 CONTINUE
!
!     Compute equality constraint equations residual length.
!
      RNORME = xDNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
!
!     Move reduced problem data upward if KRANKE.LT.ME.
!
      IF (KRANKE.LT.ME) THEN
         DO 200 J = 1,NP1
            CALL xDCOPY (M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  200    CONTINUE
      ENDIF
!
!     Compute solution of reduced problem.
!
      CALL DLSI(W(KRANKE+1, KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,               &
     &         X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
!
!     Test for consistency of equality constraints.
!
      IF (ME.GT.0) THEN
         MDEQC = 0
         XNRME = xDASUM(KRANKE,W(1,NP1),1)
         IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
         MODE = MODE + MDEQC
!
!        Check if solution to equality constraints satisfies inequality
!        constraints when there are no degrees of freedom left.
!
         IF (KRANKE.EQ.N .AND. MG.GT.0) THEN
            XNORM = xDASUM(N,X,1)
            MAPKE1 = MA + KRANKE + 1
            MEND = MA + KRANKE + MG
            DO 210 I = MAPKE1,MEND
               SIZE = xDASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
               IF (W(I,NP1).GT.TAU*SIZE) THEN
                  MODE = MODE + 2
                  GO TO 290
               ENDIF
  210       CONTINUE
         ENDIF
      ENDIF
!
!     Replace diagonal terms of lower trapezoidal matrix.
!
      IF (KRANKE.GT.0) THEN
         CALL xDCOPY (KRANKE, WS(KRANKE+1), 1, W, MDW+1)
!
!        Reapply transformation to put solution in original coordinates.
!
         DO 220 I = KRANKE,1,-1
            CALL xDH12 (2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  220    CONTINUE
!
!        Compute covariance matrix of equality constrained problem.
!
         IF (COV) THEN
            DO 270 J = MIN(KRANKE,N-1),1,-1
               RB = WS(J)*W(J,J)
               IF (RB.NE.0.D0) RB = 1.D0/RB
               JP1 = J + 1
               DO 230 I = JP1,N
                  W(I,J) = RB*xDDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)
  230          CONTINUE
!
               GAM = 0.5D0*RB*xDDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)
               CALL xDAXPY (N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
               DO 250 I = JP1,N
                  DO 240 K = I,N
                     W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
                     W(K,I) = W(I,K)
  240             CONTINUE
  250          CONTINUE
               UJ = WS(J)
               VJ = GAM*UJ
               W(J,J) = UJ*VJ + UJ*VJ
               DO 260 I = JP1,N
                  W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  260          CONTINUE
               CALL xDCOPY (N-J, W(J, JP1), MDW, W(JP1,J), 1)
  270       CONTINUE
         ENDIF
      ENDIF
!
!     Apply the scaling to the covariance matrix.
!
      IF (COV) THEN
         DO 280 I = 1,N
            CALL xDSCAL (N, WS(I+N1-1), W(I,1), MDW)
            CALL xDSCAL (N, WS(I+N1-1), W(1,I), 1)
  280    CONTINUE
      ENDIF
!
!     Rescale solution vector.
!
  290 IF (MODE.LE.1) THEN
         DO 300 J = 1,N
            X(J) = X(J)*WS(N1+J-1)
  300    CONTINUE
      ENDIF
!
      IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
      END





      SUBROUTINE DLSI (W, MDW, MA, MG, N, PRGOPT, X, RNORM, MODE, WS,               &
     &   IP)
!***BEGIN PROLOGUE  DLSI
!***SUBSIDIARY
!***PURPOSE  Subsidiary to xDLSEI
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (LSI-S, DLSI-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to xDLSEI.  The documentation for
!     xDLSEI has complete usage instructions.
!
!     Solve..
!              AX = B,  A  MA by N  (least squares equations)
!     subject to..
!
!              GX.GE.H, G  MG by N  (inequality constraints)
!
!     Input..
!
!      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
!                       (G H)
!
!     MDW,MA,MG,N
!              contain (resp) var. dimension of W(*,*),
!              and matrix dimensions.
!
!     PRGOPT(*),
!              Program option vector.
!
!     OUTPUT..
!
!      X(*),RNORM
!
!              Solution vector(unless MODE=2), length of AX-B.
!
!      MODE
!              =0   Inequality constraints are compatible.
!              =2   Inequality constraints contradictory.
!
!      WS(*),
!              Working storage of dimension K+N+(MG+2)*(N+7),
!              where K=MAX(MA+MG,N).
!      IP(MG+2*N+1)
!              Integer working storage
!
!***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDDOT, xDH12, xDHFTI,
!                    DLPDP, xDSCAL, xDSWAP
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900604  DP version created from SP version.  (RWC)
!   920422  Changed CALL to xDHFTI to include variable MA.  (WRB)
!***END PROLOGUE  DLSI
      INTEGER IP(*), MA, MDW, MG, MODE, N
      DOUBLE PRECISION PRGOPT(*), RNORM, W(MDW,*), WS(*), X(*)
!
      EXTERNAL D1MACH, xDASUM,xDAXPY,xDCOPY,xDDOT,xDH12,xDHFTI,DLPDP,               &
     &   xDSCAL, xDSWAP
      DOUBLE PRECISION D1MACH, xDASUM, xDDOT
!
      DOUBLE PRECISION ANORM, DRELPR, FAC, GAM, RB, TAU, TOL, XNORM
      INTEGER I, J, K, KEY, KRANK, KRM1, KRP1, L, LAST, LINK, M, MAP1,              &
     &   MDLPDP, MINMAN, N1, N2, N3, NEXT, NP1, MDB ! KARLINE: ADDED MDB...
      LOGICAL COV, FIRST, SCLCOV
!
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DLSI
!
!     Set the nominal tolerance used in the code.
!
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
      TOL = SQRT(DRELPR)
!
      MODE = 0
      RNORM = 0.D0
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 370
!
!     To process option vector.
!
      COV = .FALSE.
      SCLCOV = .TRUE.
      LAST = 1
      LINK = PRGOPT(1)
!
  100 IF (LINK.GT.1) THEN
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) COV = PRGOPT(LAST+2) .NE. 0.D0
         IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2) .EQ. 0.D0
         IF (KEY.EQ.5) TOL = MAX(DRELPR,PRGOPT(LAST+2))
         NEXT = PRGOPT(LINK)
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
!
!     Compute matrix norm of least squares equations.
!
      ANORM = 0.D0
      DO 110 J = 1,N
         ANORM = MAX(ANORM,xDASUM(MA,W(1,J),1))
  110 CONTINUE
!
!     Set tolerance for xDHFTI( ) rank test.
!
      TAU = TOL*ANORM
!
!     Compute Householder orthogonal decomposition of matrix.
!
      CALL xDCOPY (N, 0.D0, 0, WS, 1)
      CALL xDCOPY (MA, W(1, NP1), 1, WS, 1)
      K = MAX(M,N)
      MINMAN = MIN(MA,N)
      N1 = K + 1
      N2 = N1 + N
!      CALL xDHFTI (W, MDW, MA, N, WS, MA, 1, TAU, KRANK, RNORM, WS(N2),              &
!     &           WS(N1), IP)
        ! KARLINE:ADDED THAT ...
        MDB = MAX(MA,N)
      CALL xDHFTI (W, MDW, MA, N, WS, MDB, 1, TAU, KRANK, RNORM, WS(N2),              &
     &           WS(N1), IP)   ! and changed that...
      FAC = 1.D0
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
!
!     Reduce to DLPDP and solve.
!
      MAP1 = MA + 1
!
!     Compute inequality rt-hand side for DLPDP.
!
      IF (MA.LT.M) THEN
         IF (MINMAN.GT.0) THEN
            DO 120 I = MAP1,M
               W(I,NP1) = W(I,NP1) - xDDOT(N,W(I,1),MDW,WS,1)
  120       CONTINUE
!
!           Apply permutations to col. of inequality constraint matrix.
!
            DO 130 I = 1,MINMAN
               CALL xDSWAP (MG, W(MAP1,I), 1, W(MAP1,IP(I)), 1)
  130       CONTINUE
!
!           Apply Householder transformations to constraint matrix.
!
            IF (KRANK.GT.0 .AND. KRANK.LT.N) THEN
               DO 140 I = KRANK,1,-1
                  CALL xDH12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),             &
     &                      W(MAP1,1), MDW, 1, MG)
  140          CONTINUE
            ENDIF
!
!           Compute permuted inequality constraint matrix times r-inv.
!
            DO 160 I = MAP1,M
               DO 150 J = 1,KRANK
                W(I,J) = (W(I,J)-xDDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  150          CONTINUE
  160       CONTINUE
         ENDIF
!
!        Solve the reduced problem with DLPDP algorithm,
!        the least projected distance problem.
!
         CALL DLPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X,                  &
     &             XNORM, MDLPDP, WS(N2), IP(N+1))
!
!        Compute solution in original coordinates.
!
         IF (MDLPDP.EQ.1) THEN
            DO 170 I = KRANK,1,-1
               X(I) = (X(I)-xDDOT(KRANK-I,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170       CONTINUE
!
!           Apply Householder transformation to solution vector.
!
            IF (KRANK.LT.N) THEN
               DO 180 I = 1,KRANK
                  CALL xDH12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),             &
     &                      X, 1, 1, 1)
  180          CONTINUE
            ENDIF
!
!           Repermute variables to their input order.
!
            IF (MINMAN.GT.0) THEN
               DO 190 I = MINMAN,1,-1
                  CALL xDSWAP (1, X(I), 1, X(IP(I)), 1)
  190          CONTINUE
!
!              Variables are now in original coordinates.
!              Add solution of unconstrained problem.
!
               DO 200 I = 1,N
                  X(I) = X(I) + WS(I)
  200          CONTINUE
!
!              Compute the residual vector norm.
!
               RNORM = SQRT(RNORM**2+XNORM**2)
            ENDIF
         ELSE
            MODE = 2
         ENDIF
      ELSE
         CALL xDCOPY (N, WS, 1, X, 1)
      ENDIF
!
!     Compute covariance matrix based on the orthogonal decomposition
!     from xDHFTI( ).
!
      IF (.NOT.COV .OR. KRANK.LE.0) GO TO 370
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
!
!     Copy diagonal terms to working array.
!
      CALL xDCOPY (KRANK, W, MDW+1, WS(N2), 1)
!
!     Reciprocate diagonal terms.
!
      DO 210 J = 1,KRANK
         W(J,J) = 1.D0/W(J,J)
  210 CONTINUE
!
!     Invert the upper triangular QR factor on itself.
!
      IF (KRANK.GT.1) THEN
         DO 230 I = 1,KRM1
            DO 220 J = I+1,KRANK
               W(I,J) = -xDDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  220       CONTINUE
  230    CONTINUE
      ENDIF
!
!     Compute the inverted factor times its transpose.
!
      DO 250 I = 1,KRANK
         DO 240 J = I,KRANK
            W(I,J) = xDDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  240    CONTINUE
  250 CONTINUE
!
!     Zero out lower trapezoidal part.
!     Copy upper triangular to lower triangular part.
!
      IF (KRANK.LT.N) THEN
         DO 260 J = 1,KRANK
            CALL xDCOPY (J, W(1,J), 1, W(J,1), MDW)
  260    CONTINUE
!
         DO 270 I = KRP1,N
            CALL xDCOPY (I, 0.D0, 0, W(I,1), MDW)
  270    CONTINUE
!
!        Apply right side transformations to lower triangle.
!
         N3 = N2 + KRP1
         DO 330 I = 1,KRANK
            L = N1 + I
            K = N2 + I
            RB = WS(L-1)*WS(K-1)
!
!           If RB.GE.0.D0, transformation can be regarded as zero.
!
            IF (RB.LT.0.D0) THEN
               RB = 1.D0/RB
!
!              Store unscaled rank one Householder update in work array.
!
               CALL xDCOPY (N, 0.D0, 0, WS(N3), 1)
               L = N1 + I
               K = N3 + I
               WS(K-1) = WS(L-1)
!
               DO 280 J = KRP1,N
                  WS(N3+J-1) = W(I,J)
  280          CONTINUE
!
               DO 290 J = 1,N
                  WS(J) = RB*(xDDOT(J-I,W(J,I),MDW,WS(N3+I-1),1)+                    &
     &                    xDDOT(N-J+1,W(J,J),1,WS(N3+J-1),1))
  290          CONTINUE
!
               L = N3 + I
               GAM = 0.5D0*RB*xDDOT(N-I+1,WS(L-1),1,WS(I),1)
               CALL xDAXPY (N-I+1, GAM, WS(L-1), 1, WS(I), 1)
               DO 320 J = I,N
                  DO 300 L = 1,I-1
                     W(J,L) = W(J,L) + WS(N3+J-1)*WS(L)
  300             CONTINUE
!
                  DO 310 L = I,J
                     W(J,L) = W(J,L) + WS(J)*WS(N3+L-1)+WS(L)*WS(N3+J-1)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  330    CONTINUE
!
!        Copy lower triangle to upper triangle to symmetrize the
!        covariance matrix.
!
         DO 340 I = 1,N
            CALL xDCOPY (I, W(I,1), MDW, W(1,I), 1)
  340    CONTINUE
      ENDIF
!
!     Repermute rows and columns.
!
      DO 350 I = MINMAN,1,-1
         K = IP(I)
         IF (I.NE.K) THEN
            CALL xDSWAP (1, W(I,I), 1, W(K,K), 1)
            CALL xDSWAP (I-1, W(1,I), 1, W(1,K), 1)
            CALL xDSWAP (K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
            CALL xDSWAP (N-K, W(I, K+1), MDW, W(K, K+1), MDW)
         ENDIF
  350 CONTINUE
!
!     Put in normalized residual sum of squares scale factor
!     and symmetrize the resulting covariance matrix.
!
      DO 360 J = 1,N
         CALL xDSCAL (J, FAC, W(1,J), 1)
         CALL xDCOPY (J, W(1,J), 1, W(J,1), MDW)
  360 CONTINUE
!
  370 IP(1) = KRANK
      IP(2) = N + MAX(M,N) + (MG+2)*(N+7)
      RETURN
      END





      SUBROUTINE DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE,                  &
     &   RNORM, IDOPE, DOPE, DONE)
!***BEGIN PROLOGUE  DWNLIT
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DWNNLS
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to DWNNLS( ).
!     The documentation for DWNNLS( ) has complete usage instructions.
!
!     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
!           B as the (N+1)st col.
!
!     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
!     col interchanges.
!
!***SEE ALSO  DWNNLS
!***ROUTINES CALLED  xDCOPY, xDH12, xDROTM, xDROTMG, xDSCAL, xDSWAP, DWNLT1,
!                    DWNLT2, DWNLT3, xIDAMAX
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900604  DP version created from SP version. .  (RWC)
!***END PROLOGUE  DWNLIT
      INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
      DOUBLE PRECISION DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
      LOGICAL DONE
!
      EXTERNAL xDCOPY, xDH12, xDROTM, xDROTMG, xDSCAL, xDSWAP, DWNLT1,              &
     &   DWNLT2, DWNLT3, xIDAMAX
      INTEGER xIDAMAX
      LOGICAL DWNLT2
!
      DOUBLE PRECISION ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5),             &
     &   T, TAU
      INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME,              &
     &   MEND, NIV, NSOLN
      LOGICAL INDEP, RECALC
!
!***FIRST EXECUTABLE STATEMENT  DWNLIT
      ME    = IDOPE(1)
      NSOLN = IDOPE(2)
      L1    = IDOPE(3)
!
      ALSQ   = DOPE(1)
      EANORM = DOPE(2)
      TAU    = DOPE(3)
!
      LB     = MIN(M-1,L)
      RECALC = .TRUE.
      RNORM  = 0.D0
      KRANK  = 0
!
!     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
!     included in the test for column independence.
!
      FACTOR = 1.D0
      LEND = L
      DO 180 I=1,LB
!
!        Set IR to point to the I-th row.
!
         IR = I
         MEND = M
         CALL DWNLT1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,            &
     &                W)
!
!        Update column SS and find pivot column.
!
         CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!        Perform column interchange.
!        Test independence of incoming column.
!
  130    IF (DWNLT2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
!
!           Eliminate I-th column below diagonal using modified Givens
!           transformations applied to (A B).
!
!           When operating near the ME line, use the largest element
!           above it as the pivot.
!
            DO 160 J=M,I+1,-1
               JP = J-1
               IF (J.EQ.ME+1) THEN
                  IMAX = ME
                  AMAX = SCALE(ME)*W(ME,I)**2
                  DO 150 JP=J-1,I,-1
                     T = SCALE(JP)*W(JP,I)**2
                     IF (T.GT.AMAX) THEN
                        IMAX = JP
                        AMAX = T
                     ENDIF
  150             CONTINUE
                  JP = IMAX
               ENDIF
!
               IF (W(J,I).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(JP), SCALE(J), W(JP,I), W(J,I),                &
     &                         SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW,                 &
     &                        SPARAM)
               ENDIF
  160       CONTINUE
         ELSE IF (LEND.GT.I) THEN
!
!           Column I is dependent.  Swap with column LEND.
!           Perform column interchange,
!           and find column in remaining set with largest SS.
!
            CALL DWNLT3 (I, LEND, M, MDW, IPIVOT, H, W)
            LEND = LEND - 1
            IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
            HBAR = H(IMAX)
            GO TO 130
         ELSE
            KRANK = I - 1
            GO TO 190
         ENDIF
  180 CONTINUE
      KRANK = L1
!
  190 IF (KRANK.LT.ME) THEN
         FACTOR = ALSQ
         DO 200 I=KRANK+1,ME
            CALL xDCOPY (L, 0.D0, 0, W(I,1), MDW)
  200    CONTINUE
!
!        Determine the rank of the remaining equality constraint
!        equations by eliminating within the block of constrained
!        variables.  Remove any redundant constraints.
!
         RECALC = .TRUE.
         LB = MIN(L+ME-KRANK, N)
         DO 270 I=L+1,LB
            IR = KRANK + I - L
            LEND = N
            MEND = ME
            CALL DWNLT1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H,               &
     &                   SCALE, W)
!
!           Update col ss and find pivot col
!
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!           Perform column interchange
!           Eliminate elements in the I-th col.
!
            DO 240 J=ME,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL xDROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),               &
     &                        SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW,                 &
     &                        SPARAM)
               ENDIF
  240       CONTINUE
!
!           I=column being eliminated.
!           Test independence of incoming column.
!           Remove any redundant or dependent equality constraints.
!
            IF (.NOT.DWNLT2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
               JJ = IR
               DO 260 IR=JJ,ME
                  CALL xDCOPY (N, 0.D0, 0, W(IR,1), MDW)
                  RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
                  W(IR,N+1) = 0.D0
                  SCALE(IR) = 1.D0
!
!                 Reclassify the zeroed row as a least squares equation.
!
                  ITYPE(IR) = 1
  260          CONTINUE
!
!              Reduce ME to reflect any discovered dependent equality
!              constraints.
!
               ME = JJ - 1
               GO TO 280
            ENDIF
  270    CONTINUE
      ENDIF
!
!     Try to determine the variables KRANK+1 through L1 from the
!     least squares equations.  Continue the triangularization with
!     pivot element W(ME+1,I).
!
  280 IF (KRANK.LT.L1) THEN
         RECALC = .TRUE.
!
!        Set FACTOR=ALSQ to remove effect of heavy weight from
!        test for column independence.
!
         FACTOR = ALSQ
         DO 350 I=KRANK+1,L1
!
!           Set IR to point to the ME+1-st row.
!
            IR = ME+1
            LEND = L
            MEND = M
            CALL DWNLT1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,            &
     &                   W)
!
!           Update column SS and find pivot column.
!
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!           Perform column interchange.
!           Eliminate I-th column below the IR-th element.
!
            DO 320 J=M,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL xDROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),               &
     &                        SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW,               &
     &                        SPARAM)
               ENDIF
  320       CONTINUE
!
!           Test if new pivot element is near zero.
!           If so, the column is dependent.
!           Then check row norm test to be classified as independent.
!
            T = SCALE(IR)*W(IR,I)**2
            INDEP = T .GT. (TAU*EANORM)**2
            IF (INDEP) THEN
               RN = 0.D0
               DO 340 I1=IR,M
                  DO 330 J1=I+1,N
                     RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
               INDEP = T .GT. RN*TAU**2
            ENDIF
!
!           If independent, swap the IR-th and KRANK+1-th rows to
!           maintain the triangular form.  Update the rank indicator
!           KRANK and the equality constraint pointer ME.
!
            IF (.NOT.INDEP) GO TO 360
            CALL xDSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
            CALL xDSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
!
!           Reclassify the least square equation as an equality
!           constraint and rescale it.
!
            ITYPE(IR) = 0
            T = SQRT(SCALE(KRANK+1))
            CALL xDSCAL(N+1, T, W(KRANK+1,1), MDW)
            SCALE(KRANK+1) = ALSQ
            ME = ME+1
            KRANK = KRANK+1
  350    CONTINUE
      ENDIF
!
!     If pseudorank is less than L, apply Householder transformation.
!     from right.
!
  360 IF (KRANK.LT.L) THEN
         DO 370 J=KRANK,1,-1
            CALL xDH12 (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1,              &
     &                J-1)
  370    CONTINUE
      ENDIF
!
      NIV = KRANK + NSOLN - L
      IF (L.EQ.N) DONE = .TRUE.
!
!     End of initial triangularization.
!
      IDOPE(1) = ME
      IDOPE(2) = KRANK
      IDOPE(3) = NIV
      RETURN
      END





      SUBROUTINE DWNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,             &
     &   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
!***BEGIN PROLOGUE  DWNLSM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DWNNLS
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to DWNNLS.
!     The documentation for DWNNLS has complete usage instructions.
!
!     In addition to the parameters discussed in the prologue to
!     subroutine DWNNLS, the following work arrays are used in
!     subroutine DWNLSM  (they are passed through the calling
!     sequence from DWNNLS for purposes of variable dimensioning).
!     Their contents will in general be of no interest to the user.
!
!     Variables of type REAL are DOUBLE PRECISION.
!
!         IPIVOT(*)
!            An array of length N.  Upon completion it contains the
!         pivoting information for the cols of W(*,*).
!
!         ITYPE(*)
!            An array of length M which is used to keep track
!         of the classification of the equations.  ITYPE(I)=0
!         denotes equation I as an equality constraint.
!         ITYPE(I)=1 denotes equation I as a least squares
!         equation.
!
!         WD(*)
!            An array of length N.  Upon completion it contains the
!         dual solution vector.
!
!         H(*)
!            An array of length N.  Upon completion it contains the
!         pivot scalars of the Householder transformations performed
!         in the case KRANK.LT.L.
!
!         SCALE(*)
!            An array of length M which is used by the subroutine
!         to store the diagonal matrix of weights.
!         These are used to apply the modified Givens
!         transformations.
!
!         Z(*),TEMP(*)
!            Working arrays of length N.
!
!         D(*)
!            An array of length N that contains the
!         column scaling for the matrix (E).
!                                       (A)
!
!***SEE ALSO  DWNNLS
!***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDH12, xDNRM2, xDROTM,
!                    xDROTMG, xDSCAL, xDSWAP, DWNLIT, xIDAMAX, xXERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Fixed an error message.  (RWC)
!   900604  DP version created from SP version.  (RWC)
!   900911  Restriction on value of ALAMDA included.  (WRB)
!***END PROLOGUE  DWNLSM
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      DOUBLE PRECISION D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),             &
     &   W(MDW,*), WD(*), X(*), Z(*)
!
      EXTERNAL D1MACH,xDASUM,xDAXPY,xDCOPY,xDH12,xDNRM2,xDROTM,xDROTMG,             &
     &   xDSCAL, xDSWAP, DWNLIT, xIDAMAX, xXERMSG
      DOUBLE PRECISION D1MACH, xDASUM, xDNRM2
      INTEGER xIDAMAX
!
      DOUBLE PRECISION ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,                    &
     &   DOPE(3), DRELPR, EANORM, FAC, SM, SPARAM(5), T, TAU, WMAX, Z2,             &
     &   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,                &
     &   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,             &  
     &   NOPT, NSOLN, NTIMES
      LOGICAL DONE, FEASBL, FIRST, HITCON, POS
!
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DWNLSM
!
!     Initialize variables.
!     DRELPR is the precision for the particular machine
!     being used.  This logic avoids resetting it every entry.
!
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
!
!     Set the nominal tolerance used in the code.
!
      TAU = SQRT(DRELPR)
!
      M = MA + MME
      ME = MME
      MODE = 2
!
!     To process option vector
!
      FAC = 1.D-4
!
!     Set the nominal blow up factor used in the code.
!
      BLOWUP = TAU
!
!     The nominal column scaling used in the code is
!     the identity scaling.
!
      CALL xDCOPY (N, 1.D0, 0, D, 1)
!
!     Define bound for number of options to change.
!
      NOPT = 1000
!
!     Define bound for positive value of LINK.
!
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL xXERMSG ('SLATEC', 'DWNLSM',                                          &
     &      'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
!
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL xXERMSG ('SLATEC', 'DWNLSM',                                          &
     &      'IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.',               &
     &      3, 1)
            RETURN
         ENDIF
!
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = xDNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               D(J) = T
  110       CONTINUE
         ENDIF
!
         IF (KEY.EQ.7) CALL xDCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(DRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(DRELPR,PRGOPT(LAST+2))
!
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL xXERMSG ('SLATEC', 'DWNLSM',                                       &
     &         'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
!
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
!
      DO 120 J = 1,N
         CALL xDSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
!
!     Process option vector
!
      DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
!
!     Compute scale factor to apply to equality constraint equations.
!
      DO 130 J = 1,N
         WD(J) = xDASUM(M,W(1,J),1)
  130 CONTINUE
!
      IMAX = xIDAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = xDASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(DRELPR*FAC)
!
!     On machines, such as the VAXes using D floating, with a very
!     limited exponent range for double precision values, the previously
!     computed value of ALAMDA may cause an overflow condition.
!     Therefore, this code further limits the value of ALAMDA.
!
      ALAMDA = MIN(ALAMDA,SQRT(D1MACH(2)))
!
!     Define scaling diagonal matrix for modified Givens usage and
!     classify equation types.
!
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
!
!        When equation I is heavily weighted ITYPE(I)=0,
!        else ITYPE(I)=1.
!
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.D0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
!
!     Set the solution vector X(*) to zero and the column interchange
!     matrix to the identity.
!
      CALL xDCOPY (N, 0.D0, 0, X, 1)
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
!
!     Perform initial triangularization in the submatrix
!     corresponding to the unconstrained variables.
!     Set first L components of dual vector to zero because
!     these correspond to the unconstrained variables.
!
      CALL xDCOPY (L, 0.D0, 0, WD, 1)
!
!     The arrays IDOPE(*) and DOPE(*) are used to pass
!     information to DWNLIT().  This was done to avoid
!     a long calling sequence or the use of COMMON.
!
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
!
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,                 &
     &            IDOPE, DOPE, DONE)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
!
!     Perform WNNLS algorithm using the following steps.
!
!     Until(DONE)
!        compute search direction and feasible point
!        when (HITCON) add constraints
!        else perform multiplier test and drop a constraint
!        fin
!     Compute-Final-Solution
!
!     To compute search direction and feasible point,
!     solve the triangular system of currently non-active
!     variables and store the solution in Z(*).
!
!     To solve system
!     Copy right hand side into TEMP vector to use overwriting method.
!
  160 IF (DONE) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL xDCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
!
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL xDAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
!
!     Increment iteration counter and check against maximum number
!     of iterations.
!
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         DONE = .TRUE.
      ENDIF
!
!     Check to see if any constraints have become active.
!     If so, calculate an interpolation factor so that all
!     active constraints are removed from the basis.
!
      ALPHA = 2.D0
      HITCON = .FALSE.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.D0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .TRUE.
         ENDIF
  180 CONTINUE
!
!     Compute search direction and feasible point
!
      IF (HITCON) THEN
!
!        To add constraints, use computed ALPHA to interpolate between
!        last feasible solution X(*) and current unconstrained (and
!        infeasible) solution Z(*).
!
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .FALSE.
!
!        Remove column JCON and shift columns JCON+1 through N to the
!        left.  Swap column JCON into the N th position.  This achieves
!        upper Hessenberg form for the nonactive constraints and
!        leaves an upper Hessenberg matrix to retriangularize.
!
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL xDCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
!
!        Update permuted index vector to reflect this shift and swap.
!
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
!
!        Similarly permute X(*) vector.
!
         CALL xDCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.D0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
!
!        Retriangularize upper Hessenberg matrix after adding
!        constraints.
!
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
!
!              Zero IP1 to I in column J
!
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
!
!              Zero IP1 to I in column J
!
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL xDSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL xDSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
!
!              Swapped row was formerly a pivot element, so it will
!              be large enough to perform elimination.
!              Zero IP1 to I in column J.
!
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
!
!                 Zero IP1 to I in column J
!
                  IF (W(I+1,J).NE.0.D0) THEN
                     CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J),                     &
     &                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.D0
                     CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,             &
     &                           SPARAM)
                  ENDIF
               ELSE
                  CALL xDSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL xDSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.D0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
!
!        See if the remaining coefficients in the solution set are
!        feasible.  They should be because of the way ALPHA was
!        determined.  If any are infeasible, it is due to roundoff
!        error.  Any that are non-positive will be set to zero and
!        removed from the solution set.
!
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.D0) GO TO 250
  240    CONTINUE
         FEASBL = .TRUE.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
!
!        To perform multiplier test and drop a constraint.
!
         CALL xDCOPY (NSOLN, Z, 1, X, 1)
         IF (NSOLN.LT.N) CALL xDCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
!
!        Reclassify least squares equations as equalities as necessary.
!
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL xDSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL xDSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
!
!        Form inner product vector WD(*) of dual coefficients.
!
         DO 280 J = NSOLN+1,N
            SM = 0.D0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
!
!        Find J such that WD(J)=WMAX is maximum.  This determines
!        that the incoming column J will reduce the residual vector
!        and be positive.
!
  290    WMAX = 0.D0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.D0) GO TO 330
!
!        Set dual coefficients to zero for incoming column.
!
         WD(IWMAX) = 0.D0
!
!        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set.
!        Perform transformation to retriangularize, and test for near
!        linear dependence.
!
!        Swap column IWMAX into NSOLN-th position to maintain upper
!        Hessenberg form of adjacent columns, and add new column to
!        triangular decomposition.
!
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL xDSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.D0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
!
!        Reduce column NSOLN so that the matrix of nonactive constraints
!        variables is triangular.
!
         DO 320 J = M,NIV+1,-1
            JP = J - 1
!
!           When operating near the ME line, test to see if the pivot
!           element is near zero.  If so, use the largest element above
!           it as the pivot.  This is to maintain the sharp interface
!           between weighted and non-weighted rows in all cases.
!
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
!
            IF (W(J,NSOLN).NE.0.D0) THEN
               CALL xDROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN),                       &
     &                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.D0
               CALL xDROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),             &
     &                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
!
!        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
!        this is nonpositive or too large.  If this was true or if the
!        pivot term was zero, reject the column as dependent.
!
         IF (W(NIV,NSOLN).NE.0.D0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.D0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
!
!           Try to add row ME+1 as an additional equality constraint.
!           Check size of proposed new solution component.
!           Reject it if it is too large.
!
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.D0) THEN
            ISOL = ME + 1
            IF (POS) THEN
!
!              Swap rows ME+1 and NIV, and scale factors for these rows.
!
               CALL xDSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL xDSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .FALSE.
         ENDIF
!
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.DONE)) GO TO 290
      ENDIF
      GO TO 160
!
!     Else perform multiplier test and drop a constraint.  To compute
!     final solution.  Solve system, store results in X(*).
!
!     Copy right hand side into TEMP vector to use overwriting method.
!
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL xDCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
!
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL xDAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
!
!     Solve system.
!
      CALL xDCOPY (NSOLN, Z, 1, X, 1)
!
!     Apply Householder transformations to X(*) if KRANK.LT.L
!
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL xDH12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
!
!     Fill in trailing zeroes for constrained variables not in solution.
!
      IF (NSOLN.LT.N) CALL xDCOPY (N-NSOLN, 0.D0, 0, X(NSOLN+1), 1)
!
!     Permute solution vector to natural order.
!
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
!
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL xDSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
!
!     Rescale the solution using the column scaling.
!
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
!
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
!
      RNORM = SQRT(RNORM)
      RETURN
      END





      SUBROUTINE DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,             &
     &   SCALE, W)
!***BEGIN PROLOGUE  DWNLT1
!***SUBSIDIARY
!***PURPOSE  Subsidiary to WNLIT
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     To update the column Sum Of Squares and find the pivot column.
!     The column Sum of Squares Vector will be updated at each step.
!     When numerically necessary, these values will be recomputed.
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  xIDAMAX
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT1
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
!
      EXTERNAL xIDAMAX
      INTEGER xIDAMAX
!
      INTEGER J, K
!
!***FIRST EXECUTABLE STATEMENT  DWNLT1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
!
!        Update column SS=sum of squares.
!
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
!
!        Test for numerical accuracy.
!
         IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
!
!     If required, recalculate column SS, using rows IR through MEND.
!
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.D0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
!
!        Find column with largest SS.
!
         IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END





      LOGICAL FUNCTION DWNLT2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
!***BEGIN PROLOGUE  DWNLT2
!***SUBSIDIARY
!***PURPOSE  Subsidiary to WNLIT
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     To test independence of incoming column.
!
!     Test the column IC to determine if it is linearly independent
!     of the columns already in the basis.  In the initial tri. step,
!     we usually want the heavy weight ALAMDA to be included in the
!     test for independence.  In this case, the value of FACTOR will
!     have been set to 1.E0 before this procedure is invoked.
!     In the potentially rank deficient problem, the value of FACTOR
!     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
!     heavy weight from the test for independence.
!
!     Write new column as partitioned vector
!           (A1)  number of components in solution so far = NIV
!           (A2)  M-NIV components
!     And compute  SN = inverse weighted length of A1
!                  RN = inverse weighted length of A2
!     Call the column independent when RN .GT. TAU*SN
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT2
      DOUBLE PRECISION FACTOR, SCALE(*), TAU, WIC(*)
      INTEGER IR, ME, MEND
!
      DOUBLE PRECISION RN, SN, T
      INTEGER J
!
!***FIRST EXECUTABLE STATEMENT  DWNLT2
      SN = 0.E0
      RN = 0.E0
      DO 10 J=1,MEND
         T = SCALE(J)
         IF (J.LE.ME) T = T/FACTOR
         T = T*WIC(J)**2
!
         IF (J.LT.IR) THEN
            SN = SN + T
         ELSE
            RN = RN + T
         ENDIF
   10 CONTINUE
      DWNLT2 = RN .GT. SN*TAU**2
      RETURN
      END





      SUBROUTINE DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!***BEGIN PROLOGUE  DWNLT3
!***SUBSIDIARY
!***PURPOSE  Subsidiary to WNLIT
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Perform column interchange.
!     Exchange elements of permuted index vector and perform column
!     interchanges.
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  xDSWAP
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      DOUBLE PRECISION H(*), W(MDW,*)
!
      EXTERNAL xDSWAP
!
      DOUBLE PRECISION T
      INTEGER ITEMP
!
!***FIRST EXECUTABLE STATEMENT  DWNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
!
         CALL xDSWAP(M, W(1,IMAX), 1, W(1,I), 1)
!
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END





      SUBROUTINE DWNNLS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE,              &
     &   IWORK, WORK)
!***BEGIN PROLOGUE  DWNNLS
!***PURPOSE  Solve a linearly constrained least squares problem with
!            equality constraints and nonnegativity constraints on
!            selected variables.
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A
!***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D)
!***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
!             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
!             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Abstract
!
!     This subprogram solves a linearly constrained least squares
!     problem.  Suppose there are given matrices E and A of
!     respective dimensions ME by N and MA by N, and vectors F
!     and B of respective lengths ME and MA.  This subroutine
!     solves the problem
!
!               EX = F, (equations to be exactly satisfied)
!
!               AX = B, (equations to be approximately satisfied,
!                        in the least squares sense)
!
!               subject to components L+1,...,N nonnegative
!
!     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
!
!     The problem is reposed as problem DWNNLS
!
!               (WT*E)X = (WT*F)
!               (   A)    (   B), (least squares)
!               subject to components L+1,...,N nonnegative.
!
!     The subprogram chooses the heavy weight (or penalty parameter) WT.
!
!     The parameters for DWNNLS are
!
!     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     W(*,*),MDW,  The array W(*,*) is double subscripted with first
!     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
!                  discussion let us call M = ME + MA.  Then MDW
!                  must satisfy MDW.GE.M.  The condition MDW.LT.M
!                  is an error.
!
!                  The array W(*,*) contains the matrices and vectors
!
!                       (E  F)
!                       (A  B)
!
!                  in rows and columns 1,...,M and 1,...,N+1
!                  respectively.  Columns 1,...,L correspond to
!                  unconstrained variables X(1),...,X(L).  The
!                  remaining variables are constrained to be
!                  nonnegative. The condition L.LT.0 or L.GT.N is
!                  an error.
!
!     PRGOPT(*)    This double precision array is the option vector.
!                  If the user is satisfied with the nominal
!                  subprogram features set
!
!                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!                  Otherwise PRGOPT(*) is a linked list consisting of
!                  groups of data of the following form
!
!                  LINK
!                  KEY
!                  DATA SET
!
!                  The parameters LINK and KEY are each one word.
!                  The DATA SET can be comprised of several words.
!                  The number of items depends on the value of KEY.
!                  The value of LINK points to the first
!                  entry of the next group of data within
!                  PRGOPT(*).  The exception is when there are
!                  no more options to change.  In that
!                  case LINK=1 and the values KEY and DATA SET
!                  are not referenced. The general layout of
!                  PRGOPT(*) is as follows.
!
!               ...PRGOPT(1)=LINK1 (link to first entry of next group)
!               .  PRGOPT(2)=KEY1 (key to the option change)
!               .  PRGOPT(3)=DATA VALUE (data value for this change)
!               .       .
!               .       .
!               .       .
!               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
!               .                       next group)
!               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
!               .  PRGOPT(LINK1+2)=DATA VALUE
!               ...     .
!               .       .
!               .       .
!               ...PRGOPT(LINK)=1 (no more options to change)
!
!                  Values of LINK that are nonpositive are errors.
!                  A value of LINK.GT.NLINK=100000 is also an error.
!                  This helps prevent using invalid but positive
!                  values of LINK that will probably extend
!                  beyond the program limits of PRGOPT(*).
!                  Unrecognized values of KEY are ignored.  The
!                  order of the options is arbitrary and any number
!                  of options can be changed with the following
!                  restriction.  To prevent cycling in the
!                  processing of the option array a count of the
!                  number of options changed is maintained.
!                  Whenever this count exceeds NOPT=1000 an error
!                  message is printed and the subprogram returns.
!
!                  OPTIONS..
!
!                  KEY=6
!                         Scale the nonzero columns of the
!                  entire data matrix
!                  (E)
!                  (A)
!                  to have length one. The DATA SET for
!                  this option is a single value.  It must
!                  be nonzero if unit length column scaling is
!                  desired.
!
!                  KEY=7
!                         Scale columns of the entire data matrix
!                  (E)
!                  (A)
!                  with a user-provided diagonal matrix.
!                  The DATA SET for this option consists
!                  of the N diagonal scaling factors, one for
!                  each matrix column.
!
!                  KEY=8
!                         Change the rank determination tolerance from
!                  the nominal value of SQRT(SRELPR).  This quantity
!                  can be no smaller than SRELPR, The arithmetic-
!                  storage precision.  The quantity used
!                  here is internally restricted to be at
!                  least SRELPR.  The DATA SET for this option
!                  is the new tolerance.
!
!                  KEY=9
!                         Change the blow-up parameter from the
!                  nominal value of SQRT(SRELPR).  The reciprocal of
!                  this parameter is used in rejecting solution
!                  components as too large when a variable is
!                  first brought into the active set.  Too large
!                  means that the proposed component times the
!                  reciprocal of the parameter is not less than
!                  the ratio of the norms of the right-side
!                  vector and the data matrix.
!                  This parameter can be no smaller than SRELPR,
!                  the arithmetic-storage precision.
!
!                  For example, suppose we want to provide
!                  a diagonal matrix to scale the problem
!                  matrix and change the tolerance used for
!                  determining linear dependence of dropped col
!                  vectors.  For these options the dimensions of
!                  PRGOPT(*) must be at least N+6.  The FORTRAN
!                  statements defining these options would
!                  be as follows.
!
!                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
!                  PRGOPT(2)=7 (user-provided scaling key)
!
!                  CALL xDCOPY(N,D,1,PRGOPT(3),1) (copy the N
!                  scaling factors from a user array called D(*)
!                  into PRGOPT(3)-PRGOPT(N+2))
!
!                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
!                  PRGOPT(N+4)=8 (linear dependence tolerance key)
!                  PRGOPT(N+5)=... (new value of the tolerance)
!
!                  PRGOPT(N+6)=1 (no more options to change)
!
!
!     IWORK(1),    The amounts of working storage actually allocated
!     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
!                  respectively.  These quantities are compared with
!                  the actual amounts of storage needed for DWNNLS( ).
!                  Insufficient storage allocated for either WORK(*)
!                  or IWORK(*) is considered an error.  This feature
!                  was included in DWNNLS( ) because miscalculating
!                  the storage formulas for WORK(*) and IWORK(*)
!                  might very well lead to subtle and hard-to-find
!                  execution errors.
!
!                  The length of WORK(*) must be at least
!
!                  LW = ME+MA+5*N
!                  This test will not be made if IWORK(1).LE.0.
!
!                  The length of IWORK(*) must be at least
!
!                  LIW = ME+MA+N
!                  This test will not be made if IWORK(2).LE.0.
!
!     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     X(*)         An array dimensioned at least N, which will
!                  contain the N components of the solution vector
!                  on output.
!
!     RNORM        The residual norm of the solution.  The value of
!                  RNORM contains the residual vector length of the
!                  equality constraints and least squares equations.
!
!     MODE         The value of MODE indicates the success or failure
!                  of the subprogram.
!
!                  MODE = 0  Subprogram completed successfully.
!
!                       = 1  Max. number of iterations (equal to
!                            3*(N-L)) exceeded. Nearly all problems
!                            should complete in fewer than this
!                            number of iterations. An approximate
!                            solution and its corresponding residual
!                            vector length are in X(*) and RNORM.
!
!                       = 2  Usage error occurred.  The offending
!                            condition is noted with the error
!                            processing subprogram, xXERMSG( ).
!
!     User-designated
!     Working arrays..
!
!     WORK(*)      A double precision working array of length at least
!                  M + 5*N.
!
!     IWORK(*)     An integer-valued working array of length at least
!                  M+N.
!
!***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Report SAND77-0552, Sandia
!                 Laboratories, June 1978.
!               K. H. Haskell and R. J. Hanson, Selected algorithms for
!                 the linearly constrained least squares problem - a
!                 users guide, Report SAND78-1290, Sandia Laboratories,
!                 August 1979.
!               K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Mathematical Programming
!                 21 (1981), pp. 98-118.
!               R. J. Hanson and K. H. Haskell, Two algorithms for the
!                 linearly constrained least squares problem, ACM
!                 Transactions on Mathematical Software, September 1982.
!               !. L. Lawson and R. J. Hanson, Solving Least Squares
!                 Problems, Prentice-Hall, Inc., 1974.
!***ROUTINES CALLED  DWNLSM, xXERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
!   900510  Convert XERRWV calls to xXERMSG calls, change Prologue
!           comments to agree with WNNLS.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DWNNLS
      INTEGER IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW, ME,                &
     &     MODE, N
      DOUBLE PRECISION  PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
      CHARACTER*8 XERN1
!***FIRST EXECUTABLE STATEMENT  DWNNLS
      MODE = 0
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
!
      IF (IWORK(1).GT.0) THEN
         LW = ME + MA + 5*N
         IF (IWORK(1).LT.LW) THEN
            WRITE (XERN1, '(I8)') LW
            CALL xXERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //            &
     &         'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
!
      IF (IWORK(2).GT.0) THEN
         LIW = ME + MA + N
         IF (IWORK(2).LT.LIW) THEN
            WRITE (XERN1, '(I8)') LIW
            CALL xXERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //            &
     &         'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
!
      IF (MDW.LT.ME+MA) THEN
         CALL xXERMSG ('SLATEC', 'DWNNLS',                                          &
     &      'THE VALUE MDW.LT.ME+MA IS AN ERROR', 2, 1)
         MODE = 2
         RETURN
      ENDIF
!
      IF (L.LT.0 .OR. L.GT.N) THEN
         CALL xXERMSG ('SLATEC', 'DWNNLS',                                          &
     &      'L.GE.0 .AND. L.LE.N IS REQUIRED', 2, 1)
         MODE = 2
         RETURN
      ENDIF
!
!     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
!     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
!     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
!
      L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
!
      CALL DWNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,              &
     &            IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),                 &
     &            WORK(L4), WORK(L5))
      RETURN
      END



      SUBROUTINE xDROTM (N,DX,INCX,DY,INCY,DPARAM)
!
!     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
!
!     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
!     (DY**T)
!
!     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
!     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
!     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!
!     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!
!       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!     H=(          )    (          )    (          )    (          )
!       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
!
      DOUBLE PRECISION DFLAG,DH12,DH22,DX,TWO,Z,DH11,DH21,                          &
     & DPARAM,DY,W,ZERO
      INTEGER N,INCX,INCY,NSTEPS,I,KX,KY
      DIMENSION DX(1),DY(1),DPARAM(5)
      DATA ZERO,TWO/0.D0,2.D0/
!
      DFLAG=DPARAM(1)
      IF(N .LE. 0 .OR.(DFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
!
               NSTEPS=N*INCX
               IF(DFLAG) 50,10,30
   10          CONTINUE
               DH12=DPARAM(4)
               DH21=DPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W+Z*DH12
                    DY(I)=W*DH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               DH11=DPARAM(2)
               DH22=DPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z
                    DY(I)=-W+DH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               DH11=DPARAM(2)
               DH12=DPARAM(4)
               DH21=DPARAM(3)
               DH22=DPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z*DH12
                    DY(I)=W*DH21+Z*DH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
!
          IF(DFLAG)120,80,100
   80     CONTINUE
          DH12=DPARAM(4)
          DH21=DPARAM(3)
               DO 90 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W+Z*DH12
               DY(KY)=W*DH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          DH11=DPARAM(2)
          DH22=DPARAM(5)
               DO 110 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z
               DY(KY)=-W+DH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          DH11=DPARAM(2)
          DH12=DPARAM(4)
          DH21=DPARAM(3)
          DH22=DPARAM(5)
               DO 130 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z*DH12
               DY(KY)=W*DH21+Z*DH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
          END





      SUBROUTINE xDROTMG (DD1,DD2,DX1,DY1,DPARAM)
!
!     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
!     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
!     DY2)**T.
!     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!
!     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!
!       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!     H=(          )    (          )    (          )    (          )
!       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
!     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
!     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
!
!     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
!     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
!     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
!
      DOUBLE PRECISION GAM,ONE,RGAMSQ,DD2,DH11,DH21,DPARAM,DP2,                     &
     & DQ2,DU,DY1,ZERO,GAMSQ,DD1,DFLAG,DH12,DH22,DP1,DQ1,                           &
     & DTEMP,DX1,TWO
      INTEGER IGO
      DIMENSION DPARAM(5)
!
      DATA ZERO,ONE,TWO /0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
      IF(.NOT. DD1 .LT. ZERO) GO TO 10
!       GO ZERO-H-D-AND-DX1..
          GO TO 60
   10 CONTINUE
!     CASE-DD1-NONNEGATIVE
      DP2=DD2*DY1
      IF(.NOT. DP2 .EQ. ZERO) GO TO 20
          DFLAG=-TWO
          GO TO 260
!     REGULAR-CASE..
   20 CONTINUE
      DP1=DD1*DX1
      DQ2=DP2*DY1
      DQ1=DP1*DX1
!
      IF(.NOT. DABS(DQ1) .GT. DABS(DQ2)) GO TO 40
          DH21=-DY1/DX1
          DH12=DP2/DP1
!
          DU=ONE-DH12*DH21
!
          IF(.NOT. DU .LE. ZERO) GO TO 30
!         GO ZERO-H-D-AND-DX1..
               GO TO 60
   30     CONTINUE
               DFLAG=ZERO
               DD1=DD1/DU
               DD2=DD2/DU
               DX1=DX1*DU
!         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. DQ2 .LT. ZERO) GO TO 50
!         GO ZERO-H-D-AND-DX1..
               GO TO 60
   50     CONTINUE
               DFLAG=ONE
               DH11=DP1/DP2
               DH22=DX1/DY1
               DU=ONE+DH11*DH22
               DTEMP=DD2/DU
               DD2=DD1/DU
               DD1=DTEMP
               DX1=DY1*DU
!         GO SCALE-CHECK
               GO TO 100
!     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
          DFLAG=-ONE
          DH11=ZERO
          DH12=ZERO
          DH21=ZERO
          DH22=ZERO
!
          DD1=ZERO
          DD2=ZERO
          DX1=ZERO
!         RETURN..
          GO TO 220
!     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. DFLAG .GE. ZERO) GO TO 90
!
          IF(.NOT. DFLAG .EQ. ZERO) GO TO 80
          DH11=ONE
          DH22=ONE
          DFLAG=-ONE
          GO TO 90
   80     CONTINUE
          DH21=-ONE
          DH12=ONE
          DFLAG=-ONE
   90 CONTINUE
C      GO TO IGO,(120,150,180,210)  CHANGED INTO
      SELECT CASE (IGO)

        CASE (120)
          GOTO 120
        CASE (150)
          GOTO 150
        CASE (180)
          GOTO 180
        CASE (210)
          GOTO 210
                   
      END SELECT
      

C karline: IGO DOES NOT HAVE A VALUE ?
!     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. DD1 .LE. RGAMSQ) GO TO 130
               IF(DD1 .EQ. ZERO) GO TO 160
C               ASSIGN 120 TO IGO  CHANGED INTO:
               IGO = 120
!              FIX-H..
               GO TO 70
  120          CONTINUE
               DD1=DD1*GAM**2
               DX1=DX1/GAM
               DH11=DH11/GAM
               DH12=DH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. DD1 .GE. GAMSQ) GO TO 160
C               ASSIGN 150 TO IGO   CAHNGED INTO
               IGO = 150
!              FIX-H..
               GO TO 70
  150          CONTINUE
               DD1=DD1/GAM**2
               DX1=DX1*GAM
               DH11=DH11*GAM
               DH12=DH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. DABS(DD2) .LE. RGAMSQ) GO TO 190
               IF(DD2 .EQ. ZERO) GO TO 220
C               ASSIGN 180 TO IGO  CHANGED INTO
               IGO = 180
!              FIX-H..
               GO TO 70
  180          CONTINUE
               DD2=DD2*GAM**2
               DH21=DH21/GAM
               DH22=DH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. DABS(DD2) .GE. GAMSQ) GO TO 220
C               ASSIGN 210 TO IGO  CHANGED INTO
               IGO = 210

!              FIX-H..
               GO TO 70
  210          CONTINUE
               DD2=DD2/GAM**2
               DH21=DH21*GAM
               DH22=DH22*GAM
          GO TO 200
  220 CONTINUE
          IF(DFLAG)250,230,240
  230     CONTINUE
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               GO TO 260
  240     CONTINUE
               DPARAM(2)=DH11
               DPARAM(5)=DH22
               GO TO 260
  250     CONTINUE
               DPARAM(2)=DH11
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               DPARAM(5)=DH22
  260 CONTINUE
          DPARAM(1)=DFLAG
          RETURN
      END



      SUBROUTINE DLPDP (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS,              &
     &   IS)
!***BEGIN PROLOGUE  DLPDP
!***SUBSIDIARY
!***PURPOSE  Subsidiary to xDLSEI
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (LPDP-S, DLPDP-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!  **** Double Precision version of LPDP ****
!     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
!     where N=N1+N2.  This is a slight overestimate for WS(*).
!
!     Determine an N1-vector W, and
!               an N2-vector Z
!     which minimizes the Euclidean length of W
!     subject to G*W+H*Z .GE. Y.
!     This is the least projected distance problem, LPDP.
!     The matrices G and H are of respective
!     dimensions M by N1 and M by N2.
!
!     Called by subprogram DLSI( ).
!
!     The matrix
!                (G H Y)
!
!     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
!
!     The solution (W) is returned in X(*).
!                  (Z)
!
!     The value of MODE indicates the status of
!     the computation after returning to the user.
!
!          MODE=1  The solution was successfully obtained.
!
!          MODE=2  The inequalities are inconsistent.
!
!***SEE ALSO  xDLSEI
!***ROUTINES CALLED  xDCOPY, xDDOT, xDNRM2, xDSCAL, DWNNLS
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR section.  (WRB)
!***END PROLOGUE  DLPDP
!
      INTEGER I, IS(*), IW, IX, J, L, M, MDA, MODE, MODEW, N, N1, N2,               &
     &     NP1
      DOUBLE PRECISION A(MDA,*), xDDOT, xDNRM2, FAC, ONE,                             &
     &     PRGOPT(*), RNORM, SC, WNORM, WS(*), X(*), YNORM, ZERO
      SAVE ZERO, ONE, FAC
      DATA ZERO,ONE /0.0D0,1.0D0/, FAC /0.1D0/
!***FIRST EXECUTABLE STATEMENT  DLPDP
      N = N1 + N2
      MODE = 1
      IF (M .GT. 0) GO TO 20
         IF (N .LE. 0) GO TO 10
            X(1) = ZERO
            CALL xDCOPY(N,X,0,X,1)
   10    CONTINUE
         WNORM = ZERO
      GO TO 200
   20 CONTINUE
!        BEGIN BLOCK PERMITTING ...EXITS TO 190
            NP1 = N + 1
!
!           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
            DO 40 I = 1, M
               SC = xDNRM2(N,A(I,1),MDA)
               IF (SC .EQ. ZERO) GO TO 30
                  SC = ONE/SC
                  CALL xDSCAL(NP1,SC,A(I,1),MDA)
   30          CONTINUE
   40       CONTINUE
!
!           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
            YNORM = xDNRM2(M,A(1,NP1),1)
            IF (YNORM .EQ. ZERO) GO TO 50
               SC = ONE/YNORM
               CALL xDSCAL(M,SC,A(1,NP1),1)
   50       CONTINUE
!
!           SCALE COLS OF MATRIX H.
            J = N1 + 1
   60       IF (J .GT. N) GO TO 70
               SC = xDNRM2(M,A(1,J),1)
               IF (SC .NE. ZERO) SC = ONE/SC
               CALL xDSCAL(M,SC,A(1,J),1)
               X(J) = SC
               J = J + 1
            GO TO 60
   70       CONTINUE
            IF (N1 .LE. 0) GO TO 130
!
!              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
               IW = 0
               DO 80 I = 1, M
!
!                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
                  CALL xDCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
!
!                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
                  CALL xDCOPY(N1,A(I,1),MDA,WS(IW+1),1)
                  IW = IW + N1
!
!                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
   80          CONTINUE
               WS(IW+1) = ZERO
               CALL xDCOPY(N,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N
               WS(IW+1) = ONE
               IW = IW + 1
!
!              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
!              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
!              F = TRANSPOSE OF (0,...,0,1).
               IX = IW + 1
               IW = IW + M
!
!              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
!              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,NP1,N2,NP1-N2,M,0,PRGOPT,WS(IX),RNORM,                &
     &                     MODEW,IS,WS(IW+1))
!
!              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
               SC = ONE - xDDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)                 &
     &            GO TO 110
                  SC = ONE/SC
                  DO 90 J = 1, N1
                     X(J) = SC*xDDOT(M,A(1,J),1,WS(IX),1)
   90             CONTINUE
!
!                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
!                 VECTOR.
                  DO 100 I = 1, M
                     A(I,NP1) = A(I,NP1) - xDDOT(N1,A(I,1),MDA,X,1)
  100             CONTINUE
               GO TO 120
  110          CONTINUE
                  MODE = 2
!        .........EXIT
                  GO TO 190
  120          CONTINUE
  130       CONTINUE
            IF (N2 .LE. 0) GO TO 180
!
!              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
               IW = 0
               DO 140 I = 1, M
                  CALL xDCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
  140          CONTINUE
               WS(IW+1) = ZERO
               CALL xDCOPY(N2,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N2
               WS(IW+1) = ONE
               IW = IW + 1
               IX = IW + 1
               IW = IW + M
!
!              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
!              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
!              OF (0,...,0,1)).
!
!              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
!              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,N2+1,0,N2+1,M,0,PRGOPT,WS(IX),RNORM,MODEW,            &
     &                     IS,WS(IW+1))
!
!              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
               SC = ONE - xDDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)                 &
     &            GO TO 160
                  SC = ONE/SC
                  DO 150 J = 1, N2
                     L = N1 + J
                     X(L) = SC*xDDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150             CONTINUE
               GO TO 170
  160          CONTINUE
                  MODE = 2
!        .........EXIT
                  GO TO 190
  170          CONTINUE
  180       CONTINUE
!
!           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
            CALL xDSCAL(N,YNORM,X,1)
            WNORM = xDNRM2(N1,X,1)
  190    CONTINUE
  200 CONTINUE
      RETURN
      END




      SUBROUTINE xDH12 (MODE, LPIVOT, L1, M, U, IUE, UP, C, ICE, ICV,                &
     &   NCV)
!***BEGIN PROLOGUE  DH12
!***SUBSIDIARY
!***PURPOSE  Subsidiary to xDHFTI, xDLSEI and DWNNLS
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (H12-S, DH12-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!      &** DOUBLE PRECISION VERSION OF H12 ******
!
!     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
!     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
!
!     Construction and/or application of a single
!     Householder transformation..     Q = I + U*(U**T)/B
!
!     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
!     LPIVOT is the index of the pivot element.
!     L1,M   If L1 .LE. M   the transformation will be constructed to
!            zero elements indexed from L1 through M.   If L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
!                   IUE is the storage increment between elements.
!                                       On exit from H1 U() and UP
!                   contain quantities defining the vector U of the
!                   Householder transformation.   On entry to H2 U()
!                   and UP should contain quantities previously computed
!                   by H1.  These will not be modified by H2.
!     C()    On entry to H1 or H2 C() contains a matrix which will be
!            regarded as a set of vectors to which the Householder
!            transformation is to be applied.  On exit C() contains the
!            set of transformed vectors.
!     ICE    Storage increment between elements of vectors in C().
!     ICV    Storage increment between vectors in C().
!     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
!            no operations will be done on C().
!
!***SEE ALSO  xDHFTI, xDLSEI, DWNNLS
!***ROUTINES CALLED  xDAXPY, xDDOT, xDSWAP
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900911  Added xDDOT to DOUBLE PRECISION statement.  (WRB)
!***END PROLOGUE  DH12
      INTEGER I, I2, I3, I4, ICE, ICV, INCR, IUE, J, KL1, KL2, KLP,                 &
     &     L1, L1M1, LPIVOT, M, MML1P2, MODE, NCV
      DOUBLE PRECISION B, C, CL, CLINV, ONE, UL1M1, SM, U, UP, xDDOT
      DIMENSION U(IUE,*), C(*)
!     BEGIN BLOCK PERMITTING ...EXITS TO 140
!***FIRST EXECUTABLE STATEMENT  DH12
         ONE = 1.0D0
!
!     ...EXIT
         IF (0 .GE. LPIVOT .OR. LPIVOT .GE. L1 .OR. L1 .GT. M) GO TO 140
         CL = ABS(U(1,LPIVOT))
         IF (MODE .EQ. 2) GO TO 40
!           &***** CONSTRUCT THE TRANSFORMATION. ******
            DO 10 J = L1, M
               CL = MAX(ABS(U(1,J)),CL)
   10       CONTINUE
            IF (CL .GT. 0.0D0) GO TO 20
!     .........EXIT
               GO TO 140
   20       CONTINUE
            CLINV = ONE/CL
            SM = (U(1,LPIVOT)*CLINV)**2
            DO 30 J = L1, M
               SM = SM + (U(1,J)*CLINV)**2
   30       CONTINUE
            CL = CL*SQRT(SM)
            IF (U(1,LPIVOT) .GT. 0.0D0) CL = -CL
            UP = U(1,LPIVOT) - CL
            U(1,LPIVOT) = CL
         GO TO 50
   40    CONTINUE
!        &***** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
!
         IF (CL .GT. 0.0D0) GO TO 50
!     ......EXIT
            GO TO 140
   50    CONTINUE
!     ...EXIT
         IF (NCV .LE. 0) GO TO 140
         B = UP*U(1,LPIVOT)
!        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
!
         IF (B .LT. 0.0D0) GO TO 60
!     ......EXIT
            GO TO 140
   60    CONTINUE
         B = ONE/B
         MML1P2 = M - L1 + 2
         IF (MML1P2 .LE. 20) GO TO 80
            L1M1 = L1 - 1
            KL1 = 1 + (L1M1 - 1)*ICE
            KL2 = KL1
            KLP = 1 + (LPIVOT - 1)*ICE
            UL1M1 = U(1,L1M1)
            U(1,L1M1) = UP
            IF (LPIVOT .NE. L1M1) CALL xDSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
            DO 70 J = 1, NCV
               SM = xDDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
               SM = SM*B
               CALL xDAXPY(MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
               KL1 = KL1 + ICV
   70       CONTINUE
            U(1,L1M1) = UL1M1
!     ......EXIT
            IF (LPIVOT .EQ. L1M1) GO TO 140
            KL1 = KL2
            CALL xDSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
         GO TO 130
   80    CONTINUE
            I2 = 1 - ICV + ICE*(LPIVOT - 1)
            INCR = ICE*(L1 - LPIVOT)
            DO 120 J = 1, NCV
               I2 = I2 + ICV
               I3 = I2 + INCR
               I4 = I3
               SM = C(I2)*UP
               DO 90 I = L1, M
                  SM = SM + C(I3)*U(1,I)
                  I3 = I3 + ICE
   90          CONTINUE
               IF (SM .EQ. 0.0D0) GO TO 110
                  SM = SM*B
                  C(I2) = C(I2) + SM*UP
                  DO 100 I = L1, M
                     C(I4) = C(I4) + SM*U(1,I)
                     I4 = I4 + ICE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      RETURN
      END




      SUBROUTINE xDHFTI (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H,             &
     &   G, IP)
!***BEGIN PROLOGUE  xDHFTI
!***PURPOSE  Solve a least squares problem for banded matrices using
!            sequential accumulation of rows of the data matrix.
!            Exactly one right-hand side vector is permitted.
!***LIBRARY   SLATEC
!***CATEGORY  D9
!***TYPE      DOUBLE PRECISION (HFTI-S, xDHFTI-D)
!***KEYWORDS  CURVE FITTING, LEAST SQUARES
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
!
!     This subroutine solves a linear least squares problem or a set of
!     linear least squares problems having the same matrix but different
!     right-side vectors.  The problem data consists of an M by N matrix
!     A, an M by NB matrix B, and an absolute tolerance parameter TAU
!     whose usage is described below.  The NB column vectors of B
!     represent right-side vectors for NB distinct linear least squares
!     problems.
!
!     This set of problems can also be written as the matrix least
!     squares problem
!
!                       AX = B,
!
!     where X is the N by NB solution matrix.
!
!     Note that if B is the M by M identity matrix, then X will be the
!     pseudo-inverse of A.
!
!     This subroutine first transforms the augmented matrix (A B) to a
!     matrix (R C) using premultiplying Householder transformations with
!     column interchanges.  All subdiagonal elements in the matrix R are
!     zero and its diagonal elements satisfy
!
!                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)),
!
!                       I = 1,...,L-1, where
!
!                       L = MIN(M,N).
!
!     The subroutine will compute an integer, KRANK, equal to the number
!     of diagonal terms of R that exceed TAU in magnitude. Then a
!     solution of minimum Euclidean length is computed using the first
!     KRANK rows of (R C).
!
!     To be specific we suggest that the user consider an easily
!     computable matrix norm, such as, the maximum of all column sums of
!     magnitudes.
!
!     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
!     norm of B), it is suggested that TAU be set approximately equal to
!     EPS*(norm of A).
!
!     The user must dimension all arrays appearing in the call list..
!     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
!     permits the solution of a range of problems in the same array
!     space.
!
!     The entire set of parameters for xDHFTI are
!
!     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
!                       matrix A of the least squares problem AX = B.
!                       The first dimensioning parameter of the array
!                       A(*,*) is MDA, which must satisfy MDA.GE.M
!                       Either M.GE.N or M.LT.N is permitted.  There
!                       is no restriction on the rank of A.  The
!                       condition MDA.LT.M is considered an error.
!
!     B(*),MDB,NB       If NB = 0 the subroutine will perform the
!                       orthogonal decomposition but will make no
!                       references to the array B(*).  If NB.GT.0
!                       the array B(*) must initially contain the M by
!                       NB matrix B of the least squares problem AX =
!                       B.  If NB.GE.2 the array B(*) must be doubly
!                       subscripted with first dimensioning parameter
!                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
!                       be either doubly or singly subscripted.  In
!                       the latter case the value of MDB is arbitrary
!                       but it should be set to some valid integer
!                       value such as MDB = M.
!
!                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
!                       is considered an error.
!
!     TAU               Absolute tolerance parameter provided by user
!                       for pseudorank determination.
!
!     H(*),G(*),IP(*)   Arrays of working space used by xDHFTI.
!
!     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     A(*,*)            The contents of the array A(*,*) will be
!                       modified by the subroutine. These contents
!                       are not generally required by the user.
!
!     B(*)              On return the array B(*) will contain the N by
!                       NB solution matrix X.
!
!     KRANK             Set by the subroutine to indicate the
!                       pseudorank of A.
!
!     RNORM(*)          On return, RNORM(J) will contain the Euclidean
!                       norm of the residual vector for the problem
!                       defined by the J-th column vector of the array
!                       B(*,*) for J = 1,...,NB.
!
!     H(*),G(*)         On return these arrays respectively contain
!                       elements of the pre- and post-multiplying
!                       Householder transformations used to compute
!                       the minimum Euclidean length solution.
!
!     IP(*)             Array in which the subroutine records indices
!                       describing the permutation of column vectors.
!                       The contents of arrays H(*),G(*) and IP(*)
!                       are not generally required by the user.
!
!***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
!                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
!***ROUTINES CALLED  D1MACH, xDH12, xXERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
!   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  xDHFTI
      INTEGER I, II, IOPT, IP(*), IP1, J, JB, JJ, K, KP1, KRANK, L,                 &
     &     LDIAG, LMAX, M, MDA, MDB, N, NB, NERR
      DOUBLE PRECISION A, B, D1MACH, DZERO, FACTOR,                                 &
     &     G, H, HMAX, RELEPS, RNORM, SM, SM1, SZERO, TAU, TMP
      DIMENSION A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*)
      SAVE RELEPS
      DATA RELEPS /0.D0/
!     BEGIN BLOCK PERMITTING ...EXITS TO 360
!***FIRST EXECUTABLE STATEMENT  xDHFTI
         IF (RELEPS.EQ.0.D0) RELEPS = D1MACH(4)
         SZERO = 0.0D0
         DZERO = 0.0D0
         FACTOR = 0.001D0
!
         K = 0
         LDIAG = MIN(M,N)
         IF (LDIAG .LE. 0) GO TO 350
!           BEGIN BLOCK PERMITTING ...EXITS TO 130
!              BEGIN BLOCK PERMITTING ...EXITS TO 120
                  IF (MDA .GE. M) GO TO 10
                     NERR = 1
                     IOPT = 2
                     CALL xXERMSG ('SLATEC', 'xDHFTI',                               &
     &                  'MDA.LT.M, PROBABLE ERROR.',                                &
     &                  NERR, IOPT)
!     ...............EXIT
                     GO TO 360
   10             CONTINUE
!
                  IF (NB .LE. 1 .OR. MAX(M,N) .LE. MDB) GO TO 20
                     NERR = 2
                     IOPT = 2
                     CALL xXERMSG ('SLATEC', 'xDHFTI',                               &
     &                  'MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.',             &
     &                  NERR, IOPT)
!     ...............EXIT
                     GO TO 360
   20             CONTINUE
!
                  DO 100 J = 1, LDIAG
!                    BEGIN BLOCK PERMITTING ...EXITS TO 70
                        IF (J .EQ. 1) GO TO 40
!
!                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
!                          ..
                           LMAX = J
                           DO 30 L = J, N
                              H(L) = H(L) - A(J-1,L)**2
                              IF (H(L) .GT. H(LMAX)) LMAX = L
   30                      CONTINUE
!                    ......EXIT
                           IF (FACTOR*H(LMAX) .GT. HMAX*RELEPS) GO TO 70
   40                   CONTINUE
!
!                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
!                       ..
                        LMAX = J
                        DO 60 L = J, N
                           H(L) = 0.0D0
                           DO 50 I = J, M
                              H(L) = H(L) + A(I,L)**2
   50                      CONTINUE
                           IF (H(L) .GT. H(LMAX)) LMAX = L
   60                   CONTINUE
                        HMAX = H(LMAX)
   70                CONTINUE
!                    ..
!                     LMAX HAS BEEN DETERMINED
!
!                     DO COLUMN INTERCHANGES IF NEEDED.
!                    ..
                     IP(J) = LMAX
                     IF (IP(J) .EQ. J) GO TO 90
                        DO 80 I = 1, M
                           TMP = A(I,J)
                           A(I,J) = A(I,LMAX)
                           A(I,LMAX) = TMP
   80                   CONTINUE
                        H(LMAX) = H(J)
   90                CONTINUE
!
!                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
!                     AND B.
!                    ..
                     CALL xDH12(1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,              &
     &                         N-J)
                     CALL xDH12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
  100             CONTINUE
!
!                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
!                  TAU.
!                 ..
                  DO 110 J = 1, LDIAG
!              ......EXIT
                     IF (ABS(A(J,J)) .LE. TAU) GO TO 120
  110             CONTINUE
                  K = LDIAG
!           ......EXIT
                  GO TO 130
  120          CONTINUE
               K = J - 1
  130       CONTINUE
            KP1 = K + 1
!
!           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
!
            IF (NB .LT. 1) GO TO 170
            DO 160 JB = 1, NB
               TMP = SZERO
               IF (M .LT. KP1) GO TO 150
               DO 140 I = KP1, M
                  TMP = TMP + B(I,JB)**2
  140          CONTINUE
  150          CONTINUE
               RNORM(JB) = SQRT(TMP)
  160       CONTINUE
  170       CONTINUE
!           SPECIAL FOR PSEUDORANK = 0
            IF (K .GT. 0) GO TO 210
               IF (NB .LT. 1) GO TO 200
               DO 190 JB = 1, NB
                  DO 180 I = 1, N
                     B(I,JB) = SZERO
  180             CONTINUE
  190          CONTINUE
  200          CONTINUE
            GO TO 340
  210       CONTINUE
!
!               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
!               DECOMPOSITION OF FIRST K ROWS.
!              ..
               IF (K .EQ. N) GO TO 230
                  DO 220 II = 1, K
                     I = KP1 - II
                     CALL xDH12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  220             CONTINUE
  230          CONTINUE
!
!
               IF (NB .LT. 1) GO TO 330
               DO 320 JB = 1, NB
!
!                  SOLVE THE K BY K TRIANGULAR SYSTEM.
!                 ..
                  DO 260 L = 1, K
                     SM = DZERO
                     I = KP1 - L
                     IP1 = I + 1
                     IF (K .LT. IP1) GO TO 250
                     DO 240 J = IP1, K
                        SM = SM + A(I,J)*B(J,JB)
  240                CONTINUE
  250                CONTINUE
                     SM1 = SM
                     B(I,JB) = (B(I,JB) - SM1)/A(I,I)
  260             CONTINUE
!
!                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
!                 ..
                  IF (K .EQ. N) GO TO 290
                     DO 270 J = KP1, N
                        B(J,JB) = SZERO
  270                CONTINUE
                     DO 280 I = 1, K
                        CALL xDH12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,              &
     &                            MDB,1)
  280                CONTINUE
  290             CONTINUE
!
!                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
!                   COLUMN INTERCHANGES.
!                 ..
                  DO 310 JJ = 1, LDIAG
                     J = LDIAG + 1 - JJ
                     IF (IP(J) .EQ. J) GO TO 300
                        L = IP(J)
                        TMP = B(L,JB)
                        B(L,JB) = B(J,JB)
                        B(J,JB) = TMP
  300                CONTINUE
  310             CONTINUE
  320          CONTINUE
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
!        ..
!         THE SOLUTION VECTORS, X, ARE NOW
!         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
!
         KRANK = K
  360 CONTINUE
      RETURN
      END







!********************************************************************
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

!-------------------------------------------------------------------*
! Subroutine obtained from LINPACK (ftp-site netlib.att.com)        *
!-------------------------------------------------------------------*

      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
!
!     xdgbsl solves the double precision band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgbco or xdgbfa.
!
!     on entry
!
!        abd     double precision(lda, n)
!                the output from dgbco or xdgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from dgbco or xdgbfa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgbco has set rcond .gt. 0.0
!        or xdgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call xdgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     xblas xdaxpy,xddot
!     fortran min0
!
!     internal variables
!
      double precision xddot,t
      integer k,kb,l,la,lb,lm,m,nm1
!
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call xdaxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call xdaxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = xddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + xddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end                                        ! of dgbsl


                       
!********************************************************************
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)

!-------------------------------------------------------------------*
! Subroutine obtained from LINPACK (ftp-site netlib.att.com)        *
!-------------------------------------------------------------------*

      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)
!
!     xdgbfa factors a double precision band matrix by elimination.
!
!     xdgbfa is usually called by dgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that xdgbsl will divide by zero if
!                     called.  use  rcond  in dgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas xdaxpy,xdscal,xidamax
!     fortran max0,min0
!
!     internal variables
!
      double precision t
      integer i,xidamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
!
!        zero next fill-in column
!
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
!
!        find l = pivot index
!
         lm = min0(ml,n-k)
         l = xidamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
         if (abd(l,k) .eq. 0.0d0) go to 100
!
!           interchange if necessary
!
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
!
!           compute multipliers
!
            t = -1.0d0/abd(m,k)
            call xdscal(lm,t,abd(m+1,k),1)
!
!           row elimination with column indexing
!
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call xdaxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end                                        ! of dgbfa


!********************************************************************
      subroutine dgesl(a,MXLDA,n,ipvt,b,job)

!-------------------------------------------------------------------*
! Subroutine obtained from LINPACK (ftp-site netlib.att.com)        *
!-------------------------------------------------------------------*

      integer n,ipvt(1),job,MxLDA
!      double precision a(lda,1),b(1)
      double precision a(MXlda,1),b(1)
!
!     xdgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or xdgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or xdgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     ipvt    integer(n)
!                the pivot vector from dgeco or xdgefa.
!
!     b       double precision(n)
!                the right hand side vector.
!
!     job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or xdgefa has set info .eq. 0 .
! 
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call xdgesl(a,lda,n,ipvt,c(1,j),0) 
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     xblas, xdaxpy,xddot
!
!     internal variables
!
      double precision xddot,t
      integer k,kb,l,nm1
!
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call xdaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call xdaxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            t = xddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + xddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end                                        ! of xdgesl


!********************************************************************
      subroutine dgefa(a,MXLDA,n,ipvt,info)

!-------------------------------------------------------------------*
! Subroutine obtained from LINPACK (ftp-site netlib.att.com)        *
!-------------------------------------------------------------------*

      integer n,ipvt(1),info,MXLDA
      double precision a(MXlda,1)
!
!     xdgefa factors a double precision matrix by gaussian elimination.
!
!     xdgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for xdgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that xdgesl or xdgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas xdaxpy,xdscal,xidamax
!
!     internal variables
!
      double precision t
      integer xidamax,j,k,kp1,l,nm1
!
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = xidamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (a(l,k) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           compute multipliers
!
            t = -1.0d0/a(k,k)
            call xdscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call xdaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end                                        ! of xdgefa

*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
*DECK DSCAL
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
      INTEGER N,INCX,INCY,IX,IY,I,M,MP1,NS
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
      INTEGER N,INCX,INCY,IX,IY,I,M,MP1,NS
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
      INTEGER N,INCX,INCY,IX,IY,I,M,MP1,NS

C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END








      double precision function Xdasum(n,dx,incx)
C
C     takes the sum of the absolute values.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      IMPLICIT NONE 
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
C
      Xdasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
C
C        code for increment not equal to 1
C
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      Xdasum = dtemp
      return
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))  &
     &  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 Xdasum = dtemp
      return
      end




      SUBROUTINE xXERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)

C***PURPOSE  Process error messages for SLATEC and other libraries.

      IMPLICIT NONE
      CHARACTER(LEN=*)   LIBRAR, SUBROU, MESSG
      CHARACTER(LEN=8)   XLIBR, XSUBR
      CHARACTER(LEN=72)  TEMP
      CHARACTER(LEN=20)  LFIRST
      INTEGER :: LKNTRL,MAXMES,NERR,Lerr,LEVEL,LLEVEL,LTEMP,Kount
      INTEGER :: MKNTRL,I

CKARLINE:
       IF (Level .GE. 2) THEN 
C         WRITE(*,*) " Unrecoverable error C- program stopped"
       CALL XMESSAGE (" Unrecoverable error C- LINPACK routine stopped")
       RETURN
       ENDIF


C***FIRST EXECUTABLE STATEMENT  xXERMSG
C KARLINE: CHANGED -ERRMSG NOT WRITTEN TO FILE
C      LKNTRL = J4SAVE (2, 0, .FALSE.)
C      MAXMES = J4SAVE (4, 0, .FALSE.)
       LKNTRL = 1
       MAXMES = 1 
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING xXERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.     &
     &   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL xXERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //             &
     &      'xXERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//               &
     &      'JOB ABORT DUE TO FATAL ERROR.', 72)
C KS        CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)

C KS        CALL XERHLT (' ***xXERMSG -- INVALID INPUT')
         
C         WRITE(*,*) ' ***xXERMSG -- INVALID INPUT'
         CALL XMESSAGE                                                    &
     & (" ***xXERMSG -- INVALID INPUT- LINPACK routine stopped")
         RETURN
C
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
C KS      I = J4SAVE (1, NERR, .TRUE.)
C KS      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C KS: 
      KOUNT = 1 
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
C KS      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XxERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.                         &
     &       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL xXERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL xXERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL xXERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
C KS         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL xXERPRN (' *  ', -1, ' ', 72)
         CALL xXERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL xXERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL xXERPRN                                                  &
     &         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
         CALL xXERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
C KS      CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
C         CALL XERHLT (' ')
         STOP
C
      ELSE
C         CALL XERHLT (MESSG)
C         WRITE(*,*) MESSG
          CALL XMESSAGE (MESSG)

         STOP
C
      ENDIF
      RETURN
      END SUBROUTINE xXERMSG

C*********************************************************************

      SUBROUTINE XXERPRN (PREFIX, NPREF, MESSG, NWRAP)

CC***PURPOSE  Print error messages processed by xXERMSG.
C KARLINE: DO NOT PRINT TO SCREEN ->          CALL XMESSAGE 

      IMPLICIT NONE
      CHARACTER(LEN=*)   PREFIX, MESSG
      INTEGER            NPREF, NWRAP
      CHARACTER(LEN=148) CBUFF
      CHARACTER(LEN=2)   NEWLIN
      PARAMETER (NEWLIN = '$$')
      INTEGER :: LPREF,LWRAP,LENMSG,NEXTC,LPIECE,IDELTA,I,N

C KARLINE: RETURN
C      RETURN

C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
C            WRITE(*, '(A)') CBUFF(1:LPREF+1)
         CALL XMESSAGE (CBUFF(1:LPREF+1))

         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
C         WRITE(*, '(A)') CBUFF(1:LPREF+LPIECE)
         CALL XMESSAGE (CBUFF(1:LPREF+LPIECE))
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN

      END SUBROUTINE XXERPRN



