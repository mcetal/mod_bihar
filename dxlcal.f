*DECK DXLCAL
      SUBROUTINE DXLCAL (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
     +   WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A,
     +   ISYM)
C***BEGIN PROLOGUE  DXLCAL
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SXLCAL-S, DXLCAL-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This  routine computes the solution  XL,  the current DGMRES
C        iterate, given the  V(I)'s and  the  QR factorization of the
C        Hessenberg  matrix HES.   This routine  is  only called when
C        ITOL=11.
C
C *Usage:
C      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED)
C      INTEGER NELT, IA(NELT), JA(NELT), ISYM
C      DOUBLE PRECISION X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL),
C     $                 V(N,MAXLP1), R0NRM, WK(N), SZ(N),
C     $                 RPAR(USER DEFINED), A(NELT)
C      EXTERNAL MSOLVE
C
C      CALL DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
C     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
C     $     NELT, IA, JA, A, ISYM)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C LGMR   :IN       Integer
C         The number of iterations performed and
C         the current order of the upper Hessenberg
C         matrix HES.
C X      :IN       Double Precision X(N)
C         The current approximate solution as of the last restart.
C XL     :OUT      Double Precision XL(N)
C         An array of length N used to hold the approximate
C         solution X(L).
C         Warning: XL and ZL are the same array in the calling routine.
C ZL     :IN       Double Precision ZL(N)
C         An array of length N used to hold the approximate
C         solution Z(L).
C HES    :IN       Double Precision HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,i) and V(*,k).
C MAXLP1 :IN       Integer
C         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C         MAXL is the maximum allowable order of the matrix HES.
C Q      :IN       Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR.
C V      :IN       Double Precision V(N,MAXLP1)
C         The N by(LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C R0NRM  :IN       Double Precision
C         The scaled norm of the initial residual for the
C         current call to DPIGMR.
C WK     :IN       Double Precision WK(N)
C         A double precision work array of length N.
C SZ     :IN       Double Precision SZ(N)
C         A vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C JPRE   :IN       Integer
C         The preconditioner type flag.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as below.  RPAR is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IPAR is an integer work array
C         for the same purpose as RPAR.
C NMSL   :IN       Integer
C         The number of calls to MSOLVE.
C RPAR   :IN       Double Precision RPAR(USER DEFINED)
C         Double Precision workspace passed directly to the MSOLVE
C         routine.
C IPAR   :IN       Integer IPAR(USER DEFINED)
C         Integer workspace passed directly to the MSOLVE routine.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Double Precision A(NELT)
C         A double precision array of length NELT containing matrix
C         data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all non-zero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DXLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION R0NRM
      INTEGER ISYM, JPRE, JSCAL, LGMR, MAXLP1, N, NELT, NMSL
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), HES(MAXLP1,*), Q(*), RPAR(*), SZ(*),
     +                 V(N,*), WK(N), X(N), XL(N), ZL(N)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
C     .. Subroutine Arguments ..
      EXTERNAL MSOLVE
C     .. Local Scalars ..
      INTEGER I, K, LL, LLP1
C     .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS
C***FIRST EXECUTABLE STATEMENT  DXLCAL
      LL = LGMR
      LLP1 = LL + 1
      DO 10 K = 1,LLP1
         WK(K) = 0
 10   CONTINUE
      WK(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, WK)
      DO 20 K = 1,N
         ZL(K) = 0
 20   CONTINUE
      DO 30 I = 1,LL
         CALL DAXPY(N, WK(I), V(1,I), 1, ZL, 1)
 30   CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 40 K = 1,N
            ZL(K) = ZL(K)/SZ(K)
 40      CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, ZL, 1, WK, 1)
         CALL MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
C         calculate XL from X and ZL.
      DO 50 K = 1,N
         XL(K) = X(K) + ZL(K)
 50   CONTINUE
      RETURN
C------------- LAST LINE OF DXLCAL FOLLOWS ----------------------------
      END
