!====================================================================
!SUBROUTINE:   dareachar
!PART OF:      sctrob
!DESCRIPTION:  compute the area-specific characteristic: 0=sum; 1=mean
               !works with balanced and unbalanced data 
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER g(n), n(1), typ(1), nsize(g); DOUBLE PRECISION vec
!ON RETURN
!  DOUBLE PRECISION areachar(g)
!--------------------------------------------------------------------
SUBROUTINE dareachar(g, n, nsize, typ, vec, areachar)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: g, n, typ
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: vec(n)
   DOUBLE PRECISION, INTENT(OUT) :: areachar(g)
   !local definitions
   INTEGER :: i, j
   INTEGER :: l(g), u(g)
   DOUBLE PRECISION :: ssw
   !a check
   IF (SUM(nsize) /= n) THEN
      areachar = 0D0
   END IF
   !
   l(1) = 1
   u(1) = nsize(1)
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
      u(i) = u(i - 1) + nsize(i) 
   END DO
   !
   DO i = 1, g
      ssw = 0D0
      DO j = l(i), u(i)
         ssw = ssw + vec(j)
      END DO
      areachar(i) = ssw
   END DO
   !if characteristic == mean
   IF (typ == 1) THEN
      DO i = 1, g
         areachar(i) = areachar(i) / nsize(i)
      END DO
   END IF   
END SUBROUTINE
!====================================================================
!SUBROUTINE:   dcenter 
!PART OF:      sctrob
!DESCRIPTION:  centers a vector by the area-specific mean; works for
!              balanced and unbalanced data
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER g(n), n(1); DOUBLE PRECISION vec
!ON RETURN
!  DOUBLE PRECISION vec(n)
!--------------------------------------------------------------------
SUBROUTINE dcenter(g, n, nsize, vec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: g, n
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(INOUT) :: vec(n)
   !local declaration
   INTEGER :: i, j
   INTEGER :: l(g), u(g)
   DOUBLE PRECISION :: ssm
   !if not it returns nsize as first elements of vec
   IF (SUM(nsize) /= n) THEN
      vec = 0D0
      DO i = 1, g
         vec(i) = nsize(i)
      END DO
   ELSE
      l(1) = 1
      u(1) = nsize(1)
      DO i = 2, g
         l(i) = l(i - 1) + nsize(i - 1)
         u(i) = u(i - 1) + nsize(i) 
      END DO
      !
      DO i = 1, g
         ssm = 0D0
         DO j = l(i), u(i)
            ssm = ssm + (vec(j) / nsize(i))
         END DO
         DO j = l(i), u(i)
            vec(j) = vec(j) - ssm
         END DO
      END DO
   END IF
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   dasymhuberwgt
!PART OF:      sctrob
!DESCRIPTION:  compute huber psi-weight; NOTE: 
!              typ = 1 for sqrt(wgt)  
!              typ = 0 for wgt, 
!              typ = 2 for wgt^2
!              QUANTILE EXTENSION
!BENCHMARK:    robustbase:::huberPsi@wgt (v 0.7-6), 
!              approved, June 19, 2011 
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER n(1), typ(1)
!  REAL k(1), vec(n)
!ON RETURN:
!  REAL vec(n)
!--------------------------------------------------------------------
SUBROUTINE dasymhuberwgt(n, k, typ, vec, q)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, typ
   DOUBLE PRECISION, INTENT(IN) :: k, q
   DOUBLE PRECISION, INTENT(INOUT) :: vec(n)
   !local declarations
   INTEGER :: i
   DOUBLE PRECISION :: choice
   DOUBLE PRECISION :: work(n)
   !
   work = vec
   !
   SELECT CASE (typ) 
      CASE(1) !take the square root of the weights
         DO i = 1, n
            choice = k / ABS(vec(i))
            IF (choice < 1D0) THEN
               vec(i) = SQRT(choice)
            ELSE
               vec(i) = 1
            END IF
            ! introduce the asymmetry
            IF ( work(i) <= 0D0 ) THEN
               vec(i) = vec(i) * SQRT((1 - q) * 2D0)
            ELSE
               vec(i) = vec(i) * SQRT(q * 2D0)
            END IF
         END DO
      CASE(0) !take the weights as they are
         DO i = 1, n
            choice = k / ABS(vec(i))
            IF (choice < 1D0) THEN
               vec(i) = choice
            ELSE
               vec(i) = 1
            END IF
            ! introduce the asymmetry
            IF ( work(i) <= 0D0 ) THEN
               vec(i) = vec(i) * (1 - q) * 2D0
            ELSE
               vec(i) = vec(i) * q * 2D0
            END IF
         END DO
      CASE(2) !the weights to the power of two
         DO i = 1, n
            choice = k / ABS(vec(i))
            IF (choice < 1D0) THEN
               vec(i) = choice ** 2
            ELSE
               vec(i) = 1
            END IF
         END DO
            ! introduce the asymmetry
            IF ( work(i) <= 0D0 ) THEN
               vec(i) = 4D0 * vec(i) * (1 - q) ** 2
            ELSE
               vec(i) = 4D0 * vec(i) * q ** 2     
            END IF
      CASE DEFAULT !an errorneous call returns one 
         vec = 0 
   END SELECT
END SUBROUTINE
!  
!====================================================================
!SUBROUTINE:   drsaebeta
!PART OF:      rsaehuber
!DESCRIPTION:  compute the residual vector
!DEPENDENCY:   
!  dgemv (BLAS2 and LAPACK), dgels (LAPACK) 
!  dhuberwgt, dsqrtinvva
!ON ENTRY:
!  INTEGER n(1), p(1), k(1), nsize(g), info(1), dec(1)
!  REAL k(1), xmat(n,p) yvec(n), d(1), v(1)
!      beta(p), q(1)
!ON RETURN: 
!  INTEGER info(1)
!  REAL beta(p), sumwgt(1)
!--------------------------------------------------------------------
SUBROUTINE dbeta(n, p, g, k, xmat, yvec, v, d, nsize, beta,& 
      sumwgt, info, dec, q)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec, q
   INTEGER, INTENT(IN) :: nsize(g)
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: d, v, k
   DOUBLE PRECISION, INTENT(INOUT) :: beta(p)
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt
   !local declarations (most are used in dqrls)
   INTEGER :: i, j, lwork
   INTEGER, PARAMETER :: lworkmax = 10000
   DOUBLE PRECISION :: work(lworkmax)
   DOUBLE PRECISION :: modyvec(n), res(n) 
   DOUBLE PRECISION :: modxmat(n,p) 
   !
   res = yvec
   CALL dgemv("N", n, p, -1D0, xmat, n, beta, 1, 1D0, res, 1)
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, res)
   CALL dasymhuberwgt(n, k, 1, res, q)
   modxmat = xmat
   modyvec = yvec
   CALL dsqrtinvva(n, p, g, nsize, d, v, 1, dec, modxmat)
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 1, dec, modyvec)
   DO j = 1, p
      sumwgt = 0D0
      DO i = 1, n
         modxmat(i, j) = modxmat(i, j) * res(i)
         modyvec(i) = modyvec(i) * res(i)
         sumwgt = sumwgt + res(i) ** 2 
       END DO  
   END DO
   lwork = -1
   CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
   lwork = MIN(lworkmax, INT(work(1)))
   CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
   IF (info == 0) THEN
      beta = modyvec(1:p)
   ELSE
      beta = 0
   END IF
END SUBROUTINE 

!====================================================================
!SUBROUTINE:   drsaebetaiter
!PART OF:      rsaehuber
!DESCRIPTION:  fully iterated drsaebeta; info carries the # of iter
!DEPENDENCY:   
!  drsaebeta
!  ddelta 
!ON ENTRY:
!  INTEGER n(1), p(1), k(1), nsize(g), iter(1), dec(1)
!  REAL k(1), xmat(n,p) yvec(n), v(1), d(1), acc(1)
!      beta(p), q(1)
!ON RETURN: 
!  INTEGER converged(1), info(1)
!  REAL beta(p), sumwgt(1)
!--------------------------------------------------------------------
SUBROUTINE dbetaiter(n, p, g, k, xmat, yvec, v, d, nsize, acc, &
      beta, iter, converged, sumwgt, info, dec, q)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec, q
   INTEGER, INTENT(IN) :: nsize(g)
   INTEGER, INTENT(IN) :: iter
   INTEGER, INTENT(OUT) :: converged
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: acc, v, d, k
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(INOUT) :: beta(p)
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt
   INTEGER :: i, coinfo, niter
   DOUBLE PRECISION :: betaold(p)
   !
   niter = 0
   DO i = 1, iter
      betaold = beta
      CALL dbeta(n, p, g, k, xmat, yvec, v, d, nsize,&
         beta, sumwgt, coinfo, dec, q)
      IF (coinfo /= 0) THEN
         beta = 0
         EXIT 
      END IF
      niter = niter + 1
      CALL ddelta(p, betaold, beta, acc, converged)
      IF (converged == 1) THEN
         EXIT
      END IF
   END DO
   info = niter
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   ddelta
!PART OF:      sctrob
!DESCRIPTION:  computes the squared norm for two parameter vectors and
!              evaluates if their difference is smaller than the 
!              reference acc; return 1 if true, 0 else.
!              this function computation the termination rule for 
!              iterative algorithms (either for parameters or resids)
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER p(1), res(1)
!  REAL acc(1), oldvec(p), newvec(p)
!ON RETURN
!  INTEGER res(1)
!--------------------------------------------------------------------
SUBROUTINE ddelta(p, oldvec, newvec, acc, res)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: p
   INTEGER, INTENT(OUT) :: res
   DOUBLE PRECISION, INTENT(IN) :: acc
   DOUBLE PRECISION, INTENT(IN) :: oldvec(p), newvec(p)
   !local declaration
   DOUBLE PRECISION :: ratio
   DOUBLE PRECISION, PARAMETER :: themin = 1.0D-15
   !
   ratio = SQRT(SUM((oldvec - newvec)**2) / MAX(SUM(oldvec**2), themin))
   IF (ratio < acc) THEN
      res = 1
   ELSE
      res = 0
   END IF
END SUBROUTINE 
!

!====================================================================
!SUBROUTINE:   dsqrtinvva
!PART OF:      rsaehuber
!STATUS:       June 23, 2011; mod November 3, 2011
!DESCRIPTION:  pre-multiply a matrix (vector) by the square root of
!              the inverse covariance matrix with either
!              decomposition (dec):
!                 0: SVD
!                 1: Choleski
!
!              parametrization
!                 0: MLM
!                 1: Hartley-Rao
!                 (else): returns a zero matrix/vector
!BENCHMARK:    (self-contained testing; approved June 23, 2011; 
!              modifications: Nov 16, 2011)
!DEPENDENCY:   DPOTRF (LAPACK)
!FORTRAN:      uses dynamic allocation (only v90, not v77)
!ON ENTRY:
!  INTEGER n(1), p(1), g(1), nsize(g), par(1), dec(1)
!  REAL d(1), v(1), amat(n,p)
!ON RETURN:
!  REAL amat(n,p)
!--------------------------------------------------------------------
SUBROUTINE dsqrtinvva(n, p, g, nsize, d, v, par, dec, amat)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, par, dec
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: d, v
   DOUBLE PRECISION, INTENT(INOUT) :: amat(n, p)
   INTEGER :: i, j, k, info 
   INTEGER :: l(g), u(g)
   DOUBLE PRECISION :: fd, sqrtv
   DOUBLE PRECISION :: m(p)
   DOUBLE PRECISION, ALLOCATABLE :: winvv(:, :)
   l(1) = 1
   u(1) = nsize(1)
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
      u(i) = u(i - 1) + nsize(i) 
   END DO
   IF (dec == 0) THEN
      SELECT CASE (par)
      CASE(1)  !sqrtinvv of type Hartley-Rao as default
         DO i = 1, g
            fd = (1 / SQRT(1 + d * nsize(i))) - 1
            m = fd * (SUM(amat(l(i):u(i), :), 1) / nsize(i))
            DO k = 1, p
               DO j = 1, nsize(i)
                  amat(l(i) + j - 1, k) = m(k) + amat(l(i) + j - 1, k)
               END DO
            END DO
         END DO
      CASE(0)  !default MLM sqrtinv * amat
         sqrtv = SQRT(v)
         DO i = 1, g
            fd = (1 / SQRT(1 + d * nsize(i))) - 1
            m = fd * (SUM(amat(l(i):u(i), :), 1) / nsize(i))
            DO k = 1, p
               DO j = 1, nsize(i)
                  amat(l(i) + j - 1, k) = (1/sqrtv) * (m(k) &
                     + amat(l(i) + j - 1, k))
               END DO
            END DO
         END DO
      CASE DEFAULT
         amat = 0D0 
      END SELECT  
   ELSE
      DO i = 1, g
         ALLOCATE( winvv( nsize(i), nsize(i) ) )
         winvv =  (- d) / ( 1D0 + ( d * nsize(i) ) )
         DO j = 1, nsize(i)
            winvv(j, j) = winvv(j, j) + 1
         END DO
         CALL dpotrf("U", nsize(i), winvv, nsize(i), info) 
         CALL dtrmm("L", "U", "N", "N", nsize(i), p, 1D0, winvv, &
            nsize(i), amat( l(i) : u(i), :), nsize(i)) 
         DEALLOCATE(winvv)
      END DO
      IF (par == 0) THEN
         DO j = 1, p
            DO i = 1, n
               amat(i, j) = amat(i, j) * SQRT(v)
            END DO
         END DO
      END IF
   END IF
END SUBROUTINE 
!
