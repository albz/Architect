MODULE shapiro_wilks
  IMPLICIT NONE
  CONTAINS

SUBROUTINE swilk (init, x, n, n1, n2, a, w, pw, ifault)
!        ALGORITHM APPL. STATIST. (1995) VOL.44, NO.4
!        Calculates the Shapiro-Wilk W test and its significance level
! ARGUMENTS:
!   INIT     Set to .FALSE. on the first call so that weights A(N2) can be
!            calculated.   Set to .TRUE. on exit unless IFAULT = 1 or 3.
!   X(N1)    Sample values in ascending order.
!   N        The total sample size (including any right-censored values).
!   N1       The number of uncensored cases (N1 <= N).
!   N2       Integer part of N/2.
!   A(N2)    The calculated weights.
!   W        The Shapiro-Wilks W-statistic.
!   PW       The P-value for W.
!   IFAULT   Error indicator:
!            = 0 for no error
!            = 1 if N1 < 3
!            = 2 if N > 5000 (a non-fatal error)
!            = 3 if N2 < N/2
!            = 4 if N1 > N or (N1 < N and N < 20).
!            = 5 if the proportion censored (N - N1)/N > 0.8.
!            = 6 if the data have zero range.
!            = 7 if the X's are not sorted in increasing order

IMPLICIT NONE

LOGICAL, INTENT(IN OUT)  :: init
REAL(8), INTENT(IN)         :: x(:)
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: n1
INTEGER, INTENT(IN)      :: n2
REAL(8), INTENT(OUT)        :: a(:)
REAL(8), INTENT(OUT)        :: w
REAL(8), INTENT(OUT)        :: pw
INTEGER, INTENT(OUT)     :: ifault

! Local variables

REAL(8)    :: summ2, ssumm2, fac, rsn, an, an25, a1, a2, delta, range
REAL(8)    :: sa, sx, ssx, ssa, sax, asa, xsx, ssassx, w1, y, xx, xi
REAL(8)    :: gamma, m, s, ld, bf, z90f, z95f, z99f, zfm, zsd, zbar
INTEGER :: ncens, nn2, i, i1, j

LOGICAL, SAVE   :: upper = .TRUE.
REAL(8), PARAMETER :: c1(6) = (/ 0.0, 0.221157, -0.147981, -2.07119, 4.434685, &
                             -2.706056 /)
REAL(8), PARAMETER :: c2(6) = (/ 0.0, 0.042981, -0.293762, -1.752461, 5.682633, &
                             -3.582633 /)
REAL(8), PARAMETER :: c3(4) = (/ 0.5440, -0.39978, 0.025054, -0.6714E-3 /)
REAL(8), PARAMETER :: c4(4) = (/ 1.3822, -0.77857, 0.062767, -0.0020322 /)
REAL(8), PARAMETER :: c5(4) = (/ -1.5861, -0.31082, -0.083751, 0.0038915 /)
REAL(8), PARAMETER :: c6(3) = (/ -0.4803, -0.082676, 0.0030302 /)
REAL(8), PARAMETER :: c7(2) = (/ 0.164, 0.533 /)
REAL(8), PARAMETER :: c8(2) = (/ 0.1736, 0.315 /)
REAL(8), PARAMETER :: c9(2) = (/ 0.256, -0.00635 /)
REAL(8), PARAMETER :: g(2)  = (/ -2.273, 0.459 /)
REAL(8), PARAMETER :: z90 = 1.2816, z95 = 1.6449, z99 = 2.3263
REAL(8), PARAMETER :: zm = 1.7509, zss = 0.56268
REAL(8), PARAMETER :: bf1 = 0.8378, xx90 = 0.556, xx95 = 0.622
REAL(8), PARAMETER :: zero = 0.0, one = 1.0, two = 2.0, three = 3.0
REAL(8), PARAMETER :: sqrth = 0.70711, qtr = 0.25, th = 0.375, small = 1E-19
REAL(8), PARAMETER :: pi6 = 1.909859, stqr = 1.047198

pw  =  one
IF (w >= zero) w = one
an = n
ifault = 3
nn2 = n/2
IF (n2 < nn2) RETURN
ifault = 1
IF (n < 3) RETURN

!        If INIT is false, calculates coefficients for the test

IF (.NOT. init) THEN
  IF (n == 3) THEN
    a(1) = sqrth
  ELSE
    an25 = an + qtr
    summ2 = zero
    DO i = 1, n2
      CALL ppnd7((i - th)/an25, a(i), ifault)
      summ2 = summ2 + a(i) ** 2
    END DO
    summ2 = summ2 * two
    ssumm2 = SQRT(summ2)
    rsn = one / SQRT(an)
    a1 = poly(c1, 6, rsn) - a(1) / ssumm2

!        Normalize coefficients

    IF (n > 5) THEN
      i1 = 3
      a2 = -a(2)/ssumm2 + poly(c2,6,rsn)
      fac = SQRT((summ2 - two * a(1) ** 2  &
            - two * a(2) ** 2)/(one - two * a1 ** 2 - two * a2 ** 2))
      a(1) = a1
      a(2) = a2
    ELSE
      i1 = 2
      fac = SQRT((summ2 - two * a(1) ** 2)/ (one - two * a1 ** 2))
      a(1) = a1
    END IF
    DO i = i1, nn2
      a(i) = -a(i)/fac
    END DO
  END IF
  init = .true.
END IF
IF (n1 < 3) RETURN
ncens = n - n1
ifault = 4
IF (ncens < 0 .OR. (ncens > 0 .AND. n < 20)) RETURN
ifault = 5
delta = REAL(ncens)/an
IF (delta > 0.8) RETURN

!        If W input as negative, calculate significance level of -W

IF (w < zero) THEN
  w1 = one + w
  ifault = 0
  GO TO 70
END IF

!        Check for zero range

ifault = 6
range = x(n1) - x(1)
IF (range < small) RETURN

!        Check for correct sort order on range - scaled X

ifault = 7
xx = x(1)/range
sx = xx
sa = -a(1)
j = n - 1
DO i = 2, n1
  xi = x(i)/range
  IF (xx-xi > small) THEN
    WRITE(*, *) 'x(i)s out of order'
    RETURN
  END IF
  sx = sx + xi
  IF (i /= j) sa = sa + SIGN(1, i - j) * a(MIN(i, j))
  xx = xi
  j = j - 1
END DO
ifault = 0
IF (n > 5000) ifault = 2

!        Calculate W statistic as squared correlation
!        between data and coefficients

sa = sa/n1
sx = sx/n1
ssa = zero
ssx = zero
sax = zero
j = n
DO i = 1, n1
  IF (i /= j) THEN
    asa = SIGN(1, i - j) * a(MIN(i, j)) - sa
  ELSE
    asa = -sa
  END IF
  xsx = x(i)/range - sx
  ssa = ssa + asa * asa
  ssx = ssx + xsx * xsx
  sax = sax + asa * xsx
  j = j - 1
END DO

!        W1 equals (1-W) claculated to avoid excessive rounding error
!        for W very near 1 (a potential problem in very large samples)

ssassx = SQRT(ssa * ssx)
w1 = (ssassx - sax) * (ssassx + sax)/(ssa * ssx)
70 w = one - w1

!        Calculate significance level for W (exact for N=3)

IF (n == 3) THEN
  pw = pi6 * (ASIN(SQRT(w)) - stqr)
  RETURN
END IF
y = LOG(w1)
xx = LOG(an)
m = zero
s = one
IF (n <= 11) THEN
  gamma = poly(g, 2, an)
  IF (y >= gamma) THEN
    pw = small
    RETURN
  END IF
  y = -LOG(gamma - y)
  m = poly(c3, 4, an)
  s = EXP(poly(c4, 4, an))
ELSE
  m = poly(c5, 4, xx)
  s = EXP(poly(c6, 3, xx))
END IF
IF (ncens > 0) THEN

!        Censoring by proportion NCENS/N.  Calculate mean and sd
!        of normal equivalent deviate of W.

  ld = -LOG(delta)
  bf = one + xx * bf1
  z90f = z90 + bf * poly(c7, 2, xx90 ** xx) ** ld
  z95f = z95 + bf * poly(c8, 2, xx95 ** xx) ** ld
  z99f = z99 + bf * poly(c9, 2, xx) ** ld

!        Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
!        pseudo-mean and pseudo-sd of z as the slope and intercept

  zfm = (z90f + z95f + z99f)/three
  zsd = (z90*(z90f-zfm)+z95*(z95f-zfm)+z99*(z99f-zfm))/zss
  zbar = zfm - zsd * zm
  m = m + zbar * s
  s = s * zsd
END IF
pw = alnorm((y - m)/s, upper)

RETURN
END SUBROUTINE swilk


SUBROUTINE ppnd7 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477- 484.

! Produces the normal deviate Z corresponding to a given lower tail area of P;
! Z is accurate to about 1 part in 10**7.

! The hash sums below are the sums of the mantissas of the coefficients.
! They are included for use in checking transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

!INTEGER, PARAMETER           :: low_prec = SELECTED_REAL_KIND(6, 30)
REAL (8), INTENT(IN)  :: p
INTEGER, INTENT(OUT)         :: ifault
REAL (8), INTENT(OUT) :: normal_dev

! Local variables

REAL (8) :: zero = 0.0, one = 1.0, half = 0.5, split1 = 0.425,  &
                   split2 = 5.0, const1 = 0.180625, const2 = 1.6, q, r

! Coefficients for P close to 0.5

REAL (8) :: a0 = 3.3871327179E+00, a1 = 5.0434271938E+01, &
                   a2 = 1.5929113202E+02, a3 = 5.9109374720E+01, &
                   b1 = 1.7895169469E+01, b2 = 7.8757757664E+01, &
                   b3 = 6.7187563600E+01
! HASH SUM AB          32.3184577772

! Coefficients for P not close to 0, 0.5 or 1.

REAL (8) :: c0 = 1.4234372777E+00, c1 = 2.7568153900E+00, &
                   c2 = 1.3067284816E+00, c3 = 1.7023821103E-01, &
                   d1 = 7.3700164250E-01, d2 = 1.2021132975E-01
! HASH SUM CD          15.7614929821

! Coefficients for P near 0 or 1.

REAL (8) :: e0 = 6.6579051150E+00, e1 = 3.0812263860E+00, &
                   e2 = 4.2868294337E-01, e3 = 1.7337203997E-02, &
                   f1 = 2.4197894225E-01, f2 = 1.2258202635E-02
! HASH SUM EF          19.4052910204

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / &
               (((b3 * r + b2) * r + b1) * r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
  ELSE
    r = r - split2
    normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF

RETURN
END SUBROUTINE ppnd7



FUNCTION alnorm(x, upper) RESULT(fn_val)
!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

IMPLICIT NONE
!INTEGER, PARAMETER    :: sp = SELECTED_REAL_KIND(6, 30)
REAL (8), INTENT(IN) :: x
LOGICAL, INTENT(IN)   :: upper
REAL (8)             :: fn_val

! Local variables
REAL (8), PARAMETER :: zero = 0.0, one = 1.0, half = 0.5, &
                        con = 1.28
REAL (8) :: z, y
LOGICAL   :: up

!*** machine dependent constants
REAL (8), PARAMETER :: ltone = 7.0, utzero = 18.66

REAL (8), PARAMETER :: p = 0.398942280444, q = 0.39990348504,   &
                        r = 0.398942280385, a1 = 5.75885480458,  &
                        a2 = 2.62433121679, a3 = 5.92885724438,  &
                        b1 = -29.8213557807, b2 = 48.6959930692, &
                        c1 = -3.8052E-8, c2 = 3.98064794E-4,     &
                        c3 = -0.151679116635, c4 = 4.8385912808, &
                        c5 = 0.742380924027, c6 = 3.99019417011, &
                        d1 = 1.00000615302, d2 = 1.98615381364,  &
                        d3 = 5.29330324926, d4 = -15.1508972451, &
                        d5 = 30.789933034

up = upper
z = x
IF(z >=  zero) GO TO 10
up = .NOT. up
z = -z
10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
fn_val = zero
GO TO 40
20 y = half*z*z
IF(z > con) GO TO 30

fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
GO TO 40
30 fn_val = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40 IF(.NOT. up) fn_val = one - fn_val

RETURN
END FUNCTION alnorm



FUNCTION poly(c, nord, x) RESULT(fn_val)

!        Algorithm AS 181.2   Appl. Statist.  (1982) Vol. 31, No. 2

!        Calculates the algebraic polynomial of order nored-1 with
!        array of coefficients c.  Zero order coefficient is c(1)

REAL(8), INTENT(IN)  :: c(:)
INTEGER, INTENT(IN)  :: nord
REAL(8), INTENT(IN)  :: x
REAL(8)              :: fn_val

! Local variables

INTEGER :: i, j, n2
REAL(8) :: p

fn_val = c(1)
IF (nord == 1) RETURN
p = x*c(nord)
IF (nord == 2) GO TO 20
n2 = nord - 2
j = n2 + 1
DO i = 1,n2
  p = (p + c(j))*x
  j = j - 1
END DO
20 fn_val = fn_val + p

RETURN
END FUNCTION poly




RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
REAL(8), DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list), list, order)
RETURN
END SUBROUTINE quick_sort




RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end, list, order)

INTEGER, INTENT(IN) :: left_end, right_end
REAL(8), INTENT(INOUT) :: list(right_end)
INTEGER, INTENT(INOUT) :: order(right_end)

!     Local variables
INTEGER             :: i, j, itemp
REAL(8)             :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end, list, order)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1
  j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j, list, order)
  IF (i < right_end) CALL quick_sort_1(i, right_end, list, order)
END IF

RETURN
END SUBROUTINE quick_sort_1



SUBROUTINE interchange_sort(left_end, right_end, list, order)

INTEGER, INTENT(IN) :: left_end, right_end
REAL(8), INTENT(INOUT) :: list(right_end)
INTEGER, INTENT(INOUT) :: order(right_end)

!     Local variables
INTEGER   :: i, j, itemp
REAL(8)      :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    END IF
  END DO
END DO

RETURN
END SUBROUTINE interchange_sort


END MODULE shapiro_wilks
