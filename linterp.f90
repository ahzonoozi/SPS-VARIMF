FUNCTION LINTERP(xin, yin, xout)
  ! Routine to perform linear interpolation and extrapolation of a function yin(xin) at xout

  USE params
  USE sub_func, ONLY: locate
  IMPLICIT NONE

  REAL(KIND(1.d0)), DIMENSION(:), INTENT(IN) :: xin, yin
  REAL(KIND(1.d0)), INTENT(IN) :: xout
  REAL(KIND(1.d0)) :: linterp
  INTEGER :: klo, n
  REAL(KIND(1.d0)) :: slope

  ! Validate input array sizes
  IF (SIZE(xin) /= SIZE(yin)) THEN
    PRINT *, "Error: xin and yin must have the same size."
    STOP
  END IF

  ! Get the size of the input array
  n = SIZE(xin)

  ! Determine if xout is outside the range of xin for extrapolation
  IF (xout < xin(1)) THEN
    ! Extrapolate below the range
    slope = (yin(2) - yin(1)) / (xin(2) - xin(1))
    linterp = yin(1) + slope * (xout - xin(1))
  ELSE IF (xout > xin(n)) THEN
    ! Extrapolate above the range
    slope = (yin(n) - yin(n-1)) / (xin(n) - xin(n-1))
    linterp = yin(n) + slope * (xout - xin(n))
  ELSE
    ! Perform interpolation within the range
    klo = MAX(MIN(locate(xin, xout), n - 1), 1)
    linterp = yin(klo) + (yin(klo + 1) - yin(klo)) * &
              (xout - xin(klo)) / (xin(klo + 1) - xin(klo))
  END IF

END FUNCTION LINTERP

