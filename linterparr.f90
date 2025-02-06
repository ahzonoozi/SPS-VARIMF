FUNCTION LINTERPARR(xin, yin, xout)

  ! Routine to perform linear interpolation and extrapolation of a function yin(xin) at xout

  USE params
  USE sub_func, ONLY: locate
  IMPLICIT NONE

  REAL(KIND(1.d0)), DIMENSION(:), INTENT(IN) :: xin, yin
  REAL(KIND(1.d0)), INTENT(IN), DIMENSION(:) :: xout
  REAL(KIND(1.d0)), DIMENSION(SIZE(xout)) :: linterparr
  INTEGER :: klo, n, n2, i
  REAL(KIND(1.d0)) :: slope

  ! Validate input array sizes
  IF (SIZE(xin) /= SIZE(yin)) THEN
    PRINT *, "Error: xin and yin must have the same size."
    STOP
  END IF

  ! Get the size of input arrays
  n = SIZE(xin)
  n2 = SIZE(xout)

  ! Loop over each value in xout
  DO i = 1, n2
    IF (xout(i) < xin(1)) THEN
      ! Extrapolate below the range
      slope = (yin(2) - yin(1)) / (xin(2) - xin(1))
      linterparr(i) = yin(1) + slope * (xout(i) - xin(1))
    ELSE IF (xout(i) > xin(n)) THEN
      ! Extrapolate above the range
      slope = (yin(n) - yin(n-1)) / (xin(n) - xin(n-1))
      linterparr(i) = yin(n) + slope * (xout(i) - xin(n))
    ELSE
      ! Perform interpolation within the range
      klo = MAX(MIN(locate(xin, xout(i)), n - 1), 1)
      linterparr(i) = yin(klo) + (yin(klo + 1) - yin(klo)) * &
                      (xout(i) - xin(klo)) / (xin(klo + 1) - xin(klo))
    END IF
  ENDDO

END FUNCTION LINTERPARR

