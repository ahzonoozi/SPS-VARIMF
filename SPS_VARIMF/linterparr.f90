FUNCTION LINTERPARR(xin,yin,xout)

  !routine to linearly interpolate a function yin(xin) at xout

  USE params
  USE sub_func, ONLY: locate
  IMPLICIT NONE
  REAL(KIND(1.d0)), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(KIND(1.d0)), INTENT(in), DIMENSION(:) :: xout
  REAL(KIND(1.d0)), DIMENSION(SIZE(xout)) :: linterparr
  INTEGER :: klo,n,n2,i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  n2  = SIZE(xout)

  DO i=1,n2
     klo = MAX(MIN(locate(xin,xout(i)),n-1),1)
     linterparr(i) = yin(klo) + (yin(klo+1)-yin(klo))*&
          (xout(i)-xin(klo))/(xin(klo+1)-xin(klo))
  ENDDO

END FUNCTION LINTERPARR
