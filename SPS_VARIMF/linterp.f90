FUNCTION LINTERP(xin,yin,xout)

  !routine to linearly interpolate a function yin(xin) at xout

  USE params
  USE sub_func, ONLY: locate
  IMPLICIT NONE
  REAL(KIND(1.d0)), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(KIND(1.d0)), INTENT(in)  :: xout
  REAL(KIND(1.d0)) :: linterp
  INTEGER :: klo,n

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  klo = MAX(MIN(locate(xin,xout),n-1),1)

  linterp = yin(klo) + (yin(klo+1)-yin(klo))*&
       (xout-xin(klo))/(xin(klo+1)-xin(klo))

END FUNCTION LINTERP
