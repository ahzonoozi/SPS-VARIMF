MODULE SUB_FUNC

  INTERFACE
     FUNCTION LOCATE(xx,x)
       USE params
       REAL(KIND(1.d0)), DIMENSION(:), INTENT(IN) :: xx
       REAL(KIND(1.d0)), INTENT(IN) :: x
       INTEGER :: locate
     END FUNCTION LOCATE
  END INTERFACE

  
  INTERFACE
     FUNCTION LINTERPARR(xin,yin,xout)
        USE params
        IMPLICIT NONE
        REAL(KIND(1.d0)), DIMENSION(:), INTENT(in) :: xin,yin
        REAL(KIND(1.d0)), INTENT(in), DIMENSION(:) :: xout
        REAL(KIND(1.d0)), DIMENSION(SIZE(xout)) :: linterparr
     END FUNCTION LINTERPARR
  END INTERFACE 

 
  INTERFACE
     FUNCTION LINTERP(xin,yin,xout)
       USE params
       REAL(KIND(1.d0)), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(KIND(1.d0)), INTENT(in)  :: xout
       REAL(KIND(1.d0)) :: linterp
     END FUNCTION LINTERP
  END INTERFACE
  
 

END MODULE SUB_FUNC
