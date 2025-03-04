SUBROUTINE max_mecl(SFR,Mecl_max)


  USE PARAMS
  IMPLICIT NONE

  REAL(KIND(1.d0)), INTENT(in) :: SFR!,FeH1
  REAL(KIND(1.d0)), INTENT(out) :: Mecl_max

  INTEGER :: f
  INTEGER :: nx
  REAL, DIMENSION(1000) :: Mecl
  REAL :: Bt,Mecl_u,Mecl_l,Kecl,Mtot,M_t



 
!-----------------------------------------------------

  nx=10000
  
  Bt =-0.106*log10(SFR)+2. 
  Mecl(1)=5.0      
  M_t=0.
  Mtot=SFR*10.**7.  !Mtot=SFR*dt
  Mecl_u=10.**9.
  Mecl_l=5.
       
   

  DO f=2,900
     M_t=0.
     Mecl(f)=Mecl_l+10**(0.01*f)
     IF (Bt==2.)THEN
        M_t=(log(Mecl(f))-log(Mecl_l))/(Mecl(f)**(-1.)-Mecl_u**(-1.))
     ELSE  
        M_t=((1.-Bt)/(2.-Bt))*(Mecl(f)**(2.-Bt)-Mecl_l**(2.-Bt))/(Mecl_u**(1.-Bt)-Mecl(f)**(1.-Bt))
     ENDIF

     IF (M_t.GE.Mtot) GOTO 10

  ENDDO

10   Mecl_max=Mecl(f)
     IF(M_t.LT.Mtot)THEN
        Mecl_max=1.e9
     ENDIF

     Kecl=Mtot*(-Bt+2.)/(Mecl_max**(-Bt+2.)-Mecl_L**(-Bt+2.0))
     
    
 
  RETURN
END

