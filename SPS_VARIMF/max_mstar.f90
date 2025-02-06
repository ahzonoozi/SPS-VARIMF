SUBROUTINE max_mstar(FeH,Mecl1,m_max1)

  USE PARAMS
  IMPLICIT NONE

  !compute the maximum stellar mass for each cluster
  REAL(KIND(1.d0)), INTENT(in) :: FeH, Mecl1
  REAL(KIND(1.d0)), INTENT(out) :: m_max1
  INTEGER :: j
  INTEGER :: nx,mx
  REAL, DIMENSION(400) :: m,dm,logM
  REAL :: M_C
  REAL(KIND(1.d0)) :: Bt,delta,M_LOW,M_turn1,M_turn2,M_UP,alpha3,k1,k2,k3
  REAL :: xx1

  !==================================================
   
   M_LOW=0.08
   M_turn1=0.5
   M_turn2=1.0
   M_UP=150.

   mx=1000

   !==================================================


  IF(Z_MODE==1)THEN
   alpha1=1.3+ 79.4*(z_isoc(zmet)-0.8*0.0142)
   alpha2=2.3+ 79.4*(z_isoc(zmet)-0.8*0.0142)
  ENDIF
  
  logM(1)=-1.08
  m(1)=10.0**(logM(1))  
             
  xx1=0.99*(0.61*log10(Mecl1)+2.85-6.0)-0.14*FeH
  IF (xx1.GE.-0.87 )THEN
     alpha3=(1.0*(1.94-0.41*xx1))
  ELSEIF(xx1.LT.-0.87)THEN
     alpha3=2.3
  ENDIF

  
  M_C=0.0     
  DO j=2,400

     logM(j)=-1.08+0.01*(j-1)  
     m(j)=10.0**(logM(j))
     dm(j)=m(j)-m(j-1)
     M_C=0.0                  


     IF (m(j).GE.0.08.AND.m(j).LT.0.5) THEN
        K1=1./((M_turn1**(1.-alpha1)-m(j)**(1.-alpha1))/(1.-alpha1)+(M_turn2**(1.-alpha2)-&
           &M_turn1**(1.-alpha2))*(M_turn1**(alpha2-alpha1))/(1.-alpha2)+(M_UP**(1.-alpha3)-&
           &M_turn2**(1.-alpha3))*(M_turn1**(alpha2-alpha1))*(M_turn2**(alpha3-alpha2))/(1.-alpha3)) 
                   
        M_C=K1*(m(j)**(2.-alpha1)-M_LOW**(2.-alpha1))/(2.-alpha1)              
     ENDIF 
     IF (m(j).GE.0.5.AND.m(j).LT.1.0) THEN
        K1=1./( (M_turn2**(1.-alpha2)-m(j)**(1.-alpha2))*(M_turn1**(alpha2-alpha1))/(1.-alpha2)+&
            &(M_UP**(1.-alpha3)-M_turn2**(1.-alpha3))*(M_turn1**(alpha2-alpha1))*&
            &(M_turn2**(alpha3-alpha2))/(1.-alpha3))
                   
        M_C=K1*((M_turn1**(2.-alpha1)-M_LOW**(2.-alpha1))/(2.-alpha1)+(m(j)**(2.-alpha2)-&
            &M_turn1**(2.-alpha2))*(M_turn1**(alpha2-alpha1))/(2.-alpha2))             
     ENDIF 
     IF (m(j).GE.1..AND.m(j).LT.150.) THEN
        K1=(M_turn1**(alpha1-alpha2))*(M_turn2**(alpha2-alpha3))*(1.-alpha3)/&
            &(M_UP**(1.-alpha3)-m(j)**(1.-alpha3))
                   
        M_C=K1*((M_turn1**(2.-alpha1)-M_LOW**(2.-alpha1))/(2.-alpha1)+(M_turn2**(2.-alpha2)-&
            &M_turn1**(2.-alpha2))*(M_turn1**(alpha2-alpha1))/(2.-alpha2)+(M(j)**(2.-alpha3)-&
            &M_turn2**(2.-alpha3))*(M_turn2**(alpha3-alpha2))*(M_turn1**(alpha2-alpha1))/(2.-alpha3))             
     ENDIF 
      
 
                   
     IF(M_C.GE.Mecl1) GOTO 20 
       
    
  ENDDO

20 m_max1=m(j)
   

  IF (M_C.LT.Mecl1)THEN
     m_max1=150.
  ENDIF

  RETURN   
END  



