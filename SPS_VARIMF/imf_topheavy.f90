SUBROUTINE IMF_TOPHEAVY(mini,nmass,N_IMF,Mret,Metal_eject,Nrem)


  USE PARAMS; USE sub_func
  IMPLICIT NONE


  REAL(KIND(1.d0)), INTENT(in), DIMENSION(nm)  :: mini
  INTEGER, INTENT(in) :: nmass
  REAL(KIND(1.d0)), INTENT(out), DIMENSION(nm) :: N_IMF
  REAL(KIND(1.d0)), INTENT(out), DIMENSION(nm) :: Mret,Metal_eject,Nrem
  REAL(KIND(1.d0)), DIMENSION(2000)   :: mrem
  INTEGER :: i,j
  REAL(KIND(1.d0)) :: m1,m2,IMF,dm
  REAL(KIND(1.d0)) :: k1,k2,k3
  REAL(KIND(1.d0)) :: m_max1
  REAL(KIND(1.d0)) :: M_LOW,M_turn1,M_turn2,M_UP,alpha3,Mass_isoc,Mrem_isoc
  REAL(KIND(1.d0)) FeH, xx

  REAL(KIND(1.d0)) :: B1,KK1,KK2,delta1,delta2,MCO,MM1,qq1
  REAL(KIND(1.d0)) :: P_MCO,F_MCO,A1,A2,L1,eta,h_MCO


  !----------------------------------------------------
  !---------------------------------------------------- 

   !OPEN(55,FILE='topheavy_IMF.txt')
   M_LOW=0.08
   M_turn1=0.5
   M_turn2=1.0
   M_UP=150.
   FeH=LOG10(z_isoc(zmet)/0.02)

   CALL max_mstar(FeH,M_UCD,m_max1)
   M_UP = m_max1

   !Top_heavy IMF for SCs or UCDs that form monotically 
   xx = 0.99*(0.61*log10(M_UCD)+2.85-6.0)-0.14*FeH


   IF (xx.GE.-0.87 )THEN
      alpha3=(1.0*(1.94-0.41* xx))
   ELSEIF(xx.LT.-0.87)THEN
      alpha3=2.3
   ENDIF





   Mass_isoc=0.0
   Mrem_isoc=0.0
   

  !------------------------------------------------------

  k2= 1./( (M_turn1**(alpha1-alpha2))*(M_turn1**(-alpha1+2.)-M_LOW**(-alpha1+2.))/(-alpha1+2.) +&
         &(M_turn2**(-alpha2+2.)-M_turn1**(-alpha2+2.))/(-alpha2+2.) +&
        &(M_turn2**(alpha3-alpha2))*(M_UP**(-alpha3+2.)-M_LOW**(-alpha3+2.))/(-alpha3+2.)  )
  K1 = k2 * M_turn1**(alpha1-alpha2)
  k3 = k2 * M_turn2**(alpha3-alpha2)

  DO i=1,nmass
     dm=0.0
     IMF=0.0
     N_IMF(i)=0.0

     IF (mini(i).GE.M_LOW.AND.mini(i).LT.M_turn1) THEN
        IMF=k1*mini(i)**(-alpha1)
     ENDIF                                       
     IF (mini(i).GE.M_turn1.AND.mini(i).LE.M_turn2) THEN
        IMF= k2*mini(i)**(-alpha2)    
     ENDIF
     IF (mini(i).GE.M_turn2.AND.mini(i).LE.M_UP) THEN
        IMF= k3*mini(i)**(-alpha3)
     ENDIF



     IF (i.EQ.1) THEN
        m1 = M_LOW
     ELSE
        m1 = mini(i) - 0.5*(mini(i)-mini(i-1))
     ENDIF
     IF (i.EQ.nmass) THEN
        m2 = mini(i)
     ELSE
        m2 = mini(i) + 0.5*(mini(i+1)-mini(i))
     ENDIF

        dm=m2-m1
        N_IMF(i) = IMF*dm
        !total live stars number and mass in each isoc
               
        !Num_isoc = IMF(i)*dm(i)+Num_isoc
        Mass_isoc = N_IMF(i)*mini(i)+Mass_isoc
         
               
  ENDDO 
 


   !------------------------------------------------
   !---------------add remnants---------------------
   !------------------------------------------------ 

DO j=1,2000
   Mret(j)=0.0
   mrem(j)=0.0
   Nrem(j)=0.0
   Metal_eject(j) =0.0


ENDDO


B1  = 0.0
A1  = 0.0
L1  = 0.0
KK1 = 0.0
kk2 = 0.0
MM1 =0.0
qq1 =0.0
delta1 = 0.0
delta2 =0.0
MCO = 0.0
P_MCO = 0.0
F_MCO = 0.0
h_MCO = 0.0




IF (mini(nmass).LT.M_UP) THEN

   DO j=1,num_rem
      mrem(j)=mini(nmass)+j*(M_Up-mini(nmass))/num_rem
      Metal_eject(j) = Metal_yield(max(locate(z_yield,z_isoc(zmet)),1),locate(mass_yield,mrem(j)))

!-------------------------------
      IF(remnant_cal.EQ.1)THEN



         IF (z_isoc(zmet).GT. 4.*10.**(-3)) THEN
            B1 = 59.63-2.969*(10.**3)*z_isoc(zmet)+4.988*(10.**4)*z_isoc(zmet)**2.
            KK1 = 45.04-2.176*(10.**3)*z_isoc(zmet)+3.806*(10.**4)*(z_isoc(zmet)**2.)
            kK2 = 1.389*(10.**2)-4.664*(10.**3)*z_isoc(zmet)+5.106*(10.**4)*z_isoc(zmet)**2.
            delta1 = 2.790*10.**(-2)-1.780*(10.**(-2))*z_isoc(zmet)+77.05*z_isoc(zmet)**2.
           delta2 = 6.730*10.**(-3)+2.690*z_isoc(zmet)-52.39*z_isoc(zmet)**2.

ELSEIF (z_isoc(zmet).GE.10.**(-3) .AND. z_isoc(zmet).LE.4.*10.**(-3)) THEN
B1 = 40.98+3.415*(10.**4)*z_isoc(zmet)-8.064*(10.**6)*(z_isoc(zmet)**2.)
KK1 = 35.17+1.548*(10.**4.)*z_isoc(zmet)-3.759*(10.**6)*(z_isoc(zmet)**2.)
KK2 = 20.36+1.162*(10.**5.)*z_isoc(zmet)-2.276*(10.**7)*(z_isoc(zmet)**2.)
delta1 = 2.5*(10.**(-2))-4.346*z_isoc(zmet)+1.340*(10.**3)*(z_isoc(zmet)**2.)
delta2 = 1.750*(10.**(-2.))+11.39*z_isoc(zmet)-2.902*(10.**3)*(z_isoc(zmet)**2.)

ELSEIF (z_isoc(zmet).LT.10.**(-3.)) THEN
B1 = 67.07
KK1 = 46.89
KK2 = 1.138*(10.**2.)
delta1 = 2.199*(10.**(-2.))
delta2 = 2.602*(10.**(-2.))
ENDIF

MCO = -2.0+(B1+2.0)*(0.5/(1.+10.**((KK1-Mrem(j))*delta1)))+&
(0.5/(1.+10.**((KK2-Mrem(j))*delta2)))

IF(z_isoc(zmet) .LE.5.*10.**(-4)) THEN
MM1 = -6.476*(10.**2)*z_isoc(zmet)+1.911
qq1 = 2.3*(10.**3)*z_isoc(zmet)+11.67
P_MCO = -2.333+0.1559 * MCO+0.27*(MCO**2.)
F_MCO = MM1*MCO+qq1


IF (MCO.LE.5.) THEN
Mret(j) = MAX(p_MCO, 1.27)
ELSEIF (MCO.GT. 5. .AND. MCO.LT.10.0) THEN
Mret(j) = P_MCO
ELSEIF (MCO.GE.10.) THEN
Mret(j) = MIN(P_MCO, F_MCO)
ENDIF

ELSEIF (z_isoc(zmet).GT.5.*10.**(-4.)) THEN

IF (z_isoc(zmet).GE.10.**(-3.)) THEN
A1 = 1.34- 29.46/(1.+(z_isoc(zmet)/(1.11*10.**(-3)))**2.36)
A2 = 80.22- 74.73*(z_isoc(zmet)**0.965)/(2.72*10.**(-3)+(z_isoc(zmet)**0.965))
L1 = 5.683+3.533/(1.+(z_isoc(zmet)/(7.43*10.**(-3)))**1.993)
eta = 1.066- 1.121/(1.+(z_isoc(zmet)/(2.558*10.**(-2)))**0.609)
ELSEIF (z_isoc(zmet).LT. 10.**(-3)) THEN
A1 = 1.105*(10.**5)*z_isoc(zmet)-1.258*10.**2.
A2 = 91.56-1.957*(10.**4.)*z_isoc(zmet)-1.558*(10.**7)*(z_isoc(zmet)**2.)
L1 = 1.134*(10.**4.)*z_isoc(zmet)-2.143
eta = 3.09*(10.**(-2.))-22.3*z_isoc(zmet)+7.363*(10.**4)*z_isoc(zmet)**2.
ENDIF
IF (z_isoc(zmet).GE.2.*10.**(-3)) THEN
MM1 = 1.217
qq1 = 1.061
ELSE IF (z_isoc(zmet) .LT.2.*10.**(-3).AND.z_isoc(zmet).GE.10.**(-3.)) THEN
MM1 = -43.82*z_isoc(zmet)+1.304
qq1 = -1.296*(10.**4)*z_isoc(zmet)+26.98
ELSE IF(z_isoc(zmet).LT.10.**(-3)) THEN
MM1 = -6.476*(10.**2)*z_isoc(zmet)+1.911
qq1 = 2.3*(10.**3)*z_isoc(zmet)+11.67
ENDIF

h_MCO = A1+(A2-A1)/(1.+10.**((L1-MCO)*eta))
f_MCO = MM1*MCO + qq1


IF (MCO.LE.5.) THEN
Mret(j) = max(h_MCO, 1.27)
ELSEIF (MCO.GT. 5. .AND. MCO .LT.10.0) THEN
Mret(j) = h_MCO
ELSEIF (MCO.GE.10.) THEN
Mret(j) = max(h_MCO,F_MCO)
ENDIF


ENDIF

IF(mrem(j).LE. 2.6) THEN
Mret(j)= 0.077*mrem(j)+0.48
ENDIF


!-------------------------------

ELSE IF(remnant_cal.EQ.0)THEN

              IF (mrem(j).GE.40.) THEN
                 Mret(j) = 0.5*mrem(j)
              ENDIF
              IF (mrem(j).GE.8.5.AND.mrem(j).LT.40.) THEN
                 Mret(j) = 1.4
              ENDIF
              IF (mrem(j).LT.8.5) THEN
                 Mret(j)=0.077*mrem(j)+0.48
              ENDIF

ENDIF
ENDDO

!--------------------

           DO j=1,num_rem
              dm=0.0
              IMF=0.0
              Nrem(j)=0.0

              IF (mrem(j).GE.M_LOW.AND.mrem(j).LT.M_turn1) THEN
                  IMF = k1*mrem(j)**(-alpha1)
               ENDIF
               IF (mrem(j).GE.M_turn1.AND.mrem(j).LE.M_turn2) THEN
                  IMF = k2*mrem(j)**(-alpha2)
               ENDIF
               IF (mrem(j).GE.M_turn2.AND.mrem(j).LE.M_UP) THEN
                  IMF= k3*mrem(j)**(-alpha3)
               ENDIF



               IF (j.EQ.1) THEN
                  m1 = mini(nmass)
               ELSE
                  m1 = mrem(j) - 0.5*(mrem(j)-mrem(j-1))
               ENDIF
               IF (j.EQ. 2000) THEN
                  m2 = mrem(j)
               ELSE
                  m2 = mrem(j) + 0.5*(mrem(j+1)-mrem(j))
               ENDIF
                
               dm=m2-m1

               Nrem(j) = IMF*dm

               !total live number and mass in each isoc 
               !Num_isoc = IMF(i)*dm(i)+Num_isoc
               Mrem_isoc = Nrem(j)*mrem(j)+Mrem_isoc
           ENDDO
        ENDIF



        !The number of stars are normalized to the isochrone mass #normalization is checked
        DO i=1,nmass
       
           N_IMF(i) = N_IMF(i)/(Mrem_isoc+Mass_isoc)                     
        ENDDO 

        DO j=1,2000

           Nrem(j) = Nrem(j)/(Mrem_isoc+Mass_isoc)           
        ENDDO 

     
  RETURN
END  !SUBROUTINE IMF_TOPHEAVY


