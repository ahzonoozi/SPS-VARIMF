SUBROUTINE GETSPEC(mact,logt,lbol,logg,phase,ffco,spec)

  ! Routine to compute the spectrum of a star based on input parameters (Z,logg and logt)
  ! Inputs: 
  !   - mact: mass of the star
  !   - logt: log(Teff)
  !   - lbol: luminosity (L/Lsun)
  !   - logg: log(surface gravity)
  !   - phase: stellar evolutionary phase
  !   - ffco: chemical composition (for WR and TP-AGB stars)
  ! Outputs:
  !   - spec: computed spectrum (Fnu)
 
  USE params; USE sub_func
  IMPLICIT NONE

  ! Inputs
  REAL(KIND(1.d0)), INTENT(IN) :: mact, logt, lbol, logg, phase, ffco
  REAL(KIND(1.d0)), INTENT(INOUT), DIMENSION(nspec) :: spec

  ! Internal variables
  INTEGER :: i_t, i_g, i_z, l
  REAL(KIND(1.d0)), DIMENSION(nspec) :: spec_1, spec_2
  REAL(KIND(1.d0)) :: r2, aa_t, aa_g, aa_z

  ! Constants
  REAL(KIND(1.d0)), PARAMETER :: pi_value = 3.14159265
  !Gravitational constant [cgs]
  REAL(KIND(1.d0)), PARAMETER :: G_const = 6.67428E-8    
  !speed of light [Ang/s]
  REAL(KIND(1.d0)), PARAMETER :: clight = 2.9979E18  
  !Solar mass [g]   
  REAL(KIND(1.d0)), PARAMETER :: msun = 1.989E33  
  !Solar luminosity [erg/s]      
  REAL(KIND(1.d0)), PARAMETER :: lsun = 3.839E33         



  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  ! Initialize spectrum to a small value
  spec  = small_value

  !----------------------------------------------------------------

  !post-AGB NLTE model spectra from Rauch 2003

  IF (phase.EQ.9.0.AND.logt.GE.4.699) THEN


     i_t = locate(pagb_logt,logt)
        IF (i_t .LT. 1)   i_t = 1
        IF (i_t .GT. nt_pagb-1)  i_t = nt_pagb-1

     i_g = locate(pagb_logg,logg)
        IF (i_g .LT. 1)   i_g = 1
        IF (i_g .GT. ng_pagb-1)  i_g = ng_pagb-1

     aa_t = (logt-pagb_logt(i_t)) /(pagb_logt(i_t+1)-pagb_logt(i_t))
     aa_g = (logg-pagb_logg(i_g)) /(pagb_logg(i_g+1)-pagb_logg(i_g))

     IF (aa_t .LT. 0.0)  aa_t=0.0
     IF (aa_t .GT. 1.0)  aa_t=1.0
     IF (aa_g .LT. 0.0)  aa_g=0.0
     IF (aa_g .GT. 1.0)  aa_g=1.0
     !bb_t = 1-aa

     l=1
     IF (z_isoc(zmet)/zsun_isoc.GT.0.5) l=2

     spec_1 = (1-aa_t)*pagb_spec(:,i_t,i_g,l)+aa_t*pagb_spec(:,i_t+1,i_g,l)
     spec_2 = (1-aa_t)*pagb_spec(:,i_t,i_g+1,l)+aa_t*pagb_spec(:,i_t+1,i_g+1,l)

     spec = lbol*((1-aa_g)*spec_1+aa_g*spec_2)

 !--------------------------------------------------------------

  !WR library from Smith et al. 2002 (CMFGEN)
  ELSE IF (phase.EQ.10.0) THEN

     IF (ffco.LT.10) THEN

        i_t = locate(wrn_logt,logt)
        IF (i_t .LT. 1)   i_t = 1
        IF (i_t .GT. nt_wr-1)  i_t = nt_wr-1

        aa_t = (logt-wrn_logt(i_t))/(wrn_logt(i_t+1)-wrn_logt(i_t))
        IF (aa_t .LT. 0.0)  aa_t=0.0
        IF (aa_t .GT. 1.0)  aa_t=1.0

        spec = lbol*((1-aa_t)*wrn_spec(:,i_t,zmet)+&
             aa_t*wrn_spec(:,i_t+1,zmet))

     ELSE IF (ffco.GE.10) THEN
     

        i_t = locate(wrc_logt,logt)
            IF (i_t .LT. 1)   i_t = 1
            IF (i_t .GT. nt_wr-1)  i_t = nt_wr-1

        aa_t = (logt-wrc_logt(i_t))/(wrc_logt(i_t+1)-wrc_logt(i_t))

        IF (aa_t .LT. 0.0)  aa_t=0.0
        IF (aa_t .GT. 1.0)  aa_t=1.0

        spec = lbol*((1-aa_t)*wrc_spec(:,i_t,zmet)+&
             aa_t*wrc_spec(:,i_t+1,zmet))


     ENDIF

!---------------------------------------------------------------------------------

  !O-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.7.0.OR.phase.EQ.8.0.AND.logt.LT.3.6.AND.ffco.LE.1.0) THEN
     


     i_t =  locate(agb_logt_o(:,zmet),logt)
         IF (i_t .LT. 1)   i_t = 1
         IF (i_t .GT. n_agb_o-1)  i_t = n_agb_o-1


     aa_t  = (logt - agb_logt_o(i_t,zmet)) / &
             (agb_logt_o(i_t+1,zmet)-agb_logt_o(i_t,zmet))

     IF (aa_t .LT. 0.0)  aa_t=0.0
     IF (aa_t .GT. 1.0)  aa_t=1.0



     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.

     spec = lbol * lambda*lambda/clight * &
          ( (1.-aa_t)*agb_spec_o(:,i_t) + aa_t*(agb_spec_o(:,i_t+1)) )


  
  !C-rich TP-AGB spectra
   ELSE IF (phase.EQ.7.0.OR.phase.EQ.8.0.AND.logt.LT.3.6.AND.ffco.GT.1.0) THEN


        i_t  =  locate(agb_logt_c,logt)
            IF (i_t .LT. 1)   i_t = 1
            IF (i_t .GT. n_agb_c-1)  i_t = n_agb_c-1



        aa_t =  (logt - agb_logt_c(i_t)) / &
                (agb_logt_c(i_t+1)-agb_logt_c(i_t))

        IF (aa_t .LT. 0.0)  aa_t=0.0
        IF (aa_t .GT. 1.0)  aa_t=1.0


        !The spectra are Fdlambda, need to convert to Fdnu and
        spec = lbol* lambda*lambda/clight * &
               ( (1.-aa_t)*agb_spec_c(:,i_t) + aa_t*(agb_spec_c(:,i_t+1)) )


      
!------------------------------------------------------------------------

  !use WMBasic grid from JJ Eldridge for T>25,000K MS stars
  ELSE IF (phase.LE.1.0.AND.logt.GT.4.397) THEN



     i_t = locate(wmb_logt,logt)
         IF (i_t .LT. 1)   i_t = 1
         IF (i_t .GT. nt_wmb-1)  i_t = nt_wmb-1


    i_z = locate(wmb_z,z_isoc(zmet))
        IF (i_z .LT. 1)   i_z = 1
        IF (i_z .GT. nz_wmb-1)  i_z = nz_wmb-1




     aa_t   = (logt-wmb_logt(i_t)) / (wmb_logt(i_t+1)-wmb_logt(i_t))
     IF (aa_t .LT. 0.0)  aa_t=0.0
     IF (aa_t .GT. 1.0)  aa_t=1.0


     aa_z   = (z_isoc(zmet)-wmb_z(i_z)) / (wmb_z(i_z+1)-wmb_z(i_z))
     IF (aa_z .LT. 0.0)  aa_z=0.0
     IF (aa_z .GT. 1.0)  aa_z=1.0



    IF (logg.LT.3.5) l=1
    IF (logg.GT.3.5 .AND. logg.LT.4.25) l=2
    IF (logg.GT.4.25) l=3

    spec_1 = (1-aa_t)*wmb_spec(:,i_z,i_t,l)+aa_t*wmb_spec(:,i_z,i_t+1,l)
    spec_2 = (1-aa_t)*wmb_spec(:,i_z+1,i_t,l)+aa_t*wmb_spec(:,i_z+1,i_t+1,l)


     spec = lbol*((1-aa_z)*spec_1+aa_z*spec_2)

  ELSE

!________________________________________________________________________________

  !use the main spectral library for the rest of stars


     i_t = locate(speclib_logt,logt)
         IF (i_t .LT. 1)   i_t = 1
         IF (i_t .GT. ndim_logt-1)  i_t = ndim_logt-1

     i_g = locate(speclib_logg,logg)
         IF (i_g .LT. 1)   i_g = 1
         IF (i_g .GT. ndim_logg-1)  i_g = ndim_logg-1


     aa_t =(logt-speclib_logt(i_t))/(speclib_logt(i_t+1)-speclib_logt(i_t))

     IF (aa_t .LT. 0.0)  aa_t=0.0
     IF (aa_t .GT. 1.0)  aa_t=1.0

     !
     aa_g = (logg-speclib_logg(i_g))/(speclib_logg(i_g+1)-speclib_logg(i_g))

     IF (aa_g .LT. 0.0)  aa_g=0.0
     IF (aa_g .GT. 1.0)  aa_g=1.0


     spec_1 = (1-aa_t)*spec_flux(:,zmet,i_t,i_g)+aa_t*spec_flux(:,zmet,i_t+1,i_g)
     spec_2 = (1-aa_t)*spec_flux(:,zmet,i_t,i_g+1)+aa_t*spec_flux(:,zmet,i_t+1,i_g+1)

     r2  = mact*msun*G_const/10**logg
     spec = (((4.*pi_value)**2.)*r2/lsun)*( (1.-aa_g)*spec_1 + aa_g*spec_2 )


  ENDIF

  ! Ensure the spectrum has no negative values
  spec = MAX(spec,small_value)


END SUBROUTINE GETSPEC
