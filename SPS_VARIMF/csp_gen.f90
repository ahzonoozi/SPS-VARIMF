SUBROUTINE csp_gen(jj,zmet_time,element_time,tage,mass_csp,mrem_csp,mdwarf_csp,lbol_csp,spec_csp,&
                   n_Ostars,n_Bstars,Metal_gas_csp,Metal_star_csp,element_gas_csp,element_star_csp,Mgas_csp,SNIa_csp)



  USE params; USE sub_func
  IMPLICIT NONE


  ! Inputs
  REAL(KIND(1.d0)), INTENT(IN) :: tage
  INTEGER , INTENT(IN), DIMENSION(num_time) ::  Zmet_time
  REAL, INTENT(IN), DIMENSION(num_time+1,11) ::  element_time

   ! Outputs
  REAL(KIND(1.d0)), INTENT(OUT) :: mass_csp, mrem_csp, mdwarf_csp, lbol_csp
  REAL(KIND(1.d0)), INTENT(OUT) :: SNIa_csp
  REAL(KIND(1.d0)), INTENT(OUT) :: Metal_gas_csp, Mgas_csp, Metal_star_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(11) :: element_gas_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(11) :: element_star_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(nspec) :: spec_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(5) :: n_Ostars
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(5) :: n_Bstars
  
  
  ! Internal Variables
  INTEGER :: i, s, k, imin, imax, i_young_stars
  INTEGER, INTENT(IN)::jj  
  REAL(KIND(1.d0)) :: age, t_start
  REAL(KIND(1.d0)) :: tsfr, t1, t2
  REAL(KIND(1.d0)) :: mass_frac, sfr, Mass_norm, dt
  REAL(KIND(1.d0)) :: Mgas_ini, Metal_gas_ini
  REAL(KIND(1.d0)) :: Mgas, Metal_gas, Metal_star, Met_SNIa_csp, mass_SNIa_csp, Norm_star
  REAL(KIND(1.d0)), PARAMETER :: t_min_SNIa = 0.04
  
  
  REAL(KIND(1.d0)), DIMENSION(nspec) :: spec1, spec2
  REAL(KIND(1.d0)), DIMENSION(num_time) :: mass_ssp,lbol_ssp,mrem_ssp,mdwarf_ssp,Metal_ssp,DTD,pnorm
    REAL(KIND(1.d0)), DIMENSION(num_time,11) :: element_ssp
  REAL(KIND(1.d0)), DIMENSION(num_time, 5) :: N_O, N_B
  REAL(KIND(1.d0)), DIMENSION(nspec, num_time) :: spec_ssp
  
  REAL(KIND(1.d0)), DIMENSION(11) :: element_gas, element_star
  INTEGER, DIMENSION(num_time) :: zmet_new
  REAL, DIMENSION(num_time+1,11) :: element_new
  
  
 ! Local SNIa yield variables
  character(len=2), parameter :: element_names(18) = ['Ca', 'C ', 'Fe', 'He', 'H ', 'Mg', 'Ne',&
                                'N ', 'O ', 'S ', 'Si', 'Na', 'Al', 'Ar', 'Ti', 'Cr', 'Mn', 'Ni']
  INTEGER, PARAMETER :: i_He=4, i_H=5                              
  real(8), dimension(18) :: SNIa_yields 
  real(8), dimension(18) :: elements_SNIa_csp
   
  
  ! ------------------------------------------------
  
  ! Initialize Outputs
  mass_csp = 0.0
  mrem_csp = 0.0
  mdwarf_csp = 0.0
  lbol_csp = 0.0
  Metal_gas_csp = 0.0
  Metal_star_csp = 0.0
  Mgas_csp = 0.0 
  Mgas =0.0
  Norm_star =0.0      
  spec_csp = 0.0
  spec1 = 0.0       ! spectra of young stars
  spec2 = 0.0       ! spectra of old stars
  Mass_norm = 0.0
  n_Ostars = 0.0
  n_Bstars = 0.0
  SNIa_csp = 0.0 
  !DTD = 0.0
  mass_SNIa_csp = 0.0
  Met_SNIa_csp = 0.0
  element_gas_csp = 0.0
  element_star_csp = 0.0
  
  SNIa_yields = 0.0
  elements_SNIa_csp = 0.0
  



  ! Determine start time based on star formation history
  t_start =0.0
  IF ((sfh_type .EQ. 2) .OR.(sfh_type .EQ. 3) .OR. (sfh_type .EQ. 4) .OR. (sfh_type .EQ. 1)) THEN
     t_start = Tstart * 1e9
  ENDIF


  ! Convert units and compute galaxy age
  t_tage = tage * 1e9 - t_start

  IF ((Ttrunc .LT. small_value) .OR. (Ttrunc .LT. Tstart)) THEN
     t_trunc = 10.**(time_full(num_time)) - t_start
  ELSE
     t_trunc = Ttrunc * 1e9 - t_start
  ENDIF


  ! how long ago the galaxy was quenched
  IF ((t_trunc .LE. 0.0) .OR. (t_trunc.gt.t_tage)) THEN
     t_tquench = 0.
  ELSE
     t_tquench = t_tage - t_trunc
  ENDIF

!------------------------------------------------


  ! Find indices corresponding to time intervals
  imax = MIN(MAX(locate(time_full, LOG10(t_tage))+1, 1), num_time)

  imin = MIN(MAX(locate(time_full, LOG10(t_tquench))+1, 1), num_time)


  IF (sfh_type .EQ. 0) THEN
    imax = MIN(MAX(locate(time_full, LOG10(t_tage)) + 1, 1), num_time)
    imin = MAX(imax-1 , 1)

  ENDIF

  WRITE(*,*)imin,imax


  i_young_stars = locate(time_full, logt_young_stars)

  ! Initial gas and metallicity conditions
  Mgas_ini = 1.0 / f_star
  Metal_gas_ini = Mgas_ini * z_isoc(zmet_ini)


  Mgas= Mgas_ini
  Metal_gas = Metal_gas_ini
  Metal_star = 0.0
  
  element_gas = 0.0
  element_star =0.0
  
  element_gas(i_He) = Mgas_ini *0.25
  element_gas(i_H) = Mgas_ini *0.75
  
  ! Reference: Seitenzahl et al. 2013, MNRAS, 429, 1156
  ! Below adopt the mean value of all the model results in their table 2
  IF (yield_reference_SNIa .EQ. 1) THEN
    SNIa_yields = [ 0.012, 0.0073, 0.68935, 1e-4, 0.00, 0.01928, 0.0057, 1.0e-6, 0.11,0.0935,0.248, &
                    6.8288e-5, 0.000785, 0.0148, 2.5535e-4, 0.0072,0.0106,0.065 ]
                    
  ! Reference: Iwamoto1999   https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
  ! Below adopt the mean value of all models (W, WDD, CDD) in their table 3                  
  ELSEIF (yield_reference_SNIa .EQ. 2) THEN  
     SNIa_yields = [  0.0228, 0.0508, 0.6747, 1e-4, 0.00, 0.00727, 0.00229, 1.0e-6, 0.091, 0.0914, 0.201, &
                      6.8288e-5, 3.7214e-4, 0.0191, 5.3057e-4, 0.00773, 0.00666, 0.0834 ]
                      
                      
  ! Reference: Iwamoto1999_W70    https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
  ! Below adopt the main isotope of W70 model                    
  ELSEIF (yield_reference_SNIa .EQ. 3) THEN  
     SNIa_yields = [ 0.0181, 0.0508, 0.68, 1e-4, 0.00, 0.0158, 0.00229, 1.0e-6, 0.133,  0.0914, 0.142, &
                     6.8288e-5, 1.31e-4, 0.0191, 3.13e-4, 0.00773, 0.00666, 0.0834 ]  
                     
  
  ! Reference: Iwamoto1999_W7    https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
  ! Below adopt the main isotope of W70 model                   
  ELSEIF (yield_reference_SNIa .EQ. 4) THEN  
     SNIa_yields = [ 0.0119, 0.0483, 0.626, 1e-4, 0.00, 0.0085 * 5, 0.00202, 1.0e-6, 0.143, 0.0846, 0.154, &
                     6.8288e-5, 9.86e-4, 0.0147, 2.05e-4, 0.00636,  0.00887, 0.11 ]  
  
  
  ! Reference: Iwamoto1999_WDD3    https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract                   
  ELSEIF (yield_reference_SNIa .EQ. 5) THEN  
     SNIa_yields = [ 1.88e-2, 1.66e-2, 7.95e-1, 1e-4, 0.00, 2.62e-3, 4.55e-4, 1.0e-6, 5.58e-2, 9.37e-2, 1.58e-1, &
                     3.01e-05, 1.41e-4, 1.87e-2, 5.23e-4, 1.13e-2, 6.16e-3, 4.97e-2 ]                                   
     
                   
  endif
  
  
  


  ! Main loop over time bins
  DO i=  imin, imax
     
     ! Compute time intervals
     IF(i.EQ.1)THEN
        t2 = t_tage
     ELSE
        t2 = t_tage - 10.**(time_full(i-1))
     ENDIF
     t1 = t_tage - 10.**(time_full(i))

    
     IF(t2 .LT. small_value) t2 = 0.0
     IF(t1 .LT. small_value) t1 = 0.0


     CALL sfh_gal(t1,t2, mass_frac,tsfr,dt)

     age = t_tage*10.0**(-9)

     ! Update metallicity
     k=jj+1-i

     zmet_new(i)= zmet_time(k)
     IF(k==0)THEN
        zmet_new(i)= zmet_time(1)
     ELSEIF(K.LT.jj)THEN
        zmet_new(i)= (zmet_time(k)+zmet_time(k+1))/2.       
     ENDIF


     IF (Z_MODE==0) THEN
        zmet = zmet_ini
     ELSEIF (Z_MODE==1 ) THEN
        zmet = zmet_new(i)
     ENDIF
     
     
     element_new(i,:)=0.0
     DO s=1,11
        element_new(i,s)= element_time(k,s) 
        IF(k.LT.1)THEN
          element_new(i,s)= element_time(1,s) 
        ENDIF
     ENDDO



     ! Generate SSP properties
     CALL SSP_GEN(i,age,mass_ssp,mrem_ssp,mdwarf_ssp,lbol_ssp,spec_ssp,Metal_ssp,element_ssp,N_O,N_B,pnorm)

     IF (sfh_type .EQ. 0) THEN
        mass_frac = 1./2.0
     ENDIF

     IF(t_tage .LT. 0.0) THEN
         mass_frac =0.0
     ENDIF

    !spec1 is spectra of young stars and spec2 is spectra of old stars

     IF(dust_atten.EQ.1) THEN

        IF(i .LE. i_young_stars) THEN
           spec1 = spec1 + spec_ssp(:, i)* mass_frac
        ENDIF
        IF(i.GT.i_young_stars)THEN
           spec2 = spec2 + spec_ssp(:, i)* mass_frac
        ENDIF
     ELSE
        spec_csp  = spec_csp  +  spec_ssp(:, i) * mass_frac
     ENDIF

     mass_csp = mass_ssp(i) * mass_frac + mass_csp
     mrem_csp = mrem_ssp(i) * mass_frac + mrem_csp
     mdwarf_csp = mdwarf_ssp(i) * mass_frac + mdwarf_csp
     lbol_csp = 10.0**lbol_ssp(i) * mass_frac + lbol_csp

     Mass_norm = Mass_norm + mass_frac

     n_Ostars(:) = N_O(i,:) * mass_frac + n_Ostars(:)
     n_Bstars(:) = N_B(i,:) * mass_frac + n_Bstars(:)
     
     
     !-------------------------calculate  number of SNIa--------------------------

    IF(10.**(time_full(i)-9.) .GT. t_min_SNIa ) THEN
       DTD(i) = (4.08)* log(10.**(time_full(i)-9.)/0.04)  !*10**(-4)
    ELSEIF(10.**(time_full(i)-9.) .LE. t_min_SNIa)THEN 
       DTD(i) = 0.0
    ENDIF
         
    
     
    ! Update total SNIa contribution
    ! total number of SNIa formed
    SNIa_csp = DTD(i)* mass_frac* pnorm(i) + SNIa_csp  !*10**(-4)


    !Elemental SNIa contribution
    DO s = 1, 18
       elements_SNIa_csp(s) = SNIa_csp*SNIa_yields(s)*1e-4
    ENDDO
    
   
   
    !SNIa_csp = SNIa_csp*1e-4
    Met_SNIa_csp = SNIa_csp*0.8*1e-4
    mass_SNIa_csp = SNIa_csp*1.4*1e-4


   !-----calculate metallicity------------

     ! Update gas and metallicity
     Mgas = Mgas - (mass_ssp(i) + mrem_ssp(i)) * mass_frac

     Metal_gas = Metal_gas- z_isoc(zmet_new(i))*&
                 mass_frac + Metal_ssp(i)* mass_frac
                 
     element_gas(:) = element_gas(:) -element_new(i,:)* &
                      mass_frac + element_ssp(i,:)* mass_frac            
                 
                 
     Mgas_csp  =  Mgas + mass_SNIa_csp
     
     !IF (Mgas_csp .LT. tiny_value) Mgas_csp = tiny_value
     Metal_gas_csp  = (Metal_gas+ Met_SNIa_csp)/Mgas_csp
    
    
     Metal_star =   Metal_star + z_isoc(zmet_new(i))*mass_ssp(i)* mass_frac !z_isoc(zmet_new(imin))
     Norm_star =  Norm_star +  mass_ssp(i)* mass_frac          
                   
     DO s=1,11
        !element_gas_csp(s) = (element_gas(s)+elements_SNIa_csp(s)) /Mgas_csp    
        element_star(s) =   element_star(s)+  mass_ssp(i)* mass_frac*element_new(i,s)   !element_new(imin,s)
     ENDDO               




 ENDDO
 
 Metal_star_csp = ( Metal_star ) /Norm_star
  

     DO s=1,11
        element_gas_csp(s) = (element_gas(s)+elements_SNIa_csp(s)) /Mgas_csp    
        !element_star(s) = element_star(s)+  mass_ssp(i)* mass_frac*element_gas_csp(s)   
     ENDDO
      
    
    DO s=1,11
        element_star_csp(s) = element_star(s) / Norm_star
    ENDDO
 
  !-----------------------------------------------------------------------



  IF (dust_atten .EQ. 1) THEN
     spec_csp = spec1 * EXP(-dust_depth1*(lambda/v_lambda)**dust_power) + &
                spec2 * EXP(-dust_depth2*(lambda/v_lambda)**dust_power)
  ENDIF

  ! Normalize outputs
  IF(Mass_norm .LE. 0.0)THEN
     Mass_norm=1.
  ENDIF

     mass_csp = mass_csp
     mrem_csp = mrem_csp
     mdwarf_csp = mdwarf_csp
     lbol_csp = LOG10(lbol_csp)
     spec_csp =  spec_csp



END SUBROUTINE csp_gen


