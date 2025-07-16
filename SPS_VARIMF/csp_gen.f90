SUBROUTINE csp_gen(j, zmet_time,tage, mass_csp, mrem_csp, mdwarf_csp, lbol_csp,&
                   spec_csp, n_Ostars, n_Bstars, Metal_csp, Mgas_csp)



  USE params; USE sub_func
  IMPLICIT NONE


  ! Inputs
  REAL(KIND(1.d0)), INTENT(IN) :: tage
  INTEGER , INTENT(IN), DIMENSION(num_time) ::  Zmet_time


   ! Outputs
  REAL(KIND(1.d0)), INTENT(OUT) :: mass_csp, mrem_csp, mdwarf_csp, lbol_csp
  REAL(KIND(1.d0)), INTENT(OUT) :: Metal_csp, Mgas_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(nspec) :: spec_csp
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(5) :: n_Ostars
  REAL(KIND(1.d0)), INTENT(OUT), DIMENSION(5) :: n_Bstars
  
  
  ! Internal Variables
  INTEGER :: i, j, k, imin, imax, i_young_stars
  REAL(KIND(1.d0)) :: age, t_start
  REAL(KIND(1.d0)) :: tsfr, t1, t2
  REAL(KIND(1.d0)) :: mass_frac, sfr, Mass_norm, dt
  REAL(KIND(1.d0)) :: Mgas_ini, Metal_gas_ini
  REAL(KIND(1.d0)), DIMENSION(nspec) :: spec1, spec2
  REAL(KIND(1.d0)), DIMENSION(num_time) :: mass_ssp, lbol_ssp, mrem_ssp, mdwarf_ssp, Metal_ssp
  REAL(KIND(1.d0)), DIMENSION(num_time, 5) :: N_O, N_B
  REAL(KIND(1.d0)), DIMENSION(nspec, num_time) :: spec_ssp
  REAL(KIND(1.d0)) :: Mgas, Metal_gas, metal
  INTEGER, DIMENSION(num_time) :: zmet_new
  
  
  
  
  ! ------------------------------------------------
  
  ! Initialize Outputs
  mass_csp = 0.0
  mrem_csp = 0.0
  mdwarf_csp = 0.0
  lbol_csp = 0.0
  Metal_csp = 0.0
  Mgas_csp = 0.0    
  spec_csp = 0.0
  spec1 = 0.0       ! spectra of young stars
  spec2 = 0.0       ! spectra of old stars
  Mass_norm = 0.0
  n_Ostars = 0.0
  n_Bstars = 0.0



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


  Mgas=  Mgas_ini
  Metal_gas = Metal_gas_ini


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
     k=j+1-i

     zmet_new(i)= zmet_time(k)
     IF(k==0)THEN
        zmet_new(i)= zmet_time(1)
     ENDIF


     IF (Z_MODE==0) THEN
        zmet = zmet_ini
     ELSEIF (Z_MODE==1 ) THEN
        zmet = zmet_new(i)
     ENDIF

     ! Generate SSP properties
     CALL SSP_GEN(i,age,mass_ssp,mrem_ssp,mdwarf_ssp,lbol_ssp,spec_ssp,Metal_ssp,N_O,N_B)

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
     mdwarf_csp = mdwarf_ssp(i) * mass_frac +mdwarf_csp
     lbol_csp = 10.0**lbol_ssp(i) * mass_frac + lbol_csp

     Mass_norm = Mass_norm + mass_frac

     n_Ostars(:) = N_O(i,:) * mass_frac + n_Ostars(:)
     n_Bstars(:) = N_B(i,:) * mass_frac + n_Bstars(:)


!----------------------calculate metallicity------------------------------

     ! Update gas and metallicity
     Mgas = Mgas - mass_ssp(i)* mass_frac - mrem_ssp(i)* mass_frac

     Metal_gas = Metal_gas- z_isoc(zmet_new(i))*&
                 mass_frac + Metal_ssp(i)* mass_frac


!---------------------------------------------------

 ENDDO


  ! Final values
   Metal = Metal_gas / Mgas
   Metal_csp =  Metal
   Mgas_csp  =  Mgas


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


