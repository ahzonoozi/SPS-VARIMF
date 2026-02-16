SUBROUTINE CSP_OUT(output_name)
  
  USE params; USE sub_func                       
  IMPLICIT NONE


  ! ---- Inputs ----
  CHARACTER(100), INTENT(IN) :: output_name

  ! ---- Local variables  ----
  INTEGER :: i,k
  REAL(KIND(1.d0)), DIMENSION(num_time) :: lbol_ssp, mass_ssp, mrem_ssp, mdwarf_ssp, Metal_ssp
  REAL(KIND(1.d0)), DIMENSION(nspec, num_time) :: spec_ssp
  REAL(KIND(1.d0)) :: lbol_csp, mass_csp, mrem_csp, mdwarf_csp
  REAL(KIND(1.d0)) :: Metal_gas_csp, Metal_star_csp, Mgas_csp, SNIa_csp
  REAL(KIND(1.d0)), DIMENSION(11) :: element_gas_csp
  REAL(KIND(1.d0)), DIMENSION(11) :: element_star_csp
  REAL(KIND(1.d0)) :: age, tage, mass_frac, tsfr, dt
  REAL(KIND(1.d0)), DIMENSION(nspec) :: spec_csp
  REAL(KIND(1.d0)), DIMENSION(n_filtrs) :: mags
  REAL(KIND(1.d0)), DIMENSION(5) :: n_Ostars
  REAL(KIND(1.d0)), DIMENSION(5) :: n_Bstars
  REAL(KIND(1.d0)), DIMENSION(num_time,5) :: N_O, N_B
  
  INTEGER, DIMENSION(num_time+1) ::zmet_time
  REAL, DIMENSION(num_time+1,11) :: element_time
  REAL(KIND(1.d0)), DIMENSION(n_index) :: indices
  REAL, DIMENSION(num_time+1) :: Metal_s
  REAL, DIMENSION(num_time+1,11) :: element_s

  CHARACTER(34) :: fmt
  
  ! ------------------------------------------------------------------
  
  ! Initialize output files 
  CALL SETUP_OUTPUT(output_name, 1, num_time)

  ! Initialize arrays
  spec_ssp = 0.
  mass_ssp = 0.
  lbol_ssp = 0.
  mrem_ssp = 0.
  Metal_ssp =0.
  mdwarf_ssp = 0.
  zmet_time = 0.
  element_time = 0.



 ! --- Get CSP spectra ------------------------
 
     zmet_time(1)= zmet_ini
     element_time(1,:) = 0.
     element_time(1,4) = 0.25  !Y
     element_time(1,5) = 0.75  !X
     

  !---- Loop over output ages.
  DO i=1, num_time
  
     age = 10.**(time_full(i)-9.)    ! Convert log(age) to Gyr
     WRITE(*,*)time_full(i)
     
     

     ! Reset outputs for this time step
     Metal_s(i) = 0.0
     element_s(i,:) = 0.0
      
     mass_csp=0.
     
     
     ! Generate composite stellar population properties  
     CALL csp_gen(i,zmet_time,element_time,age,mass_csp,mrem_csp,mdwarf_csp,&
                  lbol_csp,spec_csp,n_Ostars,n_Bstars,Metal_gas_csp,Metal_star_csp,&
                  element_gas_csp,element_star_csp,Mgas_csp,SNIa_csp)

     zmet_time(i+1)= MAX(locate(z_isoc,Metal_gas_csp),1)
     element_time(i+1,:) = element_gas_csp(:)
     
     ! Truncate metallicity if age > truncation time
     Metal_s(i) = Metal_star_csp
     element_s(i,:) = element_star_csp(:)
     
     IF(age .GT. Ttrunc) THEN
        Metal_s(i) = Metal_s(i-1)
        element_s(i,:)= element_s(i-1,:)
     ENDIF
     
     

     !compute indices
     IF(index_cal .EQ. 1) THEN
        CALL compute_indx(lambda,spec_csp,indices)
     ENDIF
     
     ! Update mass and calculate magnitudes 
     mass_csp  = mass_csp  + mrem_csp
     CALL gal_mags(spec_csp, mags)

     ! Calculate age and star formation history
     tage = (age - Tstart) * 10.**9
     CALL sfh_gal(0.0,tage, mass_frac, tsfr,dt)
     
     WRITE(50,*)  LOG10(age) + 9.0,Metal_gas_csp,(element_gas_csp(k),k=1,11),Metal_s(i),(element_s(i,k),k=1,11)


 !----------------------------------------------------------------
 !---- print and save outputs in .mag and .spec files ------------

    ! Format for magnitudes file
    fmt = '(F7.4,1x,9(F9.5,1x),000(F7.3,1x))'
    WRITE(fmt(21:23),'(I3,1x,I4)') n_filtrs


    ! Write magnitudes to file
    WRITE(10,fmt) LOG10(age) + 9.0, LOG10(mass_csp + small_value), &
        lbol_csp, LOG10(tsfr + small_value), LOG10(mrem_csp + small_value), &
        LOG10(mdwarf_csp + small_value), Metal_gas_csp, Metal_s(i),&
        (Mgas_csp + small_value),SNIa_csp, mags  
        
        
    ! Write spectra to file
    WRITE(20,'(4(F8.4,1x))') LOG10(age) + 9.0,&
        LOG10(mass_csp + small_value), lbol_csp, LOG10(tsfr + small_value)
    WRITE(20,'(50000(E14.6))') MAX(spec_csp,small_value)

    !write indices to file
    WRITE(30,'(57(F8.4,4x))') LOG10(age) + 9.0, indices
    
    ! Write OB star counts if enabled 
    IF(OB_stars==1)THEN
       WRITE(25,*) LOG10(age) + 9.0, n_Ostars, n_Bstars
    ENDIF

 
  ENDDO
  
  ! Close OB stars file
  CLOSE(25)
  CLOSE(50)

END SUBROUTINE CSP_OUT


!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE SETUP_OUTPUT(output_name, imin, imax)

  USE params
  IMPLICIT NONE
  
  ! Inputs
  INTEGER, INTENT(IN) :: imin, imax
  CHARACTER(100), INTENT(IN) :: output_name

  ! Local variables
  CHARACTER(34) :: fmt
  !-----------------------------------------------------!
  IF(OB_stars==1)THEN  
     OPEN(25,FILE='../OUTPUTS/N_OBstars.txt',STATUS='REPLACE')
  ENDIF
  OPEN(50,FILE='../OUTPUTS/ELEMENTS.txt',STATUS='REPLACE')
  
  !------------------------------------------------------!

3 FORMAT('#   tau/Gyr= ',F8.3,', Tstart= 'F6.3,',Ttrunc= 'F6.3)



   !! Open magnitude output file
  OPEN(10,FILE='../OUTPUTS/'//TRIM(output_name)//'.mags',&
     STATUS='REPLACE')

  ! Write IMF type to magnitudes file
  SELECT CASE (imf_type)
  CASE (2)
     WRITE(10, '("#   IMF: canonical IMF, slopes= ",3F4.1)') alpha1, alpha2
  CASE (3)
     WRITE(10, '("#   IMF: IGIMF ")')
  CASE (5)
     WRITE(10, '("#   IMF: Top heavy IMF")')
  CASE DEFAULT
     WRITE(10, '("#   IMF: ",I1)') imf_type
  END SELECT
  
  
  
  WRITE(10,'("#   Mag Zero Point: AB ")')
  IF(Z_MODE .EQ. 0) THEN
     WRITE(10,'("#   Log(Z/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ELSE
     WRITE(10,'("#   Log(Z_ini/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ENDIF

  WRITE(10,'("#   Mgalaxy: ", E7.1)')M_galaxy

  IF(Z_MODE .EQ. 1) THEN
     WRITE(10,'("#   f_star = Mgalaxy/Mgas =: ", E7.1)')f_star
  ELSEIF(Z_MODE .EQ. 0)THEN
     WRITE(10,'("#   Metallicity remains constant over time", A)')
  ENDIF



  ! Write star formation history type
  SELECT CASE (sfh_type)
  CASE (0)
     WRITE(10, '("#   SFH: SSP")')
  CASE (1)
     WRITE(10, '("#   SFH: constant SFR")')
     WRITE(10, 3) tau, Tstart, Ttrunc
  CASE (2)
     WRITE(10, '("#   SFH: exponentially declining")')
     WRITE(10, 3) tau, Tstart, Ttrunc
  CASE (3)
     WRITE(10, '("#   SFH: Eappen et al. 2022")')
     WRITE(10, 3) tpeak, dTm, Tstart, Ttrunc
  CASE (4)
     WRITE(10, '("#   SFH: delay-tau")')
     WRITE(10, 3) tau, Tstart, Ttrunc
  END SELECT
  


  ! Write column headers for magnitudes
  WRITE(10,'("#")')
  WRITE(10,'("#   log(age) log(mass) Log(lbol) log(SFR) log(mass_rem)&
            log(mass_wd)  Z_gas Z_star  Mgas_csp  num_SNIa[1e-4]  mags (see FILTER_LIST)")')


!-----------------------------------------


   !! Open spectra output file
  OPEN(20,FILE='../OUTPUTS/'//TRIM(output_name)//'.spec',&
     STATUS='REPLACE')

  
  SELECT CASE (imf_type)
  CASE (2)
     WRITE(20, '("#   IMF: canonical IMF, slopes= ",3F4.1)') alpha1, alpha2
  CASE (3)
     WRITE(20, '("#   IMF: IGIMF ")')
  CASE (5)
     WRITE(20, '("#   IMF: Top heavy IMF")')
  CASE DEFAULT
     WRITE(20, '("#   IMF: ",I1)') imf_type
  END SELECT
  
  
  
  WRITE(20,'("#   Mag Zero Point: AB ")')
  IF(Z_MODE .EQ. 0) THEN
     WRITE(20,'("#   Log(Z/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ELSE
     WRITE(20,'("#   Log(Z_ini/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ENDIF

  WRITE(20,'("#  Mgalaxy: ", E7.1)')M_galaxy

  IF(Z_MODE .EQ. 1) THEN
     WRITE(20,'("#   f_star = Mgalaxy/Mgas =: ", E7.1)')f_star
  ELSEIF(Z_MODE .EQ. 0)THEN
     WRITE(20,'("#   Metallicity remains constant over time", A)')   
  ENDIF



  ! Write star formation history type
  SELECT CASE (sfh_type)
  CASE (0)
     WRITE(20, '("#   SFH: SSP")')
  CASE (1)
     WRITE(20, '("#   SFH: constant SFR")')
     WRITE(20, 3) tau, Tstart, Ttrunc
  CASE (2)
     WRITE(20, '("#   SFH: exponentially declining")')
     WRITE(20, 3) tau, Tstart, Ttrunc
  CASE (3)
     WRITE(20, '("#   SFH: Eappen et al. 2022")')
     WRITE(20, 3) tpeak, dTm, Tstart, Ttrunc   
  CASE (4)
     WRITE(20, '("#   SFH: delay-tau")')
     WRITE(20, 3) tau, Tstart, Ttrunc
  END SELECT



  WRITE(20,'("#")')
  WRITE(20,'("#   log(age) log(mass) Log(lbol) log(SFR) spectra ")')

  IF (imax-imin.EQ.1) WRITE(20,'(I3,1x,I6)') 1,nspec
  IF (imax-imin.GT.1) WRITE(20,'(I3,1x,I6)') num_time,nspec
  WRITE(20,'(50000(F15.4))') lambda

!-----------------------------------------


  ! Open indices output file
  IF(index_cal==1)THEN
     OPEN(30,FILE='../OUTPUTS/'//TRIM(output_name)//'.indx',&
     STATUS='REPLACE')
  ENDIF
     
     
     SELECT CASE (imf_type)
  CASE (2)
     WRITE(30, '("#   IMF: canonical IMF, slopes= ",3F4.1)') alpha1, alpha2
  CASE (3)
     WRITE(30, '("#   IMF: IGIMF ")')
  CASE (5)
     WRITE(30, '("#   IMF: Top heavy IMF")')
  CASE DEFAULT
     WRITE(30, '("#   IMF: ",I1)') imf_type
  END SELECT
  
  
  
  WRITE(30,'("#   Mag Zero Point: AB ")')
  IF(Z_MODE .EQ. 0) THEN
     WRITE(30,'("#   Log(Z/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ELSE
     WRITE(30,'("#   Log(Z_ini/Zsun): ",F6.3)') LOG10(z_isoc(zmet_ini)/zsun_isoc)
  ENDIF

  WRITE(30,'("#  Mgalaxy: ", E7.1)')M_galaxy

  IF(Z_MODE .EQ. 1) THEN
     WRITE(30,'("#   f_star = Mgalaxy/Mgas =: ", E7.1)')f_star
  ELSEIF(Z_MODE .EQ. 0)THEN
     WRITE(30,'("#   Metallicity remains constant over time", A)')   
  ENDIF



  ! Write star formation history type
  SELECT CASE (sfh_type)
  CASE (0)
     WRITE(30, '("#   SFH: SSP")')
  CASE (1)
     WRITE(30, '("#   SFH: constant SFR")')
     WRITE(30, 3) tau, Tstart, Ttrunc
  CASE (2)
     WRITE(30, '("#   SFH: exponentially declining")')
     WRITE(30, 3) tau, Tstart, Ttrunc
  CASE (4)
     WRITE(30, '("#   SFH: delay-tau")')
     WRITE(20, 3) tau, Tstart, Ttrunc
  END SELECT



  WRITE(30,'("#   log(age)  indices ")')
 
     



END SUBROUTINE SETUP_OUTPUT

!------------------------------------------------------------!
