PROGRAM MAIN

  USE params; USE sub_func
  IMPLICIT NONE

  
  INTEGER :: i 
  !SSP spectrum, stellar Mass, Lbol, remnant and dwarf mass
  REAL(KIND(1.d0)), DIMENSION(num_time,nspec)  :: spec_ssp
  REAL(KIND(1.d0)), DIMENSION(num_time) :: mass_ssp,lbol_ssp,mrem_ssp,mdwarf_ssp,Metal_ssp
  CHARACTER(100) :: file_name=''
  CHARACTER(len=:), ALLOCATABLE :: IMF_name
  CHARACTER(len=:), ALLOCATABLE :: SFR_name
  !CHARACTER(len=:), ALLOCATABLE :: tau_value

  
  !----------------------------------------
  !----------------------------------------

  isoc_type = 2
  !Isochrone:   1= padova2007 isochrone;
  !             2=parsec2022 isochrone

  spec_type =1
  !Spectra:     1= BaSel;
  !             2=Miles  !in progress, not available yet

  imf_type  = 3
  !IMF: =1 is reserved for salpeter IMF; But is not working yet
  !     =2 canonical IMF;
  !     =3 IGIMF ;
  !     =5 Top_heavy IMF for star clusters or UCDS that form monothically;


  sfh_type  = 1
  !SFH: sfh=0  simple stellar population,
  !     sfh=1  constant SFR  model;
  !     sfh=2  exponentially declining  model;
  !     sfh=4 delay-tau model;



  !e_folding time_scale. used in sfh_type=2 and 4
  !the range is between 0.001 -1000.0
  tau= 1000.0 !0.1


  !The start and truncation time of star formation
  Tstart = 0.0
  Ttrunc = 0.0


  !zmet_ini is the initial metallicity of galaxy
  !remain constant during the gallactic evolution
  !or change with time
  zmet_ini = 11  !22

  !The mass of all stars ever formed
  M_galaxy = 1.*10.**(11.)
  M_UCD = 1.*10.**(8.)
  !star transformation fraction, f_star=M_galaxy/M_gas
  f_star = 0.3



  ! zsun_isoc=0.019  for PADOVA isochron;  =0.0154 for PARSEC isochrone
  IF (isoc_type == 1) THEN
     zsun_isoc = 0.019
  ELSEIF (isoc_type == 2) THEN
     zsun_isoc = 0.0154
  ENDIF


  CALL READ_INPUT(zmet_ini) !read in the isochrones and spectral libraries


  IF(imf_type  == 2)THEN
    alpha1=1.3
    alpha2=2.3
  ENDIF


  !Note: alpha1 and alpha2 can only be defined here if the metallicity is time-independent.   
  !For time-evolving metallicity, alpha1 and alpha2 are dynamically assigned in the IMF subroutine.
     
  IF(imf_type  == 3 .OR. imf_type == 5)THEN
    !alpha1=1.3+63*(z_isoc(zmet_ini)-0.0142)
    !alpha2=2.3+63*(z_isoc(zmet_ini)-0.0142)

    alpha1=1.3+ 79.4*(z_isoc(zmet_ini)-0.8*0.0142)
    alpha2=2.3+ 79.4*(z_isoc(zmet_ini)-0.8*0.0142)
  ENDIF



 ! WRITE(tau_value,'(F8.3)') tau

     IF(imf_type==2)THEN
        ALLOCATE(CHARACTER(len=9) :: IMF_name)
        IMF_name = 'KroupaIMF'
     ELSEIF(imf_type==3)THEN
        ALLOCATE(CHARACTER(len=5) :: IMF_name)
        IMF_name = 'IGIMF'
     ELSEIF(imf_type == 5)THEN
        ALLOCATE(CHARACTER(len=8) :: IMF_name)
        IMF_name = 'Topheavy'
     END IF



     IF(sfh_type==1)THEN
        ALLOCATE(CHARACTER(len=8) :: SFR_name)
        SFR_name = 'ConstSFR'
     ELSEIF(sfh_type==2)THEN
       ALLOCATE(CHARACTER(len=11) :: SFR_name)
        SFR_name = 'exp_decline'
     ELSEIF(sfh_type==4)THEN
        ALLOCATE(CHARACTER(len=9) :: SFR_name)
        SFR_name = 'delay_tau'
     END IF




  IF(sfh_type==0)THEN
     file_name = 'SSP_'//IMF_name//'_'
  ELSE
     file_name = 'CSP_'//IMF_name//'_'//SFR_name//'_'
  END IF


  CALL CSP_OUT(file_name)


END PROGRAM MAIN






