SUBROUTINE sfh_gal(t1,t2, mfrac, sfr,dt)

  !----------------------------------------------------------------!
  ! Subroutine to calculate star formation history (SFH) of a galaxy.
  ! 
  ! Inputs:
  !   - t1, t2: Time interval (in years)
  ! Outputs:
  !   - mfrac: Fraction of mass formed in the time interval
  !   - sfr: Star formation rate in the time interval
  !   - dt: Duration of the time interval
  !----------------------------------------------------------------!
  
  !   - Tmax_gal is the maximum age of the galaxy or maximum isochrone age
  !   - Tcut_gal is the age of the galaxy that SF is truncated
  !   - Tprime is the given `age`,

  USE params
  IMPLICIT NONE


  REAL(KIND(1.d0)), INTENT(IN) :: t1, t2 
  REAL(KIND(1.d0)), INTENT(OUT) :: mfrac, sfr
  REAL(KIND(1.d0)) :: Tmax_gal, Tprime, Tcut_gal,dt
   REAL(KIND(1.d0)) :: tau_yr, const, total_mass, mass_tau, sfr_tau
  REAL(KIND(1.d0)) :: time1, time2

  REAL(KIND(1.d0)) ::  sfr_const, mass_const
  

 !------------------------------------------------------
  ! Initialize variables
  tau_yr = tau * 10.0**9.   ! Convert tau to years
  mfrac = 1.0
  sfr = 0.0
  
 !---------------------------------------------------------------!
  ! Handle SFH type 0: Constant star formation for 10 Myr
  
  IF (sfh_type.EQ.0) THEN
     mfrac = 1.0
     sfr = 1./(10.**7.)   !It is assumed that SF happens in 10 Myr

  ELSE IF ((sfh_type.EQ.1).OR.(sfh_type.EQ.2).OR.(sfh_type.EQ.4)) THEN

     Tmax_gal = (10.**(time_full(num_time) - 9.) - Tstart) *10.**9.
     Tprime = t2

     !Determine truncation age (Tcut_gal) 
     !Ttrunc=0.0 means that there is no truncation
     !Truncation should happen after star formation start time
     IF((Ttrunc .LT. small_value).OR.(Ttrunc .LT. Tstart))THEN
        Tcut_gal = Tmax_gal
     ELSE
        Tcut_gal = (Ttrunc - Tstart)* 10.**9.
     ENDIF


     ! Adjust time intervals based on truncation
     !age is in year
     time1 = MAX(t1, 0.0) 
     
     time2 = t2
     IF (time2 .GT. Tcut_gal)THEN
        time2 = Tcut_gal
     ENDIF
     IF (time2.LT.0.0) time2 = small_value !!!!
     !time2 = MIN(MAX(t2, 0.0), Tcut_gal)
     
     dt = time2-time1

     ! The mass-frac is the fraction of mass  that formed in the time interval time1-time2
     ! mass_formed (time1-time2) / mass_formed (maxtime)
     IF(sfh_type.EQ.2)THEN
        total_mass = tau_yr * (1. - exp(-1 * min(Tmax_gal,Tcut_gal) / tau_yr))

     ELSEIF(sfh_type.EQ.4)THEN
        total_mass = tau_yr * (1.-exp(-1 * min(Tmax_gal,Tcut_gal) / tau_yr) -&
            (min(Tmax_gal,Tcut_gal)/tau_yr) * exp(-1 * min(Tmax_gal,Tcut_gal) / tau_yr))
     ENDIF


     ! Compute mass and SFR within the time interval
     IF (Tprime.LT.0) THEN
        mass_tau = 0.0
        sfr_tau = 0.0
     ELSEIF (Tprime .GT. Tcut_gal) THEN
        mass_tau = 0.
        sfr_tau =0.
     ELSE
        ! The fraction of mass formed from time1 up to time2
        IF(sfh_type.EQ.2) THEN
           mass_tau = tau_yr * (exp(-1*time1/tau_yr)- exp(-1*time2/tau_yr))
           sfr_tau = exp(-min(Tprime, Tcut_gal) / tau_yr)
        ENDIF
        IF(sfh_type.eq.4)THEN
             mass_tau = tau_yr * ( exp(-1*time1/tau_yr)- exp(-1*time2/tau_yr)) +&
                   time1* exp(-1*time1/tau_yr) - time2* exp(-1*time2/tau_yr)
                   
             sfr_tau = (min(Tprime, Tcut_gal) / tau_yr) * &
                    exp(-min(Tprime, Tcut_gal) / tau_yr)
        ENDIF
     ENDIF


  ENDIF

  ! Handle constant star formation (sfh_type == 1)
  IF(sfh_type.EQ.1)THEN
     const =1.0
     total_mass = 1.0
  ELSE
     const=0.0
  ENDIF


  !Calculate mass fraction (mfrac) and SFR (sfr)
  IF ((sfh_type.EQ.1).OR.(sfh_type.EQ.2).OR.(sfh_type.EQ.4)) THEN

     IF (Tprime.gt.0) THEN
        sfr_const = 1.0 / MIN(Tmax_gal, Tcut_gal)
        mass_const = (time2-time1) /  MIN(Tmax_gal, Tcut_gal)
     ELSE
        sfr_const = 0.0
        mass_const = 0.0 
     ENDIF

     mfrac = (1. - const) * mass_tau / total_mass + const *  mass_const

     IF (Tprime.gt.Tcut_gal) THEN
        sfr = 0.0
     ELSE
        sfr = (1. - const) * sfr_tau / total_mass + const * sfr_const
     ENDIF
  
  ENDIF

  ! SFR  [M_sun/yr] . Total mass is normalized to 1 M_sun
  sfr=sfr 

  
END SUBROUTINE sfh_gal


