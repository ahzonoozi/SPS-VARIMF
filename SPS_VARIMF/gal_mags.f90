SUBROUTINE GAL_MAGS(spec,mags)

  ! Routine to calculate magnitudes in the AB system from an input spectrum.


  USE params; USE sub_func
  IMPLICIT NONE


  ! Input spectrum
  REAL(KIND(1.d0)), INTENT(INOUT), DIMENSION(nspec) :: spec  
  ! Output magnitudes
  REAL(KIND(1.d0)), INTENT(INOUT), DIMENSION(n_filtrs) :: mags  

  ! Local variables
  INTEGER :: i
  ! Spectrum flux in units of Fnu
  REAL(KIND(1.d0)), DIMENSION(nspec) :: fnu_zred  
  ! Integral of filter response
  REAL(KIND(1.d0)) :: sum_filtr      
  REAL(KIND(1.d0)), PARAMETER :: pi_value = 3.14159265
  ! Parsec in cm
  REAL(KIND(1.d0)), PARAMETER :: parsec = 3.08568E18  
  ! Solar luminosity in erg/s
  REAL(KIND(1.d0)), PARAMETER :: lsun = 3.839E33 
   


  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  ! Initialize output magnitudes to a default large value
  mags = 99.

  fnu_zred = spec

  ! Loop over all filters to calculate magnitudes
  DO i=1,n_filtrs

     sum_filtr = 0.0

     ! Compute the integral of the filter response over wavelength
     sum_filtr = SUM(ABS(lambda(2:nspec) - lambda(1:nspec-1)) * &
                     (Slm_filtr(2:nspec, i) / lambda(2:nspec) + &
                      Slm_filtr(1:nspec-1, i) / lambda(1:nspec-1))/2)
                      
     ! Avoid division by zero
     IF (sum_filtr .LT. small_value) sum_filtr = 1.0


     !compute f_mean = mean flux Fnu through filter (require for AB mags)
     !return mean flux, Fnu in  each filter(i th  filter)  /HZ
     !f_nu =int(dnu*Fnu*Snu/nu)/int(dnu*Snu/nu)
     !f_nu =int(dlm*Flm*Slm*lm/c)/int(dlm*Slm/lm)
     !f_nu =int(dlm*Fnu*Slm/lm)/int(dlm*Slm/lm)


     magsun(i) = SUM( ABS(lambda(2:nspec)-lambda(1:nspec-1))* &
                    (sun_spec(2:nspec)* Slm_filtr(2:nspec,i)/lambda(2:nspec)+ &
                     sun_spec(1:nspec-1)* Slm_filtr(1:nspec-1,i)/lambda(1:nspec-1))/2.0)

     magsun(i) = magsun(i) / sum_filtr


     ! Compute the flux of the input spectrum through the filter
     mags(i) = SUM( ABS(lambda(2:nspec)-lambda(1:nspec-1))* &
                (fnu_zred(2:nspec) * Slm_filtr(2:nspec,i)/lambda(2:nspec)+ &
                fnu_zred(1:nspec-1)* Slm_filtr(1:nspec-1,i)/lambda(1:nspec-1))/2.0)

     mags(i) = mags(i) / sum_filtr

     ! Convert solar flux to magnitude
     IF (magsun(i).LT.small_value) THEN
        magsun(i) = 99.0
     ELSE
        magsun(i) = -2.5*LOG10(magsun(i)) - 48.60
     ENDIF

     ! Convert flux to magnitude
     IF (mags(i) .LE. small_value) THEN
        mags(i) = 99.0 
     ELSE

        !m_AB = -2.5log10(f_nu)-48.60
        !Convert f_nu [Lsun/Hz] to [erg/s/HZ] by: *L_sun [erg/s]
        !Convert mags to absolute magnitude

        mags(i) = -2.5*LOG10(mags(i))-2.5*LOG10(lsun) - 48.60
        mags(i) = mags(i) + 2.5*LOG10(4.0*pi_value*100.0*parsec**2.0)

     ENDIF
  ENDDO

  
END SUBROUTINE GAL_MAGS
