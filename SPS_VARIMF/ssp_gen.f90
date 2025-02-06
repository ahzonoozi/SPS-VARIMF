SUBROUTINE SSP_GEN(l,age,mass_ssp,mrem_ssp,mdwarf_ssp,lbol_ssp,spec_ssp,Metal_ssp,N_O,N_B)

! calculate the evolution of a single stellar population


  USE params
  
  IMPLICIT NONE
  
  INTEGER :: i, j, k, l
  REAL(KIND(1.d0)), DIMENSION(nm) :: N_IMF
  REAL(KIND(1.d0)), DIMENSION(num_rem) :: Mret,Metal_eject,Nrem
  REAL(KIND(1.d0)), INTENT(inout), DIMENSION(nspec,num_time) :: spec_ssp
  REAL(KIND(1.d0)), DIMENSION(nspec,num_time) :: tspec_ssp
  REAL(KIND(1.d0)), INTENT(inout), DIMENSION(num_time) :: mass_ssp, lbol_ssp, mrem_ssp, mdwarf_ssp,Metal_ssp
  REAL(KIND(1.d0)), DIMENSION(num_time,nm) :: mini,mact,logl,logt,logg,ffco,phase,lmdot
  INTEGER, DIMENSION(num_time) :: nmass
  REAL(KIND(1.d0)), DIMENSION(num_time) :: time
  REAL(KIND(1.d0)), DIMENSION(num_time,5) :: N_O,N_B
  REAL(KIND(1.d0)), DIMENSION(nspec) :: tspec
  REAL(KIND(1.d0)) ::time1,t1,t2
  REAL(KIND(1.d0)) :: age
  

  !-----------------------------------------------------------!
  !--------------------------Setup----------------------------!
  !-----------------------------------------------------------!

    
  !reset arrays
  spec_ssp = 0.
  mass_ssp = 0.
  mrem_ssp = 0.
  Metal_ssp = 0.
  mdwarf_ssp = 0.
  lbol_ssp = 0.

  N_O = 0.0
  N_B = 0.0

  !Check metallicity range
  IF (zmet.LT.1.OR.zmet.GT.nz) THEN
     WRITE(*,*) 'SSP_GEN ERROR: metallicity outside of range',zmet
     STOP
  ENDIF

     
     !transfer isochrones into temporary arrays
     mini  = mini_isoc(zmet,:,:)  !initial mass
     mact  = mact_isoc(zmet,:,:)  !actual (present) mass
     logl  = logl_isoc(zmet,:,:)  !log(Lbol)
     logt  = logt_isoc(zmet,:,:)  !log(Teff)
     logg  = logg_isoc(zmet,:,:)  !log(g)
     ffco  = ffco_isoc(zmet,:,:)  !C-rich or O-rich flag for TP-AGB stars
     phase = phase_isoc(zmet,:,:) !phase of evolution
     nmass = nmass_isoc(zmet,:)   !number of elements per isochrone
     time  = timestep_isoc(zmet,:)!age of each isochrone in log(yr)
     
    
     
     !-----------------------------------------------------------!
     !---------------------Generate SSPs-------------------------!
     !-----------------------------------------------------------!


  !selected isochrone
  DO i=l,l
     mass_ssp(i) = 0.0
     mrem_ssp(i) = 0.0
     Metal_ssp(i) = 0.0
      
     ! Compute time intervals 
     IF(i.EQ.1)THEN
        !time1= ( age - Tstart) * 1e9
        t2 = ( age ) * 1e9
     ELSE
        t2 = ( age ) * 1e9- 10.**(time_full(i-1))+small_value 
     ENDIF
        t1 = ( age) * 1e9 - 10.**(time_full(i))+small_value


     IF(t2.LT. small_value) t2 = small_value
     IF(t1.LT. small_value) t1 = 0.0

     IF(imf_type.EQ.2)THEN
        CALL IMF_NORM(mini(i,:),mact(i,:),nmass(i),N_IMF(:),Mret(:),Metal_eject,Nrem(:))
     ELSEIF(imf_type.EQ.5)THEN
        CALL IMF_TOPHEAVY(mini(i,:),nmass(i),N_IMF(:),Mret(:),Metal_eject(:),Nrem(:))
     ELSEIF(imf_type.EQ.3)THEN
        CALL IGIMF_NORM(t1,t2,mini(i,:),nmass(i),N_IMF(:),Mret(:),Metal_eject(:),Nrem(:))
     ENDIF

        
        !compute IMF-weighted mass of the SSP
        mass_ssp(i) = SUM(N_IMF(1:nmass(i))*mact(i,1:nmass(i)))
        mrem_ssp(i) = SUM(Nrem(1:num_rem)*Mret(1:num_rem))
        Metal_ssp(i) = SUM(Nrem(1:num_rem)*Metal_eject(1:num_rem))


!-----------------------------
 ! Calculate the number of OB stars with different masses

     IF(OB_stars==1)THEN

        Do j=1,nmass(i)
           IF (mini(i,j).GE.15.0.AND.mini(i,j).LE.30.0) THEN
              N_O(i,1) = N_IMF(j)+N_O(i,1)
           ENDIF
           IF (mini(i,j).GE.30.0.AND.mini(i,j).LE.50.0) THEN
              N_O(i,2) = N_IMF(j)+N_O(i,2)
           ENDIF
           IF (mini(i,j).GE.50.0.AND.mini(i,j).LE.70.0) THEN
              N_O(i,3) = N_IMF(j)+N_O(i,3)
           ENDIF
           IF (mini(i,j).GE.70.0.AND.mini(i,j).LE.90.0) THEN
              N_O(i,4) = N_IMF(j)+N_O(i,4)
           ENDIF
           IF (mini(i,j).GE.90.0.AND.mini(i,j).LE.150.0) THEN
              N_O(i,5) = N_IMF(j)+N_O(i,5)
           ENDIF
  !-------------------B type stars-------------------

          IF (mini(i,j).GE.2.0.AND.mini(i,j).LE.5.0) THEN
             N_B(i,1) = N_IMF(j)+N_B(i,1)
          ENDIF
          IF (mini(i,j).GE.5.0.AND.mini(i,j).LE.7.0) THEN
            N_B(i,2) = N_IMF(j)+N_B(i,2)
          ENDIF
          IF (mini(i,j).GE.7.0.AND.mini(i,j).LE.9.0) THEN
            N_B(i,3) = N_IMF(j)+N_B(i,3)
          ENDIF
          IF (mini(i,j).GT.9.AND.mini(i,j).LE.11.0) THEN
            N_B(i,4) = N_IMF(j)+N_B(i,4)
          ENDIF
          IF (mini(i,j).GE.11.0.AND.mini(i,j).LE.15.0) THEN
            N_B(i,5) = N_IMF(j)+N_B(i,5)
          ENDIF

       ENDDO
     ENDIF
!-----------------------------------
!Calculate mass of white dwarf!

     DO k=1,num_rem !2000      
        IF(Mret(k) .LT. 1.13 .AND. Mret(k) .GT. 0.08) THEN
            mdwarf_ssp(i) =  mdwarf_ssp(i) + Nrem(k)*Mret(k)
        ENDIF
     ENDDO


!---------------------------
    !Compute bolometric luminosity    
    lbol_ssp(i) = LOG10(SUM(N_IMF(1:nmass(i))*10**logl(i,1:nmass(i))))

    !compute SSP spectrum
    spec_ssp(:,i) = 0.0


    DO j=1,nmass(i)
       CALL GETSPEC(mact(i,j),logt(i,j),10.**logl(i,j),logg(i,j),phase(i,j),&
                  ffco(i,j),tspec)
       spec_ssp(:,i) = N_IMF(j)*tspec + spec_ssp(:,i)
        !spectra (in each lambda :) of a simple population in different time(:)
     ENDDO
  ENDDO

END SUBROUTINE SSP_GEN


