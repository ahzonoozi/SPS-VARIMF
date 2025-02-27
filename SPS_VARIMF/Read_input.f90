SUBROUTINE READ_INPUT(Zinput)

  ! Subroutine to read isochrones, spectral libraries, and yields 
  
  USE params; USE sub_func
  IMPLICIT NONE


  INTEGER, INTENT(IN) :: Zinput
  INTEGER :: stat=1,i,j,m,jj,k,i_min,ll,kk
  INTEGER :: n_isoc,z,zmin,zmax
  INTEGER, PARAMETER ::  nline=1000000, num_yields = 4
  REAL(KIND(1.d0)), DIMENSION(nzinit):: z_speclib
  REAL(KIND(1.d0)) :: logage, aa, bb
  !REAL(KIND(1.d0)), DIMENSION(2901) ::mass_yield
  !REAL(KIND(1.d0)), DIMENSION(4,2901) :: Metal_yield
  REAL(KIND(1.0)), DIMENSION(nspec,nzinit,ndim_logt,ndim_logg) :: flux_speclib=0.

  CHARACTER(6) :: zstype,zyname
  CHARACTER(1) :: header

  !----------------------------------------------------------------!
  ! Initialize ranges for metallicity processing
  !----------------------------------------------------------------!
  
  IF (Z_MODE == 0) THEN
     zmin = Zinput
     zmax = Zinput
  ELSE
     zmin = 1
     zmax = nz
  ENDIF
  
  !----------------------------------------------------------------!
  ! Read Isochrones
  !----------------------------------------------------------------!

  !Reading different values of Zlegend from zlegend.dat file 
  IF (isoc_type == 1) THEN
     OPEN(100, FILE='../ISOCHRONES/Padova2007/zlegend.dat', STATUS='OLD', IOSTAT=stat)
  ELSE IF (isoc_type == 2) THEN
     OPEN(100, FILE='../ISOCHRONES/PARSEC2022/zlegend.dat', STATUS='OLD', IOSTAT=stat)
  ELSE
     PRINT *, "Error: Unsupported isochrone type!"
     STOP
  ENDIF
  
  ! Check file status
  IF (stat /= 0) THEN
     PRINT *, "Error: Could not open isochrone legend file!"
     STOP
  ENDIF
  
  ! Read metallicity values
  DO z = 1, nz
     READ(100, '(F6.4)', IOSTAT=stat) z_isoc(z)
     IF (stat /= 0) EXIT
  ENDDO
  CLOSE(100)
  
  !----------------------------------------------------------------!

  ! Process each metallicity
  DO z = zmin,zmax
     n_isoc = 0
     WRITE(zstype,'(F6.4)') z_isoc(z)

     !Open corresponding isochrone file 
     !open Padova isochrones
     iF (isoc_type == 1)Then
        OPEN(101,FILE='../ISOCHRONES/Padova2007/isoc_z'//&
             zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     END IF

     !open Parsec isochrones
     iF (isoc_type == 2)Then
        OPEN(101,FILE='../ISOCHRONES/PARSEC2022/isoc_z'//&
             zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     END IF
     ! Check file status
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open isochrone file for Z=", zstype
        STOP
     ENDIF

     ! Read isochrone data
     READ(101,*,IOSTAT=stat)
     DO i=1, num_isoc
        m=0
        n_isoc= i

        DO j=1,nline
           m=m+1
           READ(101,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                   mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                   logt_isoc(z,n_isoc,m),logg_isoc(z,n_isoc,m),&
                   ffco_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m)
           
           IF (stat.NE.0) GOTO 40

           timestep_isoc(z,n_isoc) = logage
           nmass_isoc(z,n_isoc) = nmass_isoc(z,n_isoc)+1 

        ENDDO
   
40   CONTINUE

     ENDDO
     CLOSE(101)
  ENDDO
  lambda   = 0.

  !----------------------------------------------------------------!
  ! Read Stellar Yields
  !----------------------------------------------------------------!
 
  OPEN(100,FILE='../YIELD/zyield.txt'&
          ,STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat /= 0) THEN
     PRINT *, "Error: Could not open yield legend file!"
     STOP
  ENDIF

  ! Read yield metallicity values
  DO i = 1,num_yields
    READ(100,'(F6.4)', IOSTAT=stat) z_yield(i)
    IF (stat /= 0) EXIT
  ENDDO
  CLOSE(100)

  ! Process each yield metallicity
  DO i = 1,num_yields 
     WRITE(zyname,'(F6.4)') z_yield(i)

     OPEN(101,FILE='../YIELD/portinari98_Z'//&
          zyname//'.txt',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open yield file for Z=", zyname
        STOP
     ENDIF 
     
     READ(101,*) header
     DO j =1, 2901
        READ(101,*,IOSTAT=stat) mass_yield(j),Metal_yield(i,j)
     ENDDO
     CLOSE(101)
  ENDDO

 
  !----------------------------------------------------------------!
  ! Read Spectral Libraries
  !----------------------------------------------------------------!

  IF (spec_type.EQ.1)THEN
     OPEN(100,FILE='../SPECTRA/BaSeL3.1/basel.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open spectral wavelength file!"
        STOP
     ENDIF
     OPEN(120,FILE='../SPECTRA/BaSeL3.1/zlegend.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open spectral legend file!"
        STOP
     ENDIF    
  ENDIF


  ! Read wavelengths
  DO i=1,nspec
     READ(100, *, IOSTAT=stat) lambda(i)
     IF (stat /= 0) EXIT
  ENDDO
  CLOSE(100)


  !read in primary logg and logt arrays
  OPEN(101,FILE='../SPECTRA/BaSeL3.1/basel_logt.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')

  DO i=1,ndim_logt
     READ(101,*) speclib_logt(i)
  ENDDO
  CLOSE(101)

  OPEN(102,FILE='../SPECTRA/BaSeL3.1/basel_logg.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')

  DO i=1,ndim_logg
     READ(102,*) speclib_logg(i)
  ENDDO
  CLOSE(102)
  
  z_speclib=-99.


  !read in each metallicity
  DO z=1,nzinit 
     READ(120,*, IOSTAT=stat) z_speclib(z)
     WRITE(zstype,'(F6.4)') z_speclib(z)
     IF (stat /= 0) EXIT

     ! Read spectral libraries for each metallicity
     IF (spec_type.EQ.1) THEN
         OPEN(110,FILE='../SPECTRA/BaSeL3.1/basel'//&
             '_z'//zstype//'.spectra.txt',FORM='FORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ENDIF

     READ(110,*, IOSTAT=stat) flux_speclib(:,z,:,:)
     IF (stat /= 0) EXIT
     CLOSE(110)
  ENDDO
 
  CLOSE(120)
 
   ! Interpolate Spectral Library to Isochrone Grid
   DO z=1,nz

     i_min = MIN(MAX(locate(LOG10(z_speclib/zsun_spec),&
          LOG10(z_isoc(z)/zsun_isoc)),1),nzinit-1)
     aa = (LOG10(z_isoc(z)/zsun_isoc)-LOG10(z_speclib(i_min)/zsun_spec)) / &
          (LOG10(z_speclib(i_min+1)/zsun_spec)-LOG10(z_speclib(i_min)/zsun_spec))

     IF (aa .LT. 0.0)  aa=0.0
     IF (aa .GT. 1.0)  aa=1.0

     bb = 1-aa

     spec_flux(:,z,:,:) = bb*LOG10(flux_speclib(:,i_min,:,:)+small_value) + &
          aa*LOG10(flux_speclib(:,i_min+1,:,:)+small_value)
     
     spec_flux(:,z,:,:) = 10**spec_flux(:,z,:,:)
  
  ENDDO

  !----------------------------------------------------------------!
  ! Post-processing for specific spectra
  !----------------------------------------------------------------!
    CALL SPEC_WMB()
    CALL SPEC_AGB()
    CALL SPEC_PAGB()
    CALL SPEC_WR()
   !subroutine to read filters
    CALL Filtr()
  !----------------------------------------------------------------!

     DO i=1,num_time

           time_full(i) = timestep_isoc(zmin,i)

     ENDDO

  !----------------------------------------------------------------!


END SUBROUTINE READ_INPUT


