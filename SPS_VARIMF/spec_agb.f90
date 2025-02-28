SUBROUTINE SPEC_AGB()

  !Read in TP-AGB spectra Library from Lancon & Mouhcine 2002.

  USE params; USE sub_func
  IMPLICIT NONE

  INTEGER :: stat=1,i,j
  INTEGER, PARAMETER :: nspec_agb=6152

  REAL(KIND(1.d0)) :: aa,bb,z_null
  REAL(KIND(1.d0)), DIMENSION(n_agb_o):: color_o
  REAL(KIND(1.d0)), DIMENSION(n_agb_c):: color_c
  
  REAL(KIND(1.d0)), DIMENSION(22) :: agb_logz_o
  REAL(KIND(1.d0)), DIMENSION(nspec_agb) :: agb_lm=0.0
  REAL(KIND(1.d0)), DIMENSION(nspec_agb,n_agb_o) :: agb_specinit_o=0.
  REAL(KIND(1.d0)), DIMENSION(nspec_agb,n_agb_c) :: agb_specinit_c=0.

  CHARACTER(1) :: header
  CHARACTER(4) :: agb_type

  !----------------------------------------------------------------!
  !read in AGB Teff array for O-rich spectra
  OPEN(101,FILE='../SPECTRA/AGB_spectra/Orich.teff',&
       STATUS='OLD',ACTION='READ', IOSTAT=stat)
  IF (stat /= 0) THEN
     PRINT *, "Error: Could not open Orich.teff file"
     STOP
  ENDIF
       
  
  ! Skip the header
  READ(101,*) header
  
  ! Read metallicity and temperature data for O-rich spectra
  READ(101,*) z_null, agb_logz_o
  DO i=1,n_agb_o
     READ(101,*) color_o(i), agb_logt_o(i,:)
  ENDDO
  CLOSE(101)

 !Convert to logarithmic scale
  agb_logt_o = LOG10(agb_logt_o)


  Do i=1,n_agb_o
     WRITE(agb_type,'(F4.2)')  color_o(i)
     ! Read O-rich spectra for the given color
     OPEN(110,FILE='../SPECTRA/AGB_spectra/o_bin_'//agb_type//'.dat',&
         STATUS='OLD',ACTION='READ', IOSTAT=stat)
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open O-rich spectra file for color ", agb_type
        STOP
     ENDIF
     
     DO j=1,nspec_agb
        READ(110,*) agb_lm(j),agb_specinit_o(j,i)
     ENDDO
     CLOSE(110)

     !interpolate to the main spectral grid
     agb_spec_o(:,i) =  MAX(linterparr(agb_lm,agb_specinit_o(:,i),&
          lambda),small_value)
         
  ENDDO

 !---------------------------------------------
  !read in AGB Teff array for C-rich spectra
  OPEN(94,FILE='../SPECTRA/AGB_spectra/Crich.teff',&
      STATUS='OLD',ACTION='READ', IOSTAT=stat)
  IF (stat /= 0) THEN
     PRINT *, "Error: Could not open Crich.teff file"
     STOP
  ENDIF
  
  !Skip the header
  READ(94,*) header
  DO i=1,n_agb_c
     READ(94,*) color_c(i), agb_logt_c(i)
  ENDDO
  CLOSE(94)
  agb_logt_c = LOG10(agb_logt_c)


  DO i=1,n_agb_c

     WRITE(agb_type,'(F4.2)')  color_c(i)
     !read in TP-AGB C-rich spectra
     OPEN(96,FILE='../SPECTRA/AGB_spectra/c_bin_'//agb_type//'.dat',&
        STATUS='OLD',ACTION='READ', IOSTAT=stat)
     IF (stat /= 0) THEN
        PRINT *, "Error: Could not open C-rich spectra file for color ", agb_type
        STOP
     ENDIF
     
     DO j=1,nspec_agb
         READ(96,*) agb_lm(i),agb_specinit_c(j,i)
     ENDDO
     CLOSE(96)


     !interpolate to the main spectral grid
        agb_spec_c(:,i) = MAX(linterparr(agb_lm,agb_specinit_c(:,i),&
            lambda),small_value)
            
    ENDDO

  !-------------------------------------------------------------

END SUBROUTINE SPEC_AGB


