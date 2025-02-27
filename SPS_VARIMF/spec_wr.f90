SUBROUTINE SPEC_WR()

  !read in Wolf-Rayet spectra from Smith et al.


  USE params; USE sub_func
  IMPLICIT NONE



  INTEGER :: i,j,i1
  INTEGER, PARAMETER :: nlamwr=1963
  REAL(KIND(1.d0)) :: aa,bb
  REAL(KIND(1.d0)), DIMENSION(5) :: wr_z=0.
  REAL(KIND(1.d0)), DIMENSION(nlamwr,nt_wr,5) :: tspec_wrc=0.,tspec_wrn=0.
  REAL(KIND(1.d0)), DIMENSION(nlamwr) :: lm_wr=0.,i_spec_wr=0.


  !-------------------------------------------------------

    
  !read in WR-N Teff and spectra 
  OPEN(110,FILE='../SPECTRA/WR_spectra/CMFGEN_WN'//&
       '.spec',STATUS='OLD',ACTION='READ')

  READ(110,*) lm_wr
  DO j=1,5
     DO i=1,nt_wr
        READ(110,*) wrn_logt(i),wr_z(j)
        READ(110,*) tspec_wrn(:,i,j)
     ENDDO
  ENDDO
  CLOSE(110)
  
 
  !interpolate to the main array
  DO j=1,nz
     
     i1 = MIN(MAX(locate(wr_z,z_isoc(j)),1),4) 
     aa = (z_isoc(j)-wr_z(i1))/(wr_z(i1+1)-wr_z(i1))


     IF (aa .LT. 0.0)  aa=0.0
     IF (aa .GT. 1.0)  aa=1.0

     bb = 1-aa

     DO i=1,nt_wr
        i_spec_wr = bb*(tspec_wrn(:,i,i1)+small_value) + &
             aa*(tspec_wrn(:,i,i1+1)+small_value)
        wrn_spec(:,i,j) = linterparr(lm_wr,i_spec_wr,lambda) 
     ENDDO
  ENDDO

  !read in WR-C Teff and spectra
  OPEN(120,FILE='../SPECTRA/WR_spectra/CMFGEN_WC'//&
       '.spec',STATUS='OLD',ACTION='READ')

  READ(120,*) lm_wr
  DO j=1,5
     DO i=1,nt_wr
        READ(120,*) wrc_logt(i),wr_z(j)
        READ(120,*) tspec_wrc(:,i,j)
     ENDDO
  ENDDO
  CLOSE(120)


  !interpolate to the main spectral array
  DO j=1,nz
     i1 = MIN(MAX(locate(wr_z,z_isoc(j)),1),4)
     aa = (z_isoc(j)-wr_z(i1))/(wr_z(i1+1)-wr_z(i1))

     IF (aa .LT. 0.0)  aa=0.0
     IF (aa .GT. 1.0)  aa=1.0

     bb = 1-aa


     DO i=1,nt_wr
        i_spec_wr = bb*(tspec_wrc(:,i,i1)+small_value) + &
             aa*(tspec_wrc(:,i,i1+1)+small_value)
        wrc_spec(:,i,j) = linterparr(lm_wr,i_spec_wr,lambda) 
     ENDDO
  ENDDO

  !------------------------------------



END SUBROUTINE SPEC_WR


