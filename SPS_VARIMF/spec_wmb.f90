SUBROUTINE SPEC_WMB()

  !Read WMBasic grid from BPASS dataset (Eldridge et al. 2017)

  USE params; USE sub_func
  IMPLICIT NONE


  INTEGER :: i,j,k,z
  INTEGER, PARAMETER :: nspec_wmb=1963
  REAL(KIND(1.d0)), DIMENSION(nspec_wmb) :: lm_wmb=0.0
  REAL(KIND(1.d0)) :: aa,bb

  REAL(KIND(1.d0)), DIMENSION(nspec_wmb,nz_wmb,nt_wmb,ng_wmb) :: i_spec_wmb=0.
  CHARACTER(6) :: z_type


 
  !------------------------------

  !read in Teff array
  OPEN(101,FILE='../SPECTRA/WMBASIC_teff.txt',&
       STATUS='OLD',ACTION='READ')

  DO i=1,nt_wmb
     READ(101,*) wmb_logt(i),wmb_t(i)
  ENDDO
  CLOSE(101)


  OPEN(102,FILE='../SPECTRA/WMBASIC_z.txt',STATUS='OLD',ACTION='READ')

  DO z=1,nz_wmb
     READ(102,*) wmb_z(z)
     WRITE(z_type,'(F6.4)') wmb_z(z)

     OPEN(110,FILE='../SPECTRA/WMBASIC_OB/sg_'//&
         z_type//'Z.dat',STATUS='OLD',ACTION='READ')
     OPEN(120,FILE='../SPECTRA/WMBASIC_OB/dw_'//&
        z_type//'Z.dat',STATUS='OLD',ACTION='READ')
     OPEN(130,FILE='../SPECTRA/WMBASIC_OB/hg_'//&
        z_type//'Z.dat',STATUS='OLD',ACTION='READ')


     DO i=1,nspec_wmb
        lm_wmb(i)=0.0

         !spectra are interpolated to the Basel spectra grid
         READ(110,*) lm_wmb(i),(i_spec_wmb(i,z,k,1),k=1,nt_wmb)
         READ(120,*) lm_wmb(i),(i_spec_wmb(i,z,k,2),k=1,nt_wmb)
         READ(130,*) lm_wmb(i),(i_spec_wmb(i,z,k,3),k=1,nt_wmb)

     ENDDO

     CLOSE(110)
     CLOSE(120)
     CLOSE(130)

  ENDDO

  CLOSE(102)


     wmb_spec(:,:,:,:) = i_spec_wmb (:,:,:,:)
    !----------------------------------------------------------------!


END SUBROUTINE SPEC_WMB


