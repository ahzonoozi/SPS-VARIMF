SUBROUTINE SPEC_PAGB()

  !read in post-AGB spectra from Rauch 2003.



  USE params; USE sub_func
  IMPLICIT NONE



  INTEGER :: stat=1,i,j,k,l, jj
  INTEGER, PARAMETER :: nspec_pagb=9281
  INTEGER :: n_isoc,z,zmin,zmax

  REAL(KIND(1.d0)) :: logage
  REAL(KIND(1.d0)), DIMENSION(nspec_pagb) :: pagb_lm=0.0
  REAL(KIND(1.d0)), DIMENSION(nspec_pagb,nt_pagb,3,2) :: i_spec_pagb=0.
  !INTEGER,DIMENSION(3) :: pagb_logg

  CHARACTER(7) :: pagb_t_type
  CHARACTER(2) :: pagb_g_type
  !Stefan-Boltzman constant [cgs unit]
  REAL(KIND(1.d0)), PARAMETER :: sigma = 5.670374E-5


  !----------------------------------------------------------------!

 !read in post-AGB NLTE model spectra from Rauch 2003, H_Ni composition
 !The library has 2 metallicities, solar and halo+0.1solar

  OPEN(101,FILE='../SPECTRA/post_AGB_spectra/pagb_teff.txt',&
       STATUS='OLD',ACTION='READ')

  DO j=1,nt_pagb
     READ(101,*) pagb_logt(j)

  ENDDO
  CLOSE(101)
 ! pagb_logt = LOG10(pagb_logt)



  !read in solar and halo metallicity post-AGB spectra
  DO j=1,nt_pagb

    WRITE(pagb_t_type,'(I7.7)') int(pagb_logt(j))


    DO k=1,ng_pagb
       pagb_logg(k) = (k+5)*1.0


       WRITE(pagb_g_type,'(I2.2)')  (k+5)*10
       !WRITE(*,*)pagb_g_type

       OPEN(110,FILE='../SPECTRA/post_AGB_spectra/halo_Z/'//pagb_t_type//'_'//pagb_g_type//'_halo',&
            STATUS='OLD',ACTION='READ')
       OPEN(120,FILE='../SPECTRA/post_AGB_spectra/solar_Z/'//pagb_t_type//'_'//pagb_g_type//'_solar',&
            STATUS='OLD',ACTION='READ')



        Do jj=1,35
           READ(110,*)
           READ(120,*)
        ENDDO



        DO i=1,nspec_pagb

           pagb_lm(i)=0.0
           !halo metallicity post-AGB spectra
           READ(110,*) pagb_lm(i),i_spec_pagb(i,j,k,1)
           i_spec_pagb(i,j,k,1) = (i_spec_pagb(i,j,k,1) * pagb_lm(i)*pagb_lm(i) /&
                    (2.* 3.*pagb_logt(j)**4.*sigma))*10.**(-26)


           !solar metallicity post-AGB spectra
           READ(120,*) pagb_lm(i),i_spec_pagb(i,j,k,2)
           i_spec_pagb(i,j,k,2) = (i_spec_pagb(i,j,k,2) * pagb_lm(i)*pagb_lm(i) /&
                   (2.*3.*pagb_logt(j)**4.*sigma))*10.**(-26)

         ENDDO


         CLOSE(110)
         CLOSE(120)
    ENDDO

  ENDDO


  DO l=1,2
     DO k=1,3

        DO j=1,nt_pagb

           pagb_spec(:,j,k,l) = MAX(linterparr(pagb_lm,i_spec_pagb(:,j,k,l),&
                lambda),small_value)


        ENDDO
     ENDDO
  ENDDO



  pagb_logt = LOG10(pagb_logt)


!================================================
 
END SUBROUTINE SPEC_PAGB


