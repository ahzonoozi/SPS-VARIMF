MODULE PARAMS
  
  ! ================================================================
  ! Module for global parameters and shared arrays/constants used 
  ! across the Stellar Population Synthesis (SPS) model.
  ! ================================================================ 

  IMPLICIT NONE
  SAVE

  

  ! =======================
  ! General Configuration
  ! ======================= 
  INTEGER :: imf_type
  INTEGER :: isoc_type
  INTEGER :: spec_type
  INTEGER :: sfh_type


  ! ==========================
  ! IMF Parameters (Kroupa 2001)
  ! ==========================
  REAL(KIND(1.d0)) :: alpha1,alpha2



  ! Dust attenuation options:
  ! 0 = Disabled, 1 = Enabled (power-law attenuation)
  INTEGER :: dust_atten=0

  ! Remnant mass calculation:
  ! 0 = Renzini et al. (1993) prescription
  ! 1 = Spera et al. (2015) prescription
  INTEGER :: remnant_cal=1

  ! Metallicity evolution mode:
  ! 0 = Constant metallicity  
  ! 1 = Metallicity computed at each step from gas enrichment  
  INTEGER :: Z_MODE=1

  ! OB stars count calculation:
  ! 0 = Do not compute OB stars
  ! 1 = Compute OB stars 
  INTEGER :: OB_stars=0
  
  !compute indices 
  !0 = Do not compute indices
  !1 = cmpute indices
  INTEGER :: index_cal=1


  ! ============================
  ! Galaxy and Star Properties
  ! ============================

  REAL(KIND(1.d0)) :: M_galaxy     ! Total stellar mass of the galaxy
  REAL(KIND(1.d0)) :: M_UCD        ! Ultra-Compact Dwarf (UCD) mass 
  REAL(KIND(1.d0)) :: f_star       ! Star formation efficiency: M_star / M_gas 
  REAL(KIND(1.d0)) :: zsun_isoc    ! Solar metallicity in isochrones
  REAL(KIND(1.d0)) :: tau          ! E-folding timescale for SFH (in Gyr)
  REAL(KIND(1.d0)) :: Tstart       ! Starting time of star formation (in Gyr)
  REAL(KIND(1.d0)) :: Ttrunc       ! Truncation time for star formation (in Gyr)


!set parameters based on the addopted  isochrone

  INTEGER, PARAMETER :: num_time=94      ! Number of time steps in isochrones
  INTEGER, PARAMETER :: nz=22            ! Number of metallicity bins
  INTEGER, PARAMETER :: num_isoc=94      ! Number of isochrones per metallicity

!set parameters based on the addopted  spectral library

  !spec_type =1  Basel spectra
  REAL(KIND(1.d0)), PARAMETER :: zsun_spec = 0.020
  INTEGER, PARAMETER :: nzinit=6
  INTEGER, PARAMETER :: nspec=1963

  INTEGER, PARAMETER :: n_index= 56
  ! ======================================
  ! Dust Parameters
  ! ======================================

  REAL(KIND(1.d0)), PARAMETER :: v_lambda = 5500.0          ! Reference wavelength (Angstroms)
  REAL(KIND(1.d0)), PARAMETER :: logt_young_stars = 7.0     ! Log(Teff) for young stars
  REAL(KIND(1.d0)), PARAMETER :: dust_power = -0.7          ! Power-law index for dust attenuation
  REAL(KIND(1.d0)), PARAMETER :: dust_depth1 = 1.0          ! Dust attenuation depth (young stars)
  REAL(KIND(1.d0)), PARAMETER :: dust_depth2 = 0.3          ! Dust attenuation depth (old stars)
  
  ! ====================
  ! Metallicity Config
  ! ====================
  INTEGER :: zmet                 ! Current metallicity index
  INTEGER :: zmet_ini             ! Initial metallicity index

  REAL(KIND(1.d0)) :: t_tage = 0.
  REAL(KIND(1.d0)) :: t_trunc = 0.
  REAL(KIND(1.d0)) :: t_tquench = 0.


  ! ================================
  ! Array Dimensions and Constants
  ! ================================
  INTEGER, PARAMETER :: n_filtrs = 143    ! Number of filters
  INTEGER, PARAMETER :: nm = 2000         ! Maximum number of stellar mass bins
  INTEGER, PARAMETER :: num_rem = 2000    ! Bin count for remnants
  
  !Dimensions of BaSeL library
  INTEGER, PARAMETER :: ndim_logt = 68    ! Log(Teff) grid size
  INTEGER, PARAMETER :: ndim_logg = 19    ! Log(g) grid size
  
  !Dimensions of O-rich, C-rich AGB spectra
  INTEGER, PARAMETER :: n_agb_o = 9       
  INTEGER, PARAMETER :: n_agb_c = 5 
  
  !Dimensions of post-AGB spectra      
  INTEGER, PARAMETER :: nt_pagb = 14      
  INTEGER, PARAMETER :: ng_pagb = 3 
        
  !Dimensions of WR spectra
  INTEGER, PARAMETER :: nt_wr = 12
  
  ! WMBasic grid dimensions      
  INTEGER, PARAMETER :: nt_wmb = 11
  INTEGER, PARAMETER :: nz_wmb = 11
  INTEGER, PARAMETER :: ng_wmb = 3
  !Numerical precision threshold  
  REAL(KIND(1.d0)), PARAMETER :: small_value = 10**(-70.d0)   

 
  ! ============================
  ! Common Variables and Arrays
  ! ============================
  ! Filter definitions
  REAL(KIND(1.d0)), DIMENSION(nspec,n_filtrs) :: Slm_filtr

  !magnitude of the Sun in all filters
  REAL(KIND(1.d0)), DIMENSION(n_filtrs) :: magsun        ! Magnitudes of the Sun

  !spectrum of Sun, for absolute mags of Sun
  REAL(KIND(1.d0)), DIMENSION(nspec)  :: sun_spec=0.

  !common wavelength and frequench arrays
  REAL(KIND(1.d0)), DIMENSION(nspec)  :: lambda=0.,spec_nu=0.0

  !arrays for stellar spectral information in HR diagram
  REAL(KIND(1.d0)), DIMENSION(ndim_logt) :: speclib_logt=0.
  REAL(KIND(1.d0)), DIMENSION(ndim_logg) :: speclib_logg=0.
  REAL(KIND(1.0)), DIMENSION(nspec,nz,ndim_logt,ndim_logg) :: spec_flux=0.

  !arrays for the WMBasic grid
  REAL(KIND(1.d0)), DIMENSION(nt_wmb) :: wmb_logt=0., wmb_t=0.
  REAL(KIND(1.d0)), DIMENSION(nz_wmb) :: wmb_z=0.
  REAL(KIND(1.d0)), DIMENSION(ng_wmb) :: wmb_logg=0.
  REAL(KIND(1.0)), DIMENSION(nspec,nz_wmb,nt_wmb,ng_wmb) :: wmb_spec=0.

  !AGB library (Lancon & Mouhcine 2002)
  REAL(KIND(1.d0)), DIMENSION(nspec,n_agb_o) :: agb_spec_o=0.
  REAL(KIND(1.d0)), DIMENSION(n_agb_o,nz) :: agb_logt_o=0.
  REAL(KIND(1.d0)), DIMENSION(nspec,n_agb_c) :: agb_spec_c=0.
  REAL(KIND(1.d0)), DIMENSION(n_agb_c) :: agb_logt_c=0.

  !post-AGB library (Rauch 2003)
  REAL(KIND(1.d0)), DIMENSION(nspec,nt_pagb,ng_pagb,2) :: pagb_spec=0.
  REAL(KIND(1.d0)), DIMENSION(nt_pagb) :: pagb_logt=0.
  REAL(KIND(1.d0)), DIMENSION(3) :: pagb_logg=0

  !WR library (Smith et al. 2002)
  REAL(KIND(1.d0)), DIMENSION(nspec,nt_wr,nz) :: wrn_spec=0.,wrc_spec=0.
  REAL(KIND(1.d0)), DIMENSION(nt_wr) :: wrn_logt=0.,wrc_logt=0.


  !arrays for the isochrone data
  REAL(KIND(1.d0)), DIMENSION(nz,num_time,nm) :: mact_isoc=0.,logl_isoc=0.,&
       logt_isoc=0.,logg_isoc=0.,ffco_isoc=0.,phase_isoc=0.,&
       mini_isoc=0. 

  INTEGER, DIMENSION(nz,num_time) :: nmass_isoc=0
  REAL(KIND(1.d0)), DIMENSION(nz,num_time) :: timestep_isoc=0.
  REAL(KIND(1.d0)), DIMENSION(nz) :: z_isoc=-99.

  REAL(KIND(1.d0)), DIMENSION(4) :: z_yield=-99.
  REAL(KIND(1.d0)), DIMENSION(2901) :: mass_yield
  REAL(KIND(1.d0)), DIMENSION(4,2901) :: Metal_yield
  
   REAL(KIND(1.d0)), DIMENSION(n_index,8) :: index_info=0.

  REAL(KIND(1.d0)), DIMENSION(num_time) :: time_full=0.

 
END MODULE PARAMS
