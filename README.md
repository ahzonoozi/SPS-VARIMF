# SPS-VARIMF
The Stellar Population Synthesis for Varying Initial Mass Function (SPS-VARIMF) model is designed to compute the spectral and photometric evolution of stellar populations across a broad parameter space. This model is explicitly developed to investigate the impact of a non-universal Initial Mass Function (IMF) on key galactic properties, including integrated luminosities, stellar masses, and the fractional contributions of stellar remnants. 

## IMF Variation in SPS-VARIMF

In SPS-VARIMF, the galaxy-wide IMF (gwIMF) is allowed to vary at each time step based on the galaxy's star formation rate (SFR) and the gas-phase metallicity at that epoch. The code includes both an invariant IMF and the Integrated Galaxy-wide IMF (IGIMF) framework, but it can also accommodate any other variable IMF prescription through private communication. This flexibility makes SPS-VARIMF a powerful tool for exploring the consequences of IMF variability in different astrophysical contexts.

## Metallicity Evolution

The model simulates metal enrichment processes under the assumption of a closed-box system, meaning no gas inflow or outflow is considered. The metallicity evolution is inherently dependent on the adopted star formation efficiency (SFE), which governs the rate of metal production and retention. However, users have the flexibility to assume a constant metallicity if desired, allowing for controlled comparisons and alternative modeling approaches.

## Dust Attenuation

An optional dust attenuation component is included in the model, implemented using a power-law formalism to account for the wavelength-dependent absorption and scattering of light by interstellar dust. This feature allows users to incorporate the effects of dust on spectral energy distributions (SEDs) and photometric properties, providing a more realistic comparison with observed galaxy spectra.


 ## Setting Up the Model:
   
Before compiling, ensure that the SPS-VARIMF directory is located in the same parent directory as the following required folders:


    • OUTPUTS
    
    • SPECTRA
    
    • ISOCHRONS
    
    • FILTERS
    
    • YIELD

 Extract both BaSeL_1.zip and BaSeL_2.zip into the following directory: 
 
    SPECTRA/BaSeL3.1/

Make sure all extracted files are placed directly inside the BaSeL3.1 folder.

## Compiling the Code
The code includes a Makefile to facilitate compilation. Before proceeding, ensure that your system has a Fortran compiler installed (e.g., gfortran). You may need to modify the Makefile to specify the correct Fortran compiler for your system

    1. Clean any previous builds:
         make clean
    2. Compile the program:
         make
        This will generate an executable file named main.exe in the SPS-VARIMF directory.
    3. Running the executable:
          ./main.exe    




## Adjust parameters :
To properly configure the SPS-VARIMF model, key simulation parameters must first be defined and initialized  in  params.f90.

In order to adjust optional cases 
       
   #### Dust attenuation (dust_atten):
         
                       = 0 : Disabled   
                       = 1 : Enabled (power-law attenuation)

                       
   #### Remnant mass calculation (remnant_cal):
   
                       = 0 :  Renzini et al. (1993) prescription      
                       = 1 :  Spera et al. (2015) prescription

                       
   #### Metallicity evolution mode (Z_MODE):
   
                       = 0 : Constant metallicity 
                       = 1 : Metallicity computed at each step from gas enrichment

                       
   #### OB stars count calculation (OB_stars):
   
                       = 0 : Do not compute the number of OB stars    
                       = 1 : Compute the number of OB stars

                       
                       
## Model Configuration Parameters

Modify the simulation parameters  in the main.f90 program to setup the model. The key parameters include:

  #### Initial Mass Function (IMF) (imf_type):
  
  Specifies the type of IMF used in the simulation:
  
    = 2  : Canonical IMF (Kroupa 2001) – A standard IMF with fixed power-law slopes.
    = 3  : Integrated Galactic IMF (IGIMF) – An SFR- and metallicity-dependent IMF.  
    = 5  :Top-Heavy IMF – A flatter IMF favoring high-mass star formation.

 #### Isochrones (isoc_type):
 
 Determines the stellar evolution tracks used:
 
    = 1 : Padova (2007)
    = 2 : PARSEC (2022)

 #### Spectral Library (spec_type):
 
 Specifies the spectral database used for stellar populations:
 
    = 1 : BaSeL spectral library (default)
    = 2 : Miles spectral library (in progress, not yet available)

 #### Star Formation History (SFH) (sfh_type):
 
 Defines how the galaxy forms stars over time:
 
    = 0 : Single Stellar Population (SSP) – All stars form in a single burst.
    = 1 : Constant Star Formation Rate (SFR) – Stars form at a steady rate over time.
    = 2 : Exponentially Declining SFR – Star formation rate decreases exponentially. 
    = 4 : Delayed-Tau SFH – A rising SFR phase followed by exponential decline.
    
 #### Additional SFH Parameters:
 
    • e-folding timescale (tau) (for sfh_type = 2 and 4):
        ◦ Range: 0.001 – 1000.0 Gyr
        
    • Star formation start and truncation time (Tstart, Ttrunc):
        ◦ Default: Tstart = 0.0, Ttrunc = 0.0
        ◦ Defined in Gyr, based on the age of the universe.

 #### Metallicity (zmet_ini):
 
Specifies the initial metallicity index:

    • Choose a value from 1–22, corresponding to metallicity values in Table1.
    • (Z_MODE=0 )  Constant metallicity throughout the simulation (zmet = zmet_ini).
    • (Z_MODE=1 )  Self-enrichment: Metallicity evolves over time due to stellar feedback.
          

 #### Galaxy Mass (M_galaxy or M_UCD):
 
 Defines the total stellar mass formed over time:
 
    • M_galaxy → Total stellar mass of a galaxy (in solar masses).
    
    • M_UCD → Monolithic formation of ultra-compact dwarf galaxies (UCDs), following a top-heavy IMF.

 #### Star Formation Efficiency (f_star):
 
 Determines how efficiently gas converts into stars:
 
  f_star=M_galaxy/M_gas​​
  
    • Valid range: 0 < f_star < 1 


## OUTPUT Folder
The OUTPUT/ directory stores the simulation results. It contains two main types of output files:

1️.mags Files:

These files store essential evolutionary parameters of the simulated stellar population, including time, stellar mass, luminosity, star formation rate, mass of remnants, mass of white dwarfs, metallicity, the ratio of gas-phase mass to total galaxy mass, and magnitudes in 143 photometric bands as listed in FILTER_LIST.dat, all recorded at each evolutionary time step.
These .mags files are crucial for analyzing the photometric evolution of the galaxy.

2️.spec Files:

These files contain the wavelength array and the corresponding spectra (fν) in L☉/Hz units at each evolutionary time step.



 ## Acknowledgments
This program is built upon extensive research in stellar population synthesis and uses publicly available codes and libraries, including Padova, PARSEC, BaSeL, GALAXEV (Bruzual & Charlot), FSPS (Flexible Stellar Population Synthesis), and others. We acknowledge the invaluable contributions of the researchers and developers who have made these resources available to the scientific community.
SPS-VARIMF was developed by Akram Hasani Zonoozi at the Institute for Advanced Studies in Basic Sciences (IASBS).  For any publications or research that utilize SPS-VARIMF , please cite the following reference:

Zonoozi, Haghi, and Kroupa et al. MNRAS, 2025 (https://doi.org/10.1093/mnras/staf091)

“Stellar population synthesis models with a physically varying IMF”

For comments, questions, bug reports, or feedback, please contact:

 a.hasani@iasbs.ac.ir
 e.h.zonoozi@gmail.com
 
This program is free software and is distributed under an open-source license. You are free to redistribute and modify it as needed.
