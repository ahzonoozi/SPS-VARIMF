# SPS-VARIMF
The Stellar Population Synthesis for Varying Initial Mass Function (SPS-VARIMF) model is designed to compute the spectral and photometric evolution of stellar populations across a broad parameter space. This model is explicitly developed to investigate the impact of a non-universal Initial Mass Function (IMF) on key galactic properties, including integrated luminosities, stellar masses, and the fractional contributions of stellar remnants. 
#IMF Variation in SPS-VARIMF
In SPS-VARIMF, the galaxy-wide IMF (gwIMF) is allowed to vary at each time step based on the galaxy's star formation rate (SFR) and the gas-phase metallicity at that epoch. The code includes both an invariant IMF and the Integrated Galaxy-wide IMF (IGIMF) framework, but it can also accommodate any other variable IMF prescription through private communication. This flexibility makes SPS-VARIMF a powerful tool for exploring the consequences of IMF variability in different astrophysical contexts.
Metallicity Evolution
The model simulates metal enrichment processes under the assumption of a closed-box system, meaning no gas inflow or outflow is considered. The metallicity evolution is inherently dependent on the adopted star formation efficiency (SFE), which governs the rate of metal production and retention. However, users have the flexibility to assume a constant metallicity if desired, allowing for controlled comparisons and alternative modeling approaches.
#Dust Attenuation
An optional dust attenuation component is included in the model, implemented using a power-law formalism to account for the wavelength-dependent absorption and scattering of light by interstellar dust. This feature allows users to incorporate the effects of dust on spectral energy distributions (SEDs) and photometric properties, providing a more realistic comparison with observed galaxy spectra.


For any publications or research that utilize SPS-VARIMF , please cite
the following reference:

Zonoozi, Haghi, and Kroupa et al. MNRAS, 2025 (https://doi.org/10.1093/mnras/staf091)
“Stellar population synthesis models with a physically varying IMF”
