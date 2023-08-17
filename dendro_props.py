import numpy as np
import astropy.constants as ac
import astropy.units as au
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

from analysis_phangs_hst import dendro_misc

def correct_ha_flux(flux_ha, props_all):
    
    # Calculates the bolometric luminosity from H-alpha flux
    correction_factor = props_all['HA6562_FLUX_CORR']/props_all['HA6562_FLUX']
    
    # Correct H-alpha flux for extinction
    corrected_flux_ha = flux_ha * correction_factor
    
    return(corrected_flux_ha)

def calculate_luminosity(flux_ha, distance):

    # Convert H-alpha flux to luminosity
    luminosity_ha = 4 * np.pi * (distance.quantity ** 2) * flux_ha
    
    return luminosity_ha.to('erg/s')

def get_lbol(lha, conv=17.684):
    """Convert Lhalpha luminoisty to Lbolometric luminosity
        conv = 138: Taken from Kennicutt & Evans (2012) - Lopez et al. (2014)
        conv = 17.684"""
    lbol=conv*lha
    return lbol


def func_reccoeff(temp):
    """Get recomination recoefficent"""
    alpha = 2.753e-14 * (315614 /  temp)**1.500 / ((1.0 + (115188 / temp)**0.407)**2.242)
    return(np.array(alpha) *au.cm**3/au.s)


def get_ne(ha_lum, radius, temp):

    volume = (4./3.)*np.pi*radius.quantity**3 #in units of cm^3
    ha_lum = ha_lum.to('erg/s')
    volume = volume.to('cm^3')
    
    #Recombination rate
    ha_photon_energy = 3.02e-12   # Energy in erg of single H-alpha photon
    ha_rate = ha_lum / ha_photon_energy   # Emission rate (in photons per s) of H-alpha
    rec_rate = ha_rate / 0.45     # Per Calzetti (2012), ~45% of recombinations result in emission of an H-alpha photon

    # rec_coeff = 2.6e-13    # Recombination rate coefficient: assumes case B, T_e = 1e4 K
    rec_coeff = func_reccoeff(temp.value)

    ne = np.sqrt(rec_rate/(rec_coeff*volume))
    ne = np.array(ne)/au.cm**3
    
    return(ne)


def get_lfrac(props_all, SB99models, showplot=True):
    """
    Calculate the luminosity fraction (lfrac) for each region based on SB99 models and optionally plot it.
    
    Parameters:
    - props_all (dict): A dictionary containing the properties of regions.
    - showplot (bool): Flag indicating if the results should be plotted.
    
    Returns:
    - ndarray: Array of luminosity fractions (lfrac) for each region.
    """
    
    # Initialize the lfrac array based on the number of regions
    lfrac = np.empty(len(props_all))*np.nan

    # Convert age and mass to the appropriate units
    age = props_all['reg_dolflux_Age_MinChiSq'].quantity.to('yr')
    mass = props_all['reg_dolflux_Mass_MinChiSq'].quantity.to('Msun')

    # Loop through each region to compute the luminosity fraction
    for regionID in tqdm(range(len(props_all)), desc='Calculate Lfrac'):

        # Find the nearest time in SB99models matching the age of the region
        _, mask = dendro_misc.find_nearest(SB99models['time'].value, age[regionID].value)

        # Assign the corresponding luminosity fraction to the region
        lfrac[regionID] = SB99models['lfrac'][mask]

    # If the showplot flag is True, visualize the results
    if showplot:
        
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot()

        # Plot the luminosity fraction from SB99models against time
        ax.plot(SB99models['time'].to('Myr'), np.log10(SB99models['lfrac']))

        # Scatter plot of the luminosity fractions for each region, with some jitter added for clarity
        ax.scatter(age.to('Myr')+np.random.normal(size=len(age))*au.Myr*0.1, 
                   np.log10(lfrac)+np.random.normal(size=len(age))*0.05, 
                   c=lfrac, zorder=10, cmap='jet', s=10)

        # Set plot limits and labels
        ax.set_xlim([0,20])
        ax.set_xlabel('time (Myr)')
        ax.set_ylabel('log(Lbol/Lha)')
        ax.grid(alpha=0.3, linestyle=':')
    
    # Return the luminosity fraction array with appropriate units
    return(lfrac *au.dimensionless_unscaled)
