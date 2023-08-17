import numpy as np
import astropy.constants as ac
import astropy.units as au

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