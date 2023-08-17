import numpy as np
import astropy.constants as ac

def get_pdir(lbol, radius):
    """Convert bolometric luminoisty to direct radiation pressure
        Taken from Lopez et al. (2014)"""
    pdir = (3 * lbol.to('Lsun')) / (4*np.pi * (radius.to('cm')**2) * ac.c.to('cm/s') ) / ac.k_B
    return pdir.to('K/cm^3')

def get_pth(ne, te, ionisation=2):
    """ionisation=2 assumes single ionised He"""
    pth = ne*te*ionisation
    return(pth.to('K/cm^3'))

def get_pwind(mdot, windvelo, radius):
    pwind = (mdot * windvelo) / ((4./3.)*np.pi * (radius.to('cm')**2))
    pwind = (pwind/ac.k_B.cgs)
    return(pwind.to('K/cm^3'))

def save_fits(props_all_final, props_all_file):
    filename = props_all_file.replace('02','03')
    print('[INFO] [save_fits] Saved to %s' %filename)
    props_all_final.write(filename, overwrite=True)