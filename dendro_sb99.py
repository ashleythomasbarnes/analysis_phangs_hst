import os
from astropy.table import Table
import astropy.units as au
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from tqdm.auto import tqdm 

from analysis_phangs_hst import dendro_misc

def get_SB99models(inputdir='/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/sb99/fiducial/fiducial_6p0linear/', tmin = 0, tmax = 1e8):
    
    """Get SB99 models"""

    t_ewidth = Table.read('%s/fiducial_6p0linear.ewidth1' %(inputdir), format='ascii', header_start=3, data_start=4)
    names = ['TIME', 'QHI', 'QHIf', 'QHeI', 'QHeIf', 'QHeII', 'QHeIIf', 'logL']
    t_quanta = Table.read('%s/fiducial_6p0linear.quanta1' %(inputdir), format='ascii', data_start=5, names=names)
    names = ['TIME', 'ALLp', 'OBp', 'RSGp', 'LBVp', 'WRp', 'ALLe', 'ALLm', 'OBm', 'RSGm', 'LBVm', 'WRm']
    t_power = Table.read('%s/fiducial_6p0linear.power1' %(inputdir), format='ascii', data_start=5, names=names)
    names = ['TIME', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Mdotwind', 'Mdotsne','Mdotboth', 'Mtot']
    t_yield = Table.read('%s/fiducial_6p0linear.yield1' %(inputdir), format='ascii', data_start=5, names=names)

    time = t_ewidth['TIME']
    mask = ((time>=tmin) & (time<=tmax))
    time = time[mask]*au.yr

    Q = t_quanta['QHI'][mask]
    Mdot = t_yield['Mdotwind'][mask]
    Lmech = t_power['ALLp'][mask]
    Lbol = t_quanta['logL'][mask]
    LHa = t_ewidth['LUM(H_A)'][mask]
    Lfrac = 10**(t_quanta['logL'][mask]-t_ewidth['LUM(H_A)'][mask])
    
    return({'time':time, 'mdot':Mdot, 'lmech':Lmech, 'lbol':Lbol, 'lha':LHa, 'lfrac':Lfrac})


def get_SB99models_arr(inputdir_='/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/sb99/fiducial/'):

    masses = ['4p0', '4p5', '5p0', '5p5', '6p0', '6p5']
    masses_ = [1e4, 5e4, 1e5, 5e5, 1e6, 5e6]
    n = len(masses)

    t_ewidth = dict.fromkeys(masses)
    t_quanta = dict.fromkeys(masses)
    t_power = dict.fromkeys(masses)
    t_yield = dict.fromkeys(masses)

    for mass, mass_,i in zip(masses, masses_, range(n)):


        inputdir = '%s/fiducial_%s' %(inputdir_, mass) 

        t_ewidth_ = Table.read('%s/fiducial_%s.ewidth1' %(inputdir, mass), format='ascii', header_start=3, data_start=4)

        names = ['TIME', 'QHI', 'QHIf', 'QHeI', 'QHeIf', 'QHeII', 'QHeIIf', 'logL']
        t_quanta_ = Table.read('%s/fiducial_%s.quanta1' %(inputdir, mass), format='ascii', data_start=5, names=names)

        names = ['TIME', 'ALLp', 'OBp', 'RSGp', 'LBVp', 'WRp', 'ALLe', 'ALLm', 'OBm', 'RSGm', 'LBVm', 'WRm']
        t_power_ = Table.read('%s/fiducial_%s.power1' %(inputdir, mass), format='ascii', data_start=5, names=names)

        names = ['TIME', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Mdotwind', 'Mdotsne','Mdotboth', 'Mtot']
        t_yield_ = Table.read('%s/fiducial_%s.yield1' %(inputdir, mass), format='ascii', data_start=5, names=names)

        t_ewidth_['TIME'] = np.log10(t_ewidth_['TIME'])
        t_quanta_['TIME'] = np.log10(t_quanta_['TIME'])
        t_power_['TIME'] = np.log10(t_power_['TIME'])
        t_yield_['TIME'] = np.log10(t_yield_['TIME'])

        t_ewidth_['MASS'] = np.array([mass_]*len(t_power_['TIME']))
        t_quanta_['MASS'] = np.array([mass_]*len(t_power_['TIME']))
        t_power_['MASS'] = np.array([mass_]*len(t_power_['TIME']))
        t_yield_['MASS'] = np.array([mass_]*len(t_power_['TIME']))

        if i == 0: 
            mass_arr = t_ewidth_['MASS']
            time_arr = t_ewidth_['TIME']
            mdot_arr = t_yield_['Mdotwind']
            lmec_arr = t_power_['ALLp']
        else: 
            mass_arr = np.vstack([mass_arr,t_ewidth_['MASS']])
            time_arr = np.vstack([time_arr,t_ewidth_['TIME']])
            mdot_arr = np.vstack([mdot_arr,t_yield_['Mdotwind']]) 
            lmec_arr = np.vstack([lmec_arr,t_power_['ALLp']]) 

    mass_arr = np.log10(mass_arr)
    
    return({'mass':mass_arr, 'time':time_arr, 'mdot':mdot_arr, 'lmec':lmec_arr})

def get_sb99props(SB99models_arr, props_all, showplots=True):

    def func(x, a, b):
        return (a*x)+b

    age = props_all['reg_dolflux_Age_MinChiSq'].quantity.to('yr')
    mass = props_all['reg_dolflux_Mass_MinChiSq'].quantity.to('Msun')

    xdata = np.unique(SB99models_arr['mass'])
    ydata = np.log10(np.nanmean((10**SB99models_arr['mdot'])/(10**SB99models_arr['mdot'][0]),axis=1))

    popt_mdot, _ = curve_fit(func, xdata, ydata)

    if showplots: 
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(1,2,1)
        ax.scatter(xdata,ydata)
        ax.plot(xdata, func(xdata, *popt_mdot), 'g--', label='fit: a=%5.3f, b=%5.3f' % tuple(popt_mdot))
        ax.grid(alpha=0.3, linestyle=':')
        ax.legend(fontsize=9)
        ax.set_xlabel('Mass log(Msun)')
        ax.set_ylabel('Mass loss rate (normalised to Mclust=1e4sun)')
        plt.show()

    xdata = np.unique(SB99models_arr['mass'])
    ydata = np.log10(np.nanmean((10**SB99models_arr['lmec'])/(10**SB99models_arr['lmec'][0]),axis=1)) #fraction in logspace

    popt_lmec, _ = curve_fit(func, xdata, ydata)

    if showplots: 
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(1,2,1)
        ax.scatter(xdata,ydata)
        ax.plot(xdata, func(xdata, *popt_lmec), 'g--', label='fit: a=%5.3f, b=%5.3f' % tuple(popt_lmec))
        ax.grid(alpha=0.3, linestyle=':')
        ax.legend(fontsize=9)
        ax.set_xlabel('Mass log(Msun)')
        ax.set_ylabel('Mech. lum (normalised to Mclust=1e4sun)')
        plt.show()


    lmec = np.empty(len(props_all['_idx']))*np.nan
    mdot = np.empty(len(props_all['_idx']))*np.nan

    n = len(props_all)
    for regionID in tqdm(range(n)):

        #Get in correct units
        clust_mass_ = np.log10(mass.value[regionID]) 
        clust_age_ = np.log10(age.to('yr').value[regionID])

        # Values for a 10^4Msun cluster, of the age of the obs cluster 
        _, id_ = dendro_misc.find_nearest(SB99models_arr['time'][0,:], clust_age_)
        lmec_ = SB99models_arr['lmec'][0,id_]
        mdot_ = SB99models_arr['mdot'][0,id_]

        # Values of offset from a 10^4Msun cluster, due to more or less mass of obs cluster 
        cf_lmec = func(clust_mass_, *popt_lmec)
        cf_mdot = func(clust_mass_, *popt_mdot)

        # Values of 10^4Msun cluster +- offset due to obs mass
        lmec[regionID] = lmec_+cf_lmec
        mdot[regionID] = mdot_+cf_mdot

    lmec = (10**lmec) *au.erg/au.s
    mdot = ((10**mdot) *au.Msun/au.yr).to('g/s')

    windvelo = np.sqrt((2.0*lmec)/mdot)
    windvelo = windvelo.to('km/s')
    
    return(lmec, mdot, windvelo)