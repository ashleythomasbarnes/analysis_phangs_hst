import os
from astropy.table import Table
import astropy.units as au
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from tqdm.auto import tqdm 


from analysis_phangs_hst import dendro_misc

def get_SB99models(inputdir='/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/sb99/fiducial/fiducial_6p0linear/', tmin = 0, tmax = 1e8):
    """
    Extracts SB99 stellar population synthesis model outputs from specified directory.
    
    Parameters:
    - inputdir (str): Directory containing the SB99 model output files.
    - tmin (float): Minimum time boundary to filter the data.
    - tmax (float): Maximum time boundary to filter the data.
    
    Returns:
    - dict: A dictionary containing time, mass loss rate, mechanical luminosity, 
            bolometric luminosity, H-alpha luminosity, and luminosity fraction.
    """
    
    # Read equivalent width table
    t_ewidth = Table.read(os.path.join(inputdir, 'fiducial_6p0linear.ewidth1'), format='ascii', header_start=3, data_start=4)
    
    # Define column names and read quanta table
    names = ['TIME', 'QHI', 'QHIf', 'QHeI', 'QHeIf', 'QHeII', 'QHeIIf', 'logL']
    t_quanta = Table.read(os.path.join(inputdir, 'fiducial_6p0linear.quanta1'), format='ascii', data_start=5, names=names)
    
    # Define column names and read power table
    names = ['TIME', 'ALLp', 'OBp', 'RSGp', 'LBVp', 'WRp', 'ALLe', 'ALLm', 'OBm', 'RSGm', 'LBVm', 'WRm']
    t_power = Table.read(os.path.join(inputdir, 'fiducial_6p0linear.power1'), format='ascii', data_start=5, names=names)
    
    # Define column names and read yield table
    names = ['TIME', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Mdotwind', 'Mdotsne','Mdotboth', 'Mtot']
    t_yield = Table.read(os.path.join(inputdir, 'fiducial_6p0linear.yield1'), format='ascii', data_start=5, names=names)

    # Filter the tables based on the specified time boundaries
    time = t_ewidth['TIME']
    mask = ((time >= tmin) & (time <= tmax))
    time = time[mask] * au.yr

    # Extract specific columns from tables after masking
    Q = t_quanta['QHI'][mask]
    Mdot = t_yield['Mdotwind'][mask]
    Lmech = t_power['ALLp'][mask]
    Lbol = t_quanta['logL'][mask]
    LHa = t_ewidth['LUM(H_A)'][mask]
    Lfrac = 10 ** (t_quanta['logL'][mask] - t_ewidth['LUM(H_A)'][mask])
    
    return {'time': time, 'mdot': Mdot, 'lmech': Lmech, 'lbol': Lbol, 'lha': LHa, 'lfrac': Lfrac}

def get_SB99models_arr(inputdir_='/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/sb99/fiducial/'):
    """
    Extracts SB99 models for different stellar masses and organizes them into arrays.
    
    Parameters:
    - inputdir_ (str): Base directory containing SB99 model directories for different masses.
    
    Returns:
    - dict: A dictionary containing mass, time, mass loss rate, and mechanical luminosity arrays.
    """
    
    # List of model masses in strings and their corresponding values
    masses = ['4p0', '4p5', '5p0', '5p5', '6p0', '6p5']
    masses_ = [1e4, 5e4, 1e5, 5e5, 1e6, 5e6]
    n = len(masses)

    # Initialize dictionaries to store tables for each mass
    t_ewidth = dict.fromkeys(masses)
    t_quanta = dict.fromkeys(masses)
    t_power = dict.fromkeys(masses)
    t_yield = dict.fromkeys(masses)

    # Loop over each mass and read corresponding tables
    for mass, mass_, i in zip(masses, masses_, range(n)):
        inputdir = os.path.join(inputdir_, f'fiducial_{mass}')

        # Read equivalent width, quanta, power, and yield tables for the current mass
        t_ewidth_ = Table.read(os.path.join(inputdir, f'fiducial_{mass}.ewidth1'), format='ascii', header_start=3, data_start=4)
        names = ['TIME', 'QHI', 'QHIf', 'QHeI', 'QHeIf', 'QHeII', 'QHeIIf', 'logL']
        t_quanta_ = Table.read(os.path.join(inputdir, f'fiducial_{mass}.quanta1'), format='ascii', data_start=5, names=names)
        names = ['TIME', 'ALLp', 'OBp', 'RSGp', 'LBVp', 'WRp', 'ALLe', 'ALLm', 'OBm', 'RSGm', 'LBVm', 'WRm']
        t_power_ = Table.read(os.path.join(inputdir, f'fiducial_{mass}.power1'), format='ascii', data_start=5, names=names)
        names = ['TIME', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Mdotwind', 'Mdotsne','Mdotboth', 'Mtot']
        t_yield_ = Table.read(os.path.join(inputdir, f'fiducial_{mass}.yield1'), format='ascii', data_start=5, names=names)

        # Convert time values to logarithmic scale
        for table in [t_ewidth_, t_quanta_, t_power_, t_yield_]:
            table['TIME'] = np.log10(table['TIME'])
            table['MASS'] = np.array([mass_] * len(table['TIME']))

        # Organize data into arrays
        if i == 0: 
            mass_arr = t_ewidth_['MASS']
            time_arr = t_ewidth_['TIME']
            mdot_arr = t_yield_['Mdotwind']
            lmec_arr = t_power_['ALLp']
        else: 
            mass_arr = np.vstack([mass_arr, t_ewidth_['MASS']])
            time_arr = np.vstack([time_arr, t_ewidth_['TIME']])
            mdot_arr = np.vstack([mdot_arr, t_yield_['Mdotwind']]) 
            lmec_arr = np.vstack([lmec_arr, t_power_['ALLp']])

    # Convert mass values to logarithmic scale
    mass_arr = np.log10(mass_arr)
    
    return {'mass': mass_arr, 'time': time_arr, 'mdot': mdot_arr, 'lmec': lmec_arr}


def get_sb99props(SB99models_arr, props_all, showplots=True):
    """
    Calculate mechanical luminosity, mass loss rate, and wind velocity based on SB99 models.
    
    Parameters:
    - SB99models_arr (dict): A dictionary containing SB99 model data arrays for different masses.
    - props_all (dict): A dictionary with properties such as age and mass.
    - showplots (bool): A flag to indicate if plots should be displayed.
    
    Returns:
    - tuple: Mechanical luminosity, mass loss rate, and wind velocity.
    """
    
    # Define a linear function to fit the data
    def func(x, a, b):
        return (a*x)+b

    # Convert age and mass to the appropriate units
    age = props_all['reg_dolflux_Age_MinChiSq'].quantity.to('yr')
    mass = props_all['reg_dolflux_Mass_MinChiSq'].quantity.to('Msun')

    # Extract unique mass values and calculate mean normalized mass loss rate
    xdata = np.unique(SB99models_arr['mass'])
    ydata = np.log10(np.nanmean((10**SB99models_arr['mdot'])/(10**SB99models_arr['mdot'][0]),axis=1))

    # Fit the linear function to the data
    popt_mdot, _ = curve_fit(func, xdata, ydata)

    # Plot the data and the fit if the flag is True
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

    # Calculate mean normalized mechanical luminosity
    ydata = np.log10(np.nanmean((10**SB99models_arr['lmec'])/(10**SB99models_arr['lmec'][0]),axis=1))

    # Fit the linear function to the data
    popt_lmec, _ = curve_fit(func, xdata, ydata)

    # Plot the data and the fit if the flag is True
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

    # Initialize arrays for results
    lmec = np.empty(len(props_all['_idx']))*np.nan
    mdot = np.empty(len(props_all['_idx']))*np.nan

    n = len(props_all)
    for regionID in tqdm(range(n), desc='Calculate Lmech, Mdot, Vwind'):
        # Convert to log scale and appropriate units
        clust_mass_ = np.log10(mass.value[regionID]) 
        clust_age_ = np.log10(age.to('yr').value[regionID])

        # Get the closest model values for a 10^4Msun cluster of the observed age
        _, id_ = dendro_misc.find_nearest(SB99models_arr['time'][0,:], clust_age_)
        lmec_ = SB99models_arr['lmec'][0,id_]
        mdot_ = SB99models_arr['mdot'][0,id_]

        # Calculate the offsets for observed cluster mass compared to 10^4Msun cluster
        cf_lmec = func(clust_mass_, *popt_lmec)
        cf_mdot = func(clust_mass_, *popt_mdot)

        # Adjust the 10^4Msun cluster values using the calculated offsets
        lmec[regionID] = lmec_+cf_lmec
        mdot[regionID] = mdot_+cf_mdot

    # Convert values to appropriate units
    lmec = (10**lmec) *au.erg/au.s
    mdot = ((10**mdot) *au.Msun/au.yr).to('g/s')

    # Calculate wind velocity using the relation v^2 = 2*L/m_dot
    windvelo = np.sqrt((2.0*lmec)/mdot)
    windvelo = windvelo.to('km/s')
    
    return(lmec, mdot, windvelo)
