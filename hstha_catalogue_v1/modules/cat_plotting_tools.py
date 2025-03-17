import numpy as np
from astropy import units as u 
import matplotlib.pyplot as plt
from astropy.table import QTable, join
from scipy.optimize import curve_fit
import numpy as np
from astropy import units as u 
import matplotlib.pyplot as plt
import matplotlib as mpl

import warnings
warnings.filterwarnings('ignore')  

# plt.style.use('paper.mplstyle')
# mpl.rcParams['xtick.minor.visible'] = True
# mpl.rcParams['ytick.minor.visible'] = True

def find_nearest(array, value):
    """find nearest value in array, and return id"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def get_hist(data, bins=None, nbins=50, logbins=False, norm=True, cum=False):

    data = data.flatten()

    if bins.all() is None:
        vmin=np.nanmin(data)
        vmax=np.nanmax(data)

        bmin = vmin - (np.absolute(vmin)*1)
        bmax = vmax + (np.absolute(vmax)*0.3)

        if logbins:
            min = np.nanmin(data[data>0])
            bins = np.logspace(np.log10(bmin), np.log10(bmax), nbins+1)
        else:
            bins = np.linspace(bmin, bmax, nbins+1)
    else:
        nbins = len(bins)-1

    bins_cent = np.empty([nbins])

    for i in range(nbins):
        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])

    hist = np.histogram(data.flatten(), bins)[0]

    if cum:
        hist = np.cumsum(hist)
    if norm:
        hist = hist/np.nanmax(hist)

    return(bins, bins_cent, hist)