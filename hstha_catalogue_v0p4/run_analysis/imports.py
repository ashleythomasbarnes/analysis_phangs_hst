import numpy as np
from astropy import units as u 
from astropy import constants as const
import matplotlib.pyplot as plt
from astropy.table import QTable, join
from scipy.optimize import curve_fit
import pyneb as pn
from matplotlib import colors
from scipy import stats
from astropy.visualization import wcsaxes

# import numpy as np
# from astropy import units as u 
# from astropy.io import fits
# import matplotlib.pyplot as plt
# import aplpy
# import colorcet
# import matplotlib as mpl
# from astropy.table import QTable
# from scipy import stats

import numpy as np
from astropy import units as u 
from astropy.io import fits
import matplotlib.pyplot as plt
import aplpy
import colorcet
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import imageio.v3 as iio

import warnings
warnings.filterwarnings('ignore')   

plt.style.use('paper.mplstyle')
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

import warnings
warnings.filterwarnings('ignore')

import sys 
sys.path.append('../../../../misc/data_cube_analysis')  
from bindata import * 
from histograms import *

sys.path.append('../../../../misc/data_cube_analysis') 
from Sun_Plot_Tools.sun_plot_tools.ax import *

sys.path.append('./../modules')
import cat_props, cat_misc

from tools import *