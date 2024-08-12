import numpy as np
from astropy import units as u 
from astropy import constants as const
import matplotlib.pyplot as plt
from astropy.table import QTable, join
from scipy.optimize import curve_fit
import pyneb as pn
from matplotlib import colors

plt.style.use('paper.mplstyle')

import warnings
warnings.filterwarnings('ignore')

# sys.path.append('../run_catalogue/')
# from imports import *  

import sys 
sys.path.append('../../../../misc/data_cube_analysis')  
from bindata import * 
from histograms import *

sys.path.append('../../../../misc/data_cube_analysis') 
from Sun_Plot_Tools.sun_plot_tools.ax import *