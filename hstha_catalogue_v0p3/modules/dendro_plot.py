import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table, hstack, vstack, join
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from spectral_cube import SpectralCube
import scipy 
from reproject import reproject_interp

import astropy.constants as ac
import astropy.units as au
from astropy import stats
from astrodendro.analysis import PPStatistic
from astrodendro import Dendrogram, pp_catalog
from astropy.wcs import WCS
from collections import Counter

from astropy.io import fits
import matplotlib as mpl
import pyregion
import aplpy
import math
import os
import pickle
from tqdm.auto import tqdm

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D


def get_scaling(data):
    data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    return(data)

def plot_sigmacomp(hdus, props_all, regions, regionID, min_value_sig, outdir='./../../analysis/dendro/indexmaps_figs_sigmacomp/'):

    """Make plot"""
    
    hdu = hdus['hst07_hdus_smooth'][regionID]
    mask_muse = hdus['musmask_hdus'][regionID]
    
    ra = regions['ra'][regionID].value
    dec = regions['dec'][regionID].value
    width = regions['width'][regionID].value
    height = regions['height'][regionID].value
    radius = max([width,height])
    center = [ra, dec, radius, radius]
    
    min_value_sig_n = len(min_value_sig)
    min_value_sig_str = ['%ssig' %min_value_sig[i] for i in range(min_value_sig_n)]

    props = []
    for key in min_value_sig_str:
        props += [props_all[key][regionID]]
    
    fig = plt.figure(figsize=(10,10))
    ax = ['']*min_value_sig_n
    
    hdu_ = hdu.copy()
    hdu_.data = get_scaling(hdu_.data)
    
    for i, key in enumerate(min_value_sig_str):
                
        ax[i] = aplpy.FITSFigure(hdu_, figure=fig, subplot=(1,len(min_value_sig_str),i+1))
        
        minmax = np.nanpercentile(hdu_.data, [0.1,99.99])
        ax[i].show_colorscale(cmap='gist_gray_r', vmin=minmax[0], vmax=minmax[-1], stretch='log')
        ax[i].recenter(center[0], center[1], width=center[2], height=center[3])
        ax[i].tick_labels.hide()
        ax[i].axis_labels.hide()
        ax[i].ticks.set_color('black')
        ax[i].set_nan_color('white')

        if np.ma.is_masked(props_all[key]['_idx'][regionID]): 
            ax[i].add_label(0.5,0.5,'NO FIT',c='r',relative=True)
            continue
                        
        ra_struc = props_all[key]['ra_cen'].quantity[regionID].value
        dec_struc = props_all[key]['dec_cen'].quantity[regionID].value
    
        ax[i].show_contour(hdus['indexmap_trunk_hdus_%s' %min_value_sig_str[i]][regionID], colors='C0', levels=[-1], linestyles='-')
        
        ax[i].show_ellipses(xw=ra_struc, yw=dec_struc, 
                 width=props_all[key]['radius_trunk'].quantity[regionID].to('deg').value*2,
                 height=props_all[key]['radius_trunk'].quantity[regionID].to('deg').value*2,
                 linestyle='--', edgecolor='C1', linewidth=2)
    
        #Ellipses
        ax[i].show_ellipses(xw=ra_struc, yw=dec_struc, 
                         width=props_all[key]['major_sigma'].quantity[regionID].to('deg').value*2,#sigma
                         height=props_all[key]['minor_sigma'].quantity[regionID].to('deg').value*2,#sigma
                         angle=props_all[key]['position_angle'].quantity[regionID].to('deg').value,
                         linestyle='-', edgecolor='C3', linewidth=2)

        ax[i].show_ellipses(xw=ra_struc, yw=dec_struc, 
                         width=props_all[key]['major_sigma'].quantity[regionID].to('deg').value*4.292,#FWTM
                         height=props_all[key]['minor_sigma'].quantity[regionID].to('deg').value*4.292,#FWTM
                         angle=props_all[key]['position_angle'].quantity[regionID].to('deg').value,
                         linestyle='--', edgecolor='C3', linewidth=2)  

        #Mean of ellipses
        ax[i].show_ellipses(xw=ra_struc, yw=dec_struc, 
                         width=props_all[key]['mean_sigma'].quantity[regionID].to('deg').value*2,#sigma
                         height=props_all[key]['mean_sigma'].quantity[regionID].to('deg').value*2,#sigma
                         linestyle=':', edgecolor='C3', linewidth=2)

        ax[i].show_ellipses(xw=ra_struc, yw=dec_struc, 
                         width=props_all[key]['mean_sigma'].quantity[regionID].to('deg').value*4.292,#FWTM
                         height=props_all[key]['mean_sigma'].quantity[regionID].to('deg').value*4.292,#FWTM
                         linestyle=':', edgecolor='C3', linewidth=2) 

        ax[i].add_label(0.02,0.95,'min_value_sig=%isigma' %min_value_sig[i],relative=True,ha='left')
        
    if mask_muse!=None:
        for i, key in enumerate(min_value_sig_str):
            ax[i].show_contour(mask_muse, colors='black', levels=[1], linestyles='-')
            
        fig.get_axes()[0].plot([0,0],[0,0],ls='-',c='black',label='MUSE mask')
        
    
    fig.tight_layout()
    fig.savefig('%s/dendro_regionoutput_%i.pdf' %(outdir, regionID), bbox_inches='tight', dpi=300)
    plt.close('all')
    
    return()

def plot_sigmacomp_decorator(hdus, props_all, regions, min_value_sig, outdir='./../../analysis/dendro/indexmaps_figs_sigmacomp/'):

    n = len(hdus[list(hdus.keys())[0]])
    for regionID in tqdm(range(n), desc='Dendrogram (with sigma ranges)', position=0):
        
        # if regionID > 1: 
        #     continue 

        plot_sigmacomp(hdus, props_all, regions, regionID, min_value_sig)
    
    return()


