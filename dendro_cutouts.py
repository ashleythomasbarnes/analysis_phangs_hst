import warnings
warnings.filterwarnings('ignore')

import numpy as np
from astropy.table import Table
import astropy.constants as ac
import astropy.units as au
from spectral_cube import SpectralCube

from astropy.io import fits
import pyregion
import os
import pickle
from tqdm.auto import tqdm 
import reproject

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

import gc

def get_regions(regions_file, hdu):
    """
    Extracts region properties from a regions file and converts coordinates to image coordinates.

    Args:
        regions_file (str): Path to the regions file.
        hdu (astropy.io.fits.ImageHDU): Image HDU containing the data, or list of HDUs 

    Returns:
        dict: Dictionary containing region properties.

    """
    # Open the regions file and convert the coordinates to image coordinates

    try:
        header = hdu.header
    except:
        header = hdu[0].header

    print(f'[INFO] [get_regions] Opening regions file (this may take a min)... ')
    regions = pyregion.open(regions_file).as_imagecoord(header)
    n = len(regions)

    # Initialize empty arrays for storing region properties
    ra = np.empty(n) * np.nan * au.deg
    dec = np.empty(n) * np.nan * au.deg
    width = np.empty(n) * np.nan * au.deg
    height = np.empty(n) * np.nan * au.deg
    position = [''] * n

    print(f'[INFO] [get_regions] Getting info for {n} regions...')

    for i, region in enumerate(regions):
        # Extract the region properties and store them in the respective arrays
        ra[i] = float(region.params[0].text) * au.deg
        dec[i] = float(region.params[1].text) * au.deg
        width[i] = region.params[2].degree * au.deg
        height[i] = region.params[3].degree * au.deg

    # Create a SkyCoord object for the positions
    position = SkyCoord(ra=ra, dec=dec, frame='icrs')

    # Return a dictionary containing the region properties
    return {'ra': ra, 'dec': dec, 'width': width, 'height': height, 'position': position}


def get_croppeddata(hdu, i, regions, square=True):
    """
    Function to crop data from an HDU object based on the provided region properties.

    Parameters:
        hdu (astropy.io.fits.ImageHDU): The input HDU object.
        i (int): Index of the region in the regions dictionary.
        regions (dict): Dictionary containing region properties.
        square (bool): Determines whether to create a square cutout (default: True).

    Returns:
        astropy.io.fits.ImageHDU: The cropped HDU object.
    """

    # print('now', i, regions['position'][i])
    
    hdu_crop = hdu.copy()  # Create a copy of the input HDU object
    del hdu  # Delete the original HDU object to free up memory
    _ = gc.collect()  # Perform garbage collection

    wcs = WCS(hdu_crop)  # Create a WCS object from the HDU header

    if square:
        radius = max([regions['width'][i], regions['height'][i]]) * 1.3  # Calculate the radius for the square cutout with a buffer
        cutout = Cutout2D(hdu_crop.data, regions['position'][i], radius, wcs=wcs)  # Create a square cutout
    else:
        cutout = Cutout2D(hdu_crop.data, regions['position'][i], [regions['width'][i], regions['height'][i]], wcs=wcs)  # Create a rectangular cutout

    hdu_crop.data = cutout.data  # Update the data in the cropped HDU object
    hdu_crop.header.update(cutout.wcs.to_header())  # Update the header of the cropped HDU object with the cutout's WCS information

    del cutout  # Delete the cutout to free up memory
    _ = gc.collect()  # Perform garbage collection

    return hdu_crop  # Return the cropped HDU object


def get_croppeddata_parallel(hdu, regions, n_jobs=-1):
    """
    Function to crop data from an HDU object based on the provided region properties in parallel.

    Parameters:
        hdu (astropy.io.fits.ImageHDU): The input HDU object.
        regions (dict): Dictionary containing region properties.
        n_jobs (int): Number of parallel jobs to run (-1 for using all available processors, default: -1).

    Returns:
        list: List of cropped HDU objects in the same order as the input regions.
    """

    n = len(regions['ra'])  # Get the number of regions

    # Use parallel processing to crop the data from multiple regions with progress bar
    from joblib import Parallel, delayed
    with tqdm(total=n, desc="Cropping regions") as pbar:
        def crop_data(i):
            cropped_hdu = get_croppeddata(hdu, i, regions)
            pbar.update(1)
            return i, cropped_hdu

        results = Parallel(n_jobs=n_jobs, backend='threading')(delayed(crop_data)(i) for i in range(n))

    # Sort the results based on the original order of the regions
    sorted_results = sorted(results, key=lambda x: x[0])
    hdus = [result[1] for result in sorted_results]

    return hdus  # Return the list of cropped HDU objects in the same order as the input regions