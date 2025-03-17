import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from regions import Regions, RectangleSkyRegion, CircleSkyRegion, EllipseSkyRegion
import astropy.units as u
from astropy.io import fits

def get_croppeddata(hdu, ra, dec, width, height):
    """
    Function to crop data from an HDU object based on the provided region properties.

    Parameters:
        hdu (astropy.io.fits.ImageHDU): The input HDU object.
        i (int): Index of the region in the regions dictionary.
        regions (dict): Dictionary containing region properties.

    Returns:
        astropy.io.fits.ImageHDU: The cropped HDU object.
    """
    
    hdu_crop = hdu.copy()  # Copy the HDU object to avoid modifying the original data
    wcs = WCS(hdu_crop)  # Create a WCS object from the HDU header

    centre = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')  # Create a SkyCoord object for the center of the region
    cutout = Cutout2D(hdu_crop.data, centre, [width*u.deg*1.2, height*u.deg*1.2], wcs=wcs)  # Create a rectangular cutout
    hdu_crop = fits.PrimaryHDU(cutout.data, cutout.wcs.to_header())
    
    del hdu 
    del cutout  # Delete the cutout to free up memory

    return hdu_crop  # Return the cropped HDU object


def get_regions(regions_file):
    """
    Extracts region properties from a regions file and converts coordinates to image coordinates.

    Args:
        regions_file (str): Path to the regions file.

    Returns:
        dict: Dictionary containing region properties.

    """
    print(f'[INFO] [get_regions] Opening regions file (this may take a min)... ')
    # regions = pyregion.open(regions_file).as_imagecoord(header)
    regions = Regions.read(regions_file)
    n = len(regions)

    # Initialize empty arrays for storing region properties
    ra = np.empty(n) * np.nan * u.deg
    dec = np.empty(n) * np.nan * u.deg
    radius = np.empty(n) * np.nan * u.deg

    print(f'[INFO] [get_regions] Getting info for {n} regions...')

    for i, region in enumerate(regions):
        # Extract the region properties and store them in the respective arrays
        ra[i] = float(region.center.ra.deg) * u.deg
        dec[i] = float(region.center.dec.deg) * u.deg
        radius[i] = (region.width + region.height) / 2

    # Create a SkyCoord object for the positions
    position = SkyCoord(ra=ra, dec=dec, frame='icrs')

    # Return a dictionary containing the region properties
    return ra.value, dec.value, radius.value, radius.value, n


def get_hist(data, bins=None, nbins=50, logbins=False, norm=True, cum=False):

    """Get histogram

    Parameters
    ----------
    data : np.array
        input data
    bins : list
        input bin edges for histogram calculaiton; default=''
    nbins : int
        number of bins to determine if bins is not given; defult=50
    logbins : bool
        logarithmically spaced bins if bins is not given
    norm : bool
        normalise such that max is equal to unity; default=True
    cum : bool
        cumulative distorbution; otherwise probability distorbution
    Returns
    -------
    bins : list
        bin edges for histogram calculaiton
    bin_cent : np.array
        bin centres, for easy plotting in matplotlib
    hist : np.array
        histogram data for each bin centre
    """

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