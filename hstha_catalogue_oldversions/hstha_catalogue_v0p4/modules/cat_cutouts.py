import sys
sys.path.append('../modules/')
from cat_imports import *

def get_ds9region_sq(ra, dec, radius):
    """
    Create a list of circular regions for DS9 using given RA, Dec, and radii.

    Parameters:
    - ra (array): Right Ascension coordinates of the circle centers in degrees.
    - dec (array): Declination coordinates of the circle centers in degrees.
    - radius (array): Radii of the circles.

    Returns:
    - Regions: A list of circular regions for DS9.
    """

    # Convert the given RA and Dec into SkyCoord object
    center_sky = SkyCoord(ra, dec, unit='deg', frame='fk5')
    
    # Create circular regions using provided center coordinates and radii
    region_sky = [RectangleSkyRegion(center, rad, rad) for center, rad in zip(center_sky, radius)]
    
    # Convert the list into Regions object
    region_sky = Regions(region_sky)
    
    return region_sky

def get_ds9regions_all(props_all, outputfile='ds9.reg'):
    """
    Generate a DS9 region file for circular regions using properties from an input dictionary.

    Parameters:
    - props_all (dict): A dictionary containing the properties of regions such as ra_cen, dec_cen, and radius_trunkclose.
    - outputfile (str): The filename where the DS9 regions will be saved. Default is 'ds9.reg'.

    Returns:
    - None
    """

    # Create circular regions using provided properties
    region_sky = get_ds9region_sq(props_all['cen_ra'], props_all['cen_dec'], props_all['region_circ_rad'].quantity*4)
    
    # Write the regions to a DS9 region file
    region_sky.write(outputfile, format='ds9', overwrite=True)

    return (region_sky)


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
    ra = np.empty(n) * np.nan * au.deg
    dec = np.empty(n) * np.nan * au.deg
    width = np.empty(n) * np.nan * au.deg
    height = np.empty(n) * np.nan * au.deg
    position = [''] * n

    print(f'[INFO] [get_regions] Getting info for {n} regions...')

    for i, region in enumerate(regions):
        # Extract the region properties and store them in the respective arrays
        ra[i] = float(region.center.ra.deg) * au.deg
        dec[i] = float(region.center.dec.deg) * au.deg
        width[i] = region.width
        height[i] = region.height

    # Create a SkyCoord object for the positions
    position = SkyCoord(ra=ra, dec=dec, frame='icrs')

    # Return a dictionary containing the region properties
    return {'ra': ra, 'dec': dec, 'width': width, 'height': height, 'position': position}


def get_croppeddata(hdu, i, regions):
    """
    Function to crop data from an HDU object based on the provided region properties.

    Parameters:
        hdu (astropy.io.fits.ImageHDU): The input HDU object.
        i (int): Index of the region in the regions dictionary.
        regions (dict): Dictionary containing region properties.

    Returns:
        astropy.io.fits.ImageHDU: The cropped HDU object.
    """

    wcs = WCS(hdu)  # Create a WCS object from the HDU header

    try: 
        cutout = Cutout2D(hdu.data, regions['position'][i], [regions['width'][i], regions['height'][i]], wcs=wcs)  # Create a rectangular cutout
        hdu_crop = fits.PrimaryHDU(cutout.data, cutout.wcs.to_header())
    except:
        # print('[INFO] [get_croppeddata] NoOverlapError: filling with nans for Region %i' %i)
        cutout = np.ones([5,5])*np.nan
        hdu_crop = fits.PrimaryHDU(cutout.data, hdu.header)
    
    del hdu 
    del cutout  # Delete the cutout to free up memory
    _ = gc.collect()  # Perform garbage collection

    return hdu_crop  # Return the cropped HDU object

def get_croppeddata_all(hdu, regions):
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
    hdus = []
    for i in tqdm(range(n), desc="Cropping regions"):
        hdus += [get_croppeddata(hdu, i, regions)]
                        
    return hdus 

def get_croppeddata_all_parallel(hdu, regions, n_jobs=-1):
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