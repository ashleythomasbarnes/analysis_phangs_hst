from astropy.wcs import WCS
from astropy.io import fits
from tqdm.auto import tqdm
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion, EllipseSkyRegion, Regions
import numpy as np

def get_ds9regions_circ(ra, dec, radius):
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
    region_sky = [CircleSkyRegion(center, rad) for center, rad in zip(center_sky, radius)]
    
    # Convert the list into Regions object
    region_sky = Regions(region_sky)
    
    return region_sky


def get_ds9regions_ellip(ra, dec, height, width, angle):
    """
    Create a list of elliptical regions for DS9 using given RA, Dec, height, width, and angle.

    Parameters:
    - ra (array): Right Ascension coordinates of the ellipse centers in degrees.
    - dec (array): Declination coordinates of the ellipse centers in degrees.
    - height (array): The height (or major axis) of the ellipses.
    - width (array): The width (or minor axis) of the ellipses.
    - angle (array): The orientation angles of the ellipses.

    Returns:
    - Regions: A list of elliptical regions for DS9.
    """

    # Convert the given RA and Dec into SkyCoord object
    center_sky = SkyCoord(ra, dec, unit='deg', frame='fk5')
    
    # Create elliptical regions using provided parameters
    region_sky = [EllipseSkyRegion(center, h, w, a)
                  for center, h, w, a in zip(center_sky, height, width, angle)]
    
    # Convert the list into Regions object
    region_sky = Regions(region_sky)
    
    return region_sky


def get_ds9regions_circ_decor(props_all, outputfile='ds9.reg'):
    """
    Generate a DS9 region file for circular regions using properties from an input dictionary.

    Parameters:
    - props_all (dict): A dictionary containing the properties of regions such as ra_cen, dec_cen, and radius_trunkclose.
    - outputfile (str): The filename where the DS9 regions will be saved. Default is 'ds9.reg'.

    Returns:
    - None
    """

    # Create circular regions using provided properties
    region_sky = get_ds9regions_circ(props_all['ra_cen'], props_all['dec_cen'], props_all['radius_trunkclose'].quantity)
    
    # Write the regions to a DS9 region file
    region_sky.write(outputfile, format='ds9', overwrite=True)

    print(f'[INFO] [get_ds9regions_circ_decor] Saved regions to {outputfile}')


def get_ds9regions_ellip_decor(props_all, outputfile='ds9.reg'):
    """
    Generate a DS9 region file for elliptical regions using properties from an input dictionary.

    Parameters:
    - props_all (dict): A dictionary containing the properties of regions such as ra_cen, dec_cen, major_sigma, minor_sigma, and position_angle.
    - outputfile (str): The filename where the DS9 regions will be saved. Default is 'ds9.reg'.

    Returns:
    - None
    """

    # Create elliptical regions using provided properties
    region_sky = get_ds9regions_ellip(props_all['ra_cen'], props_all['dec_cen'], 
                                      props_all['major_sigma'].quantity*2, 
                                      props_all['minor_sigma'].quantity*2, 
                                      props_all['position_angle'].quantity)
    
    # Write the regions to a DS9 region file
    region_sky.write(outputfile, format='ds9', overwrite=True)

    print(f'[INFO] [get_ds9regions_ellip_decor] Saved regions to {outputfile}')


def get_hdumask(hdu_ha, hdu_mask, outputfile='tmp.fits'):
    """
    Generate a mask HDU from a list of HDUs, overlaying it on a reference H-alpha HDU.
    
    Parameters:
    - hdus (dict): Dictionary containing HDUs of interest, indexed by 'indexmap_trunk_hdus_3sig'.
    - hdu_ha (HDU): The reference H-alpha HDU on which the masks will be overlaid.
    - outputfile (str, optional): Name of the output FITS file for the resulting mask. Default is 'tmp.fits'.
    
    Returns:
    - None
    """
    
    # Set the shape of the reference H-alpha data and initialize its data array to -1
    shape_ha = hdu_ha.data.shape
    hdu_ha.data = np.ones(shape_ha) * -1
    wcs_ha = WCS(hdu_ha.header)

    # Loop over the HDUs provided in the 'indexmap_trunk_hdus_3sig' key
    for hdu_id in tqdm(range(len(hdu_mask)), desc='Masking regions'):

        # Extract the current HDU mask and its WCS
        hdu_mask_ = hdu_mask[hdu_id]
        data_mask = (hdu_mask_.data != -1)  # Convert to a Boolean mask
        shape_mask = data_mask.shape
        n_mask = np.sum(data_mask*1)  # Count of non-negative values in the mask
        wcs_mask = WCS(hdu_mask_.header)

        # Overlay the mask onto the reference H-alpha HDU
        for i in range(shape_mask[0]):
            for j in range(shape_mask[1]):
                if data_mask[i, j]:

                    # Convert the pixel coordinate from mask HDU to a sky coordinate
                    sky = wcs_mask.pixel_to_world(j, i)
                    # Convert the sky coordinate to pixel coordinate in the reference H-alpha HDU
                    pix = wcs_ha.world_to_pixel(sky)
                    pix = np.array(pix, dtype=np.float64)
                    x = round(pix[1])
                    y = round(pix[0])

                    # Update the reference H-alpha HDU with the ID of the mask HDU
                    hdu_ha.data[x, y] = hdu_id

    # hdu_ha.data = np.int16(hdu_ha.data)
    hdu_ha.data = hdu_ha.data.astype('int16')

    # Write the resulting mask to an output file
    hdu_ha.writeto(outputfile, overwrite=True)
