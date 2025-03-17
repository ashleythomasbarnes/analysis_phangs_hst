import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table, join, vstack
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from reproject import reproject_interp

import astropy.units as au
from astropy import stats
from astrodendro import Dendrogram, pp_catalog
from astropy.wcs import WCS
from astropy.table import Column

from astropy.io import fits
import aplpy
from tqdm.auto import tqdm

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.ndimage import binary_dilation

from astropy.wcs import WCS

def get_pcperarcsec(sampletable):

    dist = sampletable['dist'].quantity[0]
    pcperarcsec = dist.to('pc').value/206265
    
    return(pcperarcsec)

def get_trunkindexmap(dendro):
    """Get index map of trunk idx"""
    
    # print('[INFO] Getting indexmap of dendrogram trunks')
    mask = 0
    for i in range(len(dendro.trunk)):

        trunk = dendro.trunk[i]
        idx = trunk.idx
        mask = trunk.get_mask()
        mask = mask*1.0
        mask[mask==0] = np.nan
        mask[mask==1] = idx

        if i==0: 
            indexmap_trunk = mask
        else: 
            indexmap_trunk = np.nanmean(np.array([indexmap_trunk,mask]), axis=0)  

    indexmap_trunk[np.isnan(indexmap_trunk)] = -1
    
    return(indexmap_trunk)

def get_circrad(data, pixsize):
    """Returns circularised radius in arcsec"""
    area = np.nansum((data > -1)*1)
    radius = (area/np.pi)**0.5
    radius = radius * pixsize.to('arcsec')
    return(radius)

def get_maskeddata(muscat_hdu, hst_hdu, muscat_id): 
    
    hst_masked_hdu = hst_hdu.copy() #only ID in mask
    hst_masked_ones_hdu = hst_hdu.copy() #only ID in mask
    hst_maskedall_hdu = hst_hdu.copy() #all MUSE regions out of mask - for noise measure
    
    #get mask 
    shape = muscat_hdu.data.shape
    # muscat_id = muscat_hdu.data[int(shape[0]/2),int(shape[0]/2)]
    
    mask1 = muscat_hdu.data!=muscat_id #if catalouge is not ID
    mask2 = np.isnan(muscat_hdu.data) #if catalouge is not a number
    
    muscat_hdu1 = muscat_hdu.copy()
    muscat_hdu2 = muscat_hdu.copy()
    
    muscat_hdu1.data[mask1] = np.nan
    muscat_hdu2.data[mask2] = 1
    muscat_hdu2.data[~mask2] = np.nan

    #regrid
    data1 = reproject_interp(muscat_hdu1, hst_hdu.header, return_footprint=False, order='nearest-neighbor')
    data2 = reproject_interp(muscat_hdu2, hst_hdu.header, return_footprint=False, order='nearest-neighbor')

    #mask 
    hst_masked_hdu.data[np.isnan(data1)] = np.nan
    hst_masked_ones_hdu.data[np.isnan(data1)] = 1
    hst_maskedall_hdu.data[np.isnan(data2)] = np.nan
    
    return(hst_masked_hdu, hst_masked_ones_hdu, hst_maskedall_hdu, muscat_id)

def get_maskedhdus(hdus, regions, muscat_regionIDs, hstha_hdu_name='hstha_hdu'):

    hstha_hdu_smooth_arr = []
    hstha_masked_hdu_arr = []
    hst_masked_ones_hdu_arr = []
    hstha_maskedall_hdu_arr = []
    muscat_id_arr = []
    mask_muse_arr = []

    print(f'[INFO] [get_maskedhdus] Getting HST maps masked by MUSE catalouge...')
    for i in tqdm(range(len(regions['ra'])), desc='Masking regions', position=0):
    # for i in tqdm(range(50), desc='Masking regions', position=0):

        # Load
        muscat_hdu = hdus['muscat_hdu'][i]
        hstha_hdu_smooth = hdus[hstha_hdu_name][i].copy()
        muscat_regionID = muscat_regionIDs[i]

        # Smooth the data
        kernel = Gaussian2DKernel(x_stddev=0.5)
        hstha_hdu_smooth.data = convolve(hstha_hdu_smooth.data, kernel)

        # Mask the data
        output = get_maskeddata(muscat_hdu, hstha_hdu_smooth, muscat_regionID)
        hstha_masked_hdu, hst_masked_ones_hdu, hstha_maskedall_hdu, muscat_id = output

        # Make HDU
        mask_muse = fits.PrimaryHDU(~np.isnan(hstha_masked_hdu.data) * 1, hstha_hdu_smooth.header)

        # Extract the results in the correct order
        hstha_hdu_smooth_arr += [hstha_hdu_smooth]
        hstha_masked_hdu_arr += [hstha_masked_hdu]
        hst_masked_ones_hdu_arr += [hst_masked_ones_hdu]
        hstha_maskedall_hdu_arr += [hstha_maskedall_hdu]
        muscat_id_arr += [muscat_regionID]
        mask_muse_arr += [mask_muse]

    # Assign the processed data to the corresponding keys in the hdus dictionary
    hdus['%s_smooth' % hstha_hdu_name] = hstha_hdu_smooth_arr
    hdus['%s_smooth_masked' % hstha_hdu_name] = hstha_masked_hdu_arr
    hdus['%s_smooth_masked_ones' % hstha_hdu_name] = hst_masked_ones_hdu_arr
    hdus['%s_smooth_maskedall' % hstha_hdu_name] = hstha_maskedall_hdu_arr
    hdus['musmask_hdu'] = mask_muse_arr

    # Return the modified hdus dictionary and muscat_id
    return hdus

def all_nan_check(hdus, regionID):
    
    data = hdus['hstha_hdu_smooth'][regionID].data
    shape = data.shape

    # xmin = shape[0]/2 - shape[0]/4
    # xmax = shape[0]/2 + shape[0]/4
    # ymin = shape[1]/2 - shape[1]/4
    # ymax = shape[1]/2 + shape[1]/4
    # data = data[int(xmin):int(xmax),int(ymin):int(ymax)]
    
    if np.isnan(data).all():
        # print(f'[INFO] [props_all] [regionID={regionID}] All nan values in map...')
        return (True)
    else:
        return (False)

def get_dedro(hdu, hdu_outmask, hdu_nomask, sampletable, min_npix=9, min_value_sig=3, min_delta_sig=3):
    """
    Perform dendrogram analysis on the input data and extract properties of the identified structures.
    
    Parameters:
        hdu (HDU): Input data HDU.
        hdu_outmask (HDU): Output mask HDU.
        hdu_nomask (HDU): Non-masked data HDU.
        min_npix (int, optional): Minimum number of pixels within a structure. Defaults to 9.
        min_value_sig (int, optional): Minimum value within a structure in terms of standard deviation. Defaults to 3.
        min_delta_sig (int, optional): Minimum values between levels within a structure in terms of standard deviation. Defaults to 3.
    
    Returns:
        tuple: A tuple containing the following:
            - dendro (Dendrogram): The computed dendrogram.
            - props (Table): Properties of the identified structure.
            - indexmap_trunk_hdu (HDU): Index map of the trunk structure.
            - indexmap_trunk_close_hdu (HDU): Index map of the trunk structure with closing operation applied.
    """

    # Replace variable names with appropriate names
    wcs = WCS(hdu.header)
    header = hdu.header

    data = hdu.data
    data_nomask = hdu_nomask.data
    data_outmask = hdu_outmask.data

    # Calculate statistics for dendrogram
    std = stats.mad_std(data_outmask, ignore_nan=True)  # Get noise
    std = stats.mad_std(data_outmask[data_outmask<20*std], ignore_nan=True)  # Get noise below threshold

    try:
        pixsize = np.array([np.abs(header['CDELT1']), np.abs(header['CDELT2'])]).mean() * au.degree
    except: 
        pixsize = np.array([np.abs(header['CD1_1']), np.abs(header['CD2_2'])]).mean() * au.degree

    if pixsize.value==1: 
        if 'CD1_1' in np.array(header.cards)[:,0]: 
            pixsize = np.array([np.abs(header['CD1_1']), np.abs(header['CD2_2'])]).mean() * au.degree
        elif 'PC1_1' in np.array(header.cards)[:,0]:
            pixsize = np.array([np.abs(header['PC1_1']), np.abs(header['PC2_2'])]).mean() * au.degree

    bmaj = bmin = 0.05 * au.arcsec  # Dummy values 

    min_value = std * min_value_sig  # Minimum value within structure
    min_delta = std * min_delta_sig  # Minimum values between levels within structure (sort of)

    # Running the dendrogram
    dendro = Dendrogram.compute(data,
                                min_delta=min_delta,
                                min_value=min_value,
                                min_npix=min_npix,
                                wcs=wcs)

    # Provide metadata for table output
    metadata = {}
    metadata['data_unit'] = au.Jy / au.beam  # Dummy unit
    metadata['spatial_scale'] = pixsize.to('arcsec')
    metadata['beam_major'] = bmaj.to('arcsec')
    metadata['beam_minor'] = bmin.to('arcsec')

    if len(np.unique(dendro.index_map)) != 1: 
        props = pp_catalog(dendro, metadata, verbose=False)  # Get table
    else: 
        return([None] * 4)

    indexmap = dendro.index_map  # Get index map
    indexmap_hdu = fits.PrimaryHDU(indexmap, header)
    indexmap_hdu.data = np.array(indexmap_hdu.data, dtype=float)
    indexmap_hdu.data[indexmap_hdu.data == -1] = np.nan

    indexmap_trunk = get_trunkindexmap(dendro)
    indexmap_trunk_hdu = fits.PrimaryHDU(indexmap_trunk, header)
    indexmap_trunk_hdu.data = np.array(indexmap_trunk_hdu.data, dtype=float)

    props_trunk = pp_catalog(dendro.trunk, metadata, verbose=False)  # Get table
    
    # Get the ID of the structure with maximum flux
    arg = np.argmax(props_trunk['flux'])
    midpix_id = props_trunk['_idx'][arg]

    # Remove all other structures from the catalogue
    indexmap_trunk_hdu.data[indexmap_trunk_hdu.data != midpix_id] = -1
    props = props[props['_idx'] == midpix_id]

    # Get single trunk index map
    idx = np.unique(indexmap_trunk_hdu.data[indexmap_trunk_hdu.data != -1])
    mask = indexmap_trunk_hdu.data != -1
    mask1 = scipy.ndimage.binary_closing(mask, iterations=10) * 1
    mask1[mask1 != 1.0] = -1.
    mask1[mask1 == 1.0] = idx
    indexmap_trunk_close = np.copy(mask1)
    indexmap_trunk_close_hdu = fits.PrimaryHDU(indexmap_trunk_close, indexmap_trunk_hdu.header)

    radius_trunk = get_circrad(indexmap_trunk_hdu.data, pixsize)
    radius_trunkclose = get_circrad(indexmap_trunk_close_hdu.data, pixsize)

    # Add properties to table    
    props['radius_trunk'] = radius_trunk
    props['radius_trunkclose'] = radius_trunkclose
    props['major_fwtm'] = props['major_sigma'] * 4.292
    props['minor_fwtm'] = props['minor_sigma'] * 4.292
    props['mean_fwtm'] = np.nanmean([props['major_fwtm'], props['minor_fwtm']]) * props['major_fwtm'].unit
    props['mean_hwtm'] = props['mean_fwtm'] / 2.
    
    props['major_fwhm'] = props['major_sigma'] * 2*np.sqrt(2*np.log(2))
    props['minor_fwhm'] = props['minor_sigma'] * 2*np.sqrt(2*np.log(2))
    props['mean_fwhm'] = np.nanmean([props['major_fwhm'], props['minor_fwhm']]) * props['major_fwhm'].unit
    props['mean_hwhm'] = props['mean_fwhm'] / 2.

    props['min_npix'] = min_npix
    props['min_value_sig'] = min_value_sig
    props['min_delta_sig'] = min_delta_sig
    
    # Convert to parcsec
    pcperarcsec = get_pcperarcsec(sampletable)
    props['radius_trunk_pc'] = props['radius_trunk'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['radius_trunkclose_pc'] = props['radius_trunkclose'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['major_fwtm_pc'] = props['major_fwtm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['minor_fwtm_pc'] = props['minor_fwtm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['mean_fwtm_pc'] = props['mean_fwtm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['mean_hwtm_pc'] = props['mean_fwtm_pc'] / 2.

    props['major_fwhm_pc'] = props['major_fwhm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['minor_fwhm_pc'] = props['minor_fwhm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['mean_fwhm_pc'] = props['mean_fwhm'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    props['mean_hwhm_pc'] = props['mean_fwhm_pc'] / 2.

    # Get max x,y positions
    mask = indexmap_trunk_close_hdu.data.copy() == -1
    data[mask] = 0
    y_max, x_max = np.where(data==np.nanmax(data))

    ax = aplpy.FITSFigure(hdu)
    ra_struc = ax.pixel2world(props['x_cen'], props['y_cen'])[0][0]
    dec_struc = ax.pixel2world(props['x_cen'], props['y_cen'])[1][0]
    ra_max = ax.pixel2world(x_max, y_max)[0][0]
    dec_max = ax.pixel2world(x_max, y_max)[1][0]
    plt.close('all')

    props['ra_cen'] = ra_struc * au.deg
    props['dec_cen'] = dec_struc * au.deg
    props.rename_column('radius', 'mean_sigma')
    props['mean_sigma_pc'] = props['mean_sigma'].quantity[0].to('arcsec').value * pcperarcsec * au.pc
    
    props['x_max'] = x_max *props['x_cen'].unit
    props['y_max'] = y_max *props['x_cen'].unit
    props['ra_max'] = ra_max *au.deg
    props['dec_max'] = dec_max *au.deg

    flux = np.nansum(hdu.data[indexmap_trunk_close_hdu.data != -1])
    props['flux'] = flux
    props['flux'].unit = au.erg / au.s / au.cm**2

    props = Table(props, masked=True, copy=False)

    return dendro, props, indexmap_trunk_hdu, indexmap_trunk_close_hdu


def get_dedro_all(hdus, sampletable, muscat_regionIDs, min_npix=9, min_value_sig=3, min_delta_sig=3):
    
    props_all = []
    indexmap_trunk_hdu_all = []
    indexmap_trunk_close_hdu_all = []
    first_no_found_run = 0

    for i in tqdm(range(len(muscat_regionIDs)), desc='Dendrogram', position=0):

        # if i >10: 
        #     continue
        
        muscat_regionID = muscat_regionIDs[i]

        if not all_nan_check(hdus, i): 

            output = get_dedro(hdus['hstha_hdu_smooth_masked'][i],
                                hdus['hstha_hdu_smooth_maskedall'][i], 
                                hdus['hstha_hdu_smooth'][i],
                                sampletable=sampletable,
                                min_npix=min_npix, 
                                min_value_sig=min_value_sig, 
                                min_delta_sig=min_delta_sig)

            dendro, props, indexmap_trunk_hdu, indexmap_trunk_close_hdu = output

            try:
                # Change index to MUSE indexing in indexmaps
                indexmap_trunk_hdu.data[indexmap_trunk_hdu.data!=-1] = muscat_regionID
                indexmap_trunk_close_hdu.data[indexmap_trunk_close_hdu.data!=-1] = muscat_regionID
            except:
                print(f'[INFO] [get_dedro] [regionID={muscat_regionID}] No dendro found...')
                None
        else: 
            print(f'[INFO] [get_dedro] [regionID={muscat_regionID}] No dendro found, all nan values...')
            props = None
        
        # Replace values when no region is found due to all nans or other
        if props is None: 

            # Only run this one once...
            # Replace values with some dummy dendro table so we have all the correct columns, will be masked later anyway... 
            if first_no_found_run == 0: 
                first_no_found_run = 1
                
                print(f'[INFO] [get_dedro] Running dummy dendro...')
                j=0
                while props is None: 
                    j+=1 
                    output_ = get_dedro(hdus['hstha_hdu_smooth_masked'][j],
                                        hdus['hstha_hdu_smooth_maskedall'][j], 
                                        hdus['hstha_hdu_smooth'][j],
                                        sampletable=sampletable,
                                        min_npix=0, 
                                        min_value_sig=-100, 
                                        min_delta_sig=-100)
                    dendro, props_, indexmap_trunk_hdu, indexmap_trunk_close_hdu = output_
                    if props_ is not None:
                        props = props_.copy()
                    else:
                        props = None
                    
            else:           
                dendro, props, indexmap_trunk_hdu, indexmap_trunk_close_hdu = output_
                props = props_.copy()

            for colname in props.colnames: 
                props[colname] = np.ma.masked

            # Blanks array in index map if nothing is found 
            shape = indexmap_trunk_hdu.data.shape
            header = indexmap_trunk_hdu.header
            indexmap_trunk_hdu = fits.PrimaryHDU(np.ones(shape)*-1, header)
            indexmap_trunk_close_hdu = fits.PrimaryHDU(np.ones(shape)*-1, header)
        
        # Add good and masked to properties 
        props.add_column(Column(muscat_regionID, name='region_ID'), index=0)  
        props.add_column(Column(i, name='hstcat_region_ID'), index=1)

        # Edge flag
        flag_edge = np.isnan(hdus['hstha_hdu_smooth_masked_ones'][i].data).any()*1
        props.add_column(Column(flag_edge, name='flag_edge_hst')) 

        props_all += [props]
        indexmap_trunk_hdu_all += [indexmap_trunk_hdu]
        indexmap_trunk_close_hdu_all += [indexmap_trunk_close_hdu]

    props_all = vstack(props_all)
    hdus['indexmap_trunk_hdu'] = indexmap_trunk_hdu_all
    hdus['indexmap_trunk_close_hdu'] = indexmap_trunk_close_hdu_all

    return props_all, hdus

def add_muse_info(props_all, muscat_table):

    print(f'[INFO] [add_muse_info] Adding MUSE catalouge info to final table...')
    keys = list(props_all.keys())
    # Add in MUSE informaiton
    for key in keys: 
        props_all[key] = join(props_all[key], muscat_table, 'muscat_id')
        props_all[key].sort('id')

    return(props_all)

def get_hdulists(hdus, regionIDs, outputdir='./'):

    keys = list(hdus.keys())
    for i in tqdm(range(len(hdus[keys[0]]))): 
        
        regionID = regionIDs[i]
        data = []
        header = []
        hdu = []

        for key in keys:
            data = np.array(hdus[key][i].data)
            header = hdus[key][i].header
            if key == keys[0]: 
                hdu += [fits.PrimaryHDU((), header=header)]
                hdu += [fits.ImageHDU(data, header, key)]
            else: 
                hdu += [fits.ImageHDU(data, header, key)]
        
        
        hdu_list = fits.HDUList(hdu)
        regionID_string = '%i' %regionID
        regionID_string = regionID_string.zfill(4)
        hdu_list.writeto('%s/hdus_%s.fits' %(outputdir, regionID_string), overwrite=True)


def remove_padding(data1, data2, buffer=5):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data1, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data1, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return (np.array([0]), np.array([0]))

    # Crop the data array
    cropped_data1 = data1[valid_y[0]-buffer:valid_y[-1]+1+buffer, valid_x[0]-buffer:valid_x[-1]+1+buffer]
    cropped_data2 = data2[valid_y[0]-buffer:valid_y[-1]+1+buffer, valid_x[0]-buffer:valid_x[-1]+1+buffer]
    
    # cropped_data1 = data1[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    # cropped_data2 = data2[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]

    return (cropped_data1, cropped_data2)


def get_flag_touch(props, hdu_mask):

    region_IDs = props['region_ID']
    flag_touch = np.ones(len(region_IDs))*0

    for i in tqdm(range(len(region_IDs))):

        region_ID = region_IDs[i]
        data_mask = hdu_mask.data

        region_mask = data_mask == region_ID
        region_invmask = (data_mask != region_ID) & (data_mask != -1)

        region_mask, region_invmask = remove_padding(region_mask, region_invmask)

        if len(region_mask) == 1: 
            region_mask_dilate = False
        else: 
            region_mask_dilate = binary_dilation(region_mask, iterations=2)

        if True in np.unique(region_mask_dilate&region_invmask): 
            flag_touch[i] = 1
        else: 
            flag_touch[i] = 0 

    props.add_column(Column(flag_touch, name='flag_nearby_hst')) 

    return(props)

# def get_flag_touch(props, hdu_mask):

#     region_IDs = props['region_ID']
#     flag_touch = np.ones(len(region_IDs))*0

#     for i in tqdm(range(len(region_IDs))):

#         region_ID = region_IDs[i]
#         data_mask = hdu_mask.data

#         region_mask = data_mask == region_ID
#         region_invmask = (data_mask != region_ID) & (data_mask != -1)

#         region_mask_dilate = binary_dilation(region_mask, iterations=2)

#         if True in np.unique(region_mask_dilate&region_invmask): 
#             flag_touch[i] = 1
#         else: 
#             flag_touch[i] = 0 

#     props.add_column(Column(flag_touch, name='flag_touch_hst')) 

#     return(props)
