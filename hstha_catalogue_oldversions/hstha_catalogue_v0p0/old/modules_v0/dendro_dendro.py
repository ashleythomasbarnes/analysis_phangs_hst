import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table, join
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

def get_maskeddata(muscat_hdu, hst_hdu): 
    
    hst_masked_hdu = hst_hdu.copy() #only ID in mask
    hst_maskedall_hdu = hst_hdu.copy() #all MUSE regions out of mask - for noise measure
    
    #get mask 
    shape = muscat_hdu.data.shape
    muscat_id = muscat_hdu.data[int(shape[0]/2),int(shape[0]/2)]
    
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
    hst_maskedall_hdu.data[np.isnan(data2)] = np.nan
    
    return(hst_masked_hdu, hst_maskedall_hdu, muscat_id)

def get_maskedhdus(hdus, regions, hstha_hdu_name='hstha_hdu'):
    # Initialize lists to store the processed data
    hstha_hdu_smooth = []
    hstha_masked_hdu = []
    hstha_maskedall_hdu = []
    muscat_id = []
    mask_muse = []

    def process_region(regionID):
        # Load the necessary data
        muscat_hdu = hdus['muscat_hdu'][regionID]
        hstha_hdu = hdus[hstha_hdu_name][regionID].copy()

        # Smooth the data
        kernel = Gaussian2DKernel(x_stddev=0.5)
        hstha_hdu.data = convolve(hstha_hdu.data, kernel)

        # Mask the data
        hstha_masked_hdu, hstha_maskedall_hdu, muscat_id = get_maskeddata(muscat_hdu, hstha_hdu)
        mask_muse = fits.PrimaryHDU(~np.isnan(hstha_masked_hdu.data) * 1, hstha_hdu.header)

        # Return the processed data along with the regionID
        return regionID, hstha_hdu, hstha_masked_hdu, hstha_maskedall_hdu, muscat_id, mask_muse

    print(f'[INFO] [get_maskedhdus] Getting HST maps masked by MUSE catalouge...')
    for regionID in tqdm(range(len(regions['ra'])), desc='Masking regions', position=0):

        # Extract the results in the correct order
        results = process_region(regionID)
        hstha_hdu_smooth += [results[1]]
        hstha_masked_hdu += [results[2]]
        hstha_maskedall_hdu += [results[3]]
        muscat_id += [results[4]]
        mask_muse += [results[5]]

    # Assign the processed data to the corresponding keys in the hdus dictionary
    hdus['%s_smooth' % hstha_hdu_name] = hstha_hdu_smooth
    hdus['%s_smooth_masked' % hstha_hdu_name] = hstha_masked_hdu
    hdus['%s_smooth_maskedall' % hstha_hdu_name] = hstha_maskedall_hdu
    hdus['musmask_hdu'] = mask_muse

    # Return the modified hdus dictionary and muscat_id
    return hdus, muscat_id

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
    # hdu = hst06_masked_hdu
    wcs = WCS(hdu.header)
    header = hdu.header
    data = hdu.data
    data_nomask = hdu_nomask.data
    data_outmask = hdu_outmask.data

    # Calculate statistics for dendrogram
    std = stats.mad_std(data_outmask, ignore_nan=True)  # Get noise
    std = stats.mad_std(data_outmask[data_outmask<20*std], ignore_nan=True)  # Get noise below threshold

    pixsize = np.array([np.abs(header['CDELT1']), np.abs(header['CDELT2'])]).mean() * au.degree
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
    props_trunk_ = pp_catalog(dendro.trunk, metadata, verbose=False)  # Get table

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

    ax = aplpy.FITSFigure(hdu)
    ra_struc = ax.pixel2world(props['x_cen'], props['y_cen'])[0][0]
    dec_struc = ax.pixel2world(props['x_cen'], props['y_cen'])[1][0]
    plt.close('all')

    props['ra_cen'] = ra_struc * au.deg
    props['dec_cen'] = dec_struc * au.deg
    props.rename_column('radius', 'mean_sigma')
    props['mean_sigma_pc'] = props['mean_sigma'].quantity[0].to('arcsec').value * pcperarcsec * au.pc

    flux = np.nansum(hdu.data[indexmap_trunk_close_hdu.data != -1])
    props['flux'] = flux
    props['flux'].unit = au.erg / au.s / au.cm**2

    props = Table(props, masked=True, copy=False)

    return dendro, props, indexmap_trunk_hdu, indexmap_trunk_close_hdu


def all_nan_check(hdus, regionID):
    
    data = hdus['hstha_hdu_smooth'][regionID].data
    shape = data.shape
    xmin = shape[0]/2 - shape[0]/4
    xmax = shape[0]/2 + shape[0]/4
    ymin = shape[1]/2 - shape[1]/4
    ymax = shape[1]/2 + shape[1]/4
    
    data = data[int(xmin):int(xmax),int(ymin):int(ymax)]
    
    if np.isnan(data).all():
        # print(f'[INFO] [props_all] [regionID={regionID}] All nan values in map...')
        return (True)
    else:
        return (False)


def get_dedro_all(hdus, regionID, min_value_sig, sampletable, muscat_ids):
    
    # Initialize variables for storing the results
    min_value_sig_n = len(min_value_sig)
    props = ['']*min_value_sig_n
    indexmap_trunk_hdu = ['']*min_value_sig_n
    indexmap_trunk_close_hdu = ['']*min_value_sig_n

    # Process each sigma value
    for i in range(min_value_sig_n):
        
        output = get_dedro(hdus['hstha_hdu_smooth_masked'][regionID],
                            hdus['hstha_hdu_smooth_maskedall'][regionID], 
                            hdus['hstha_hdu_smooth'][regionID],
                            min_npix=9, 
                            sampletable=sampletable,
                            min_value_sig=min_value_sig[i], 
                            min_delta_sig=3)

        dendro, props[i], indexmap_trunk_hdu[i], indexmap_trunk_close_hdu[i] = output
    
    # Add more info for each sigma value
    for i in range(min_value_sig_n):
        
        # Replace values
        if props[i] is None: 
            
            if all_nan_check(hdus, regionID): 

                print(f'[INFO] [get_dedro] [regionID={regionID}] [min_value_sig={min_value_sig[i]}] No dendro found, all nan values... ')
            else: 
                print(f'[INFO] [get_dedro] [regionID={regionID}] [min_value_sig={min_value_sig[i]}] No dendro found...')
            
            # Replace values with some dummy dendro table so we have all the correct columns, will be masked later anyway... 
            j=0
            while props[i] is None: 
                j+=1 

                print(hdus['hstha_hdu_smooth_masked'][j])
                print(hdus['hstha_hdu_smooth_masked'][j])
                print(hdus['hstha_hdu_smooth_maskedall'][j])
                print(hdus['hstha_hdu_smooth_maskedall'][j])
                

                props[i] = get_dedro(hdus['hstha_hdu_smooth_masked'][j],
                                    hdus['hstha_hdu_smooth_maskedall'][j], 
                                    hdus['hstha_hdu_smooth'][j],
                                    sampletable=sampletable,
                                    min_npix=0, 
                                    min_value_sig=-100, 
                                    min_delta_sig=-100)[1]

            for colname in props[i].colnames: 
                props[i][colname] = np.ma.masked
        
        props[i].add_column(Column(regionID, name='id'), index=0)
        props[i].add_column(Column(muscat_ids[regionID], name='muscat_id'), index=1)  

    # Return the results along with the regionID
    return regionID, props, indexmap_trunk_hdu

def get_dedro_all_decorator(hdus, min_value_sig, sampletable, muscat_ids):
    
    min_value_sig_n = len(min_value_sig)
    min_value_sig_str = ['%ssig' %min_value_sig[i] for i in range(min_value_sig_n)]
    props_all = dict.fromkeys(min_value_sig_str)
    n = len(hdus[list(hdus.keys())[0]])

    for key in min_value_sig_str: 
        hdus['indexmap_trunk_hdu_%s' %key] = []

    for regionID in tqdm(range(n), desc='Dendrogram (with sigma ranges)', position=0):

        # if regionID!=550:
        #     continue

        regionID, props, indexmap_trunk_hdu = get_dedro_all(hdus, regionID, min_value_sig, sampletable, muscat_ids)

        try: 
            for i in range(min_value_sig_n): 
                # print(f'[INFO] [props_all] [regionID={regionID}] [min_value_sig={min_value_sig[i]}] Appending props_all table...')
                props_all['%isig'%min_value_sig[i]].add_row(props[i][0])
        except: 
            for i in range(min_value_sig_n): 
                print(f'[INFO] [props_all] [regionID={regionID}] [min_value_sig={min_value_sig[i]}] Initialising props_all table...')
                props_all['%isig'%min_value_sig[i]] = props[i].copy()        

        for i in range(min_value_sig_n): 
            try: 
                hdus['indexmap_trunk_hdu_%ssig' %min_value_sig[i]] += [indexmap_trunk_hdu[i].copy()]
            except: 
                print(f'[INFO] [props_all] [regionID={regionID}] [min_value_sig={min_value_sig[i]}] No sigma mask, all blank added to indexmap_trunk_hdu...')
                shape = hdus['hst07_hdu_smooth_masked'][regionID].data.shape
                header = hdus['hst07_hdu_smooth_masked'][regionID].header
                hdu_empty = fits.PrimaryHDU(np.ones(shape)*-1, header)
                hdus['indexmap_trunk_hdu_%ssig' %min_value_sig[i]] += [hdu_empty.copy()]
    
    return(hdus, props_all)

def add_muse_info(props_all, muscat_table):

    print(f'[INFO] [add_muse_info] Adding MUSE catalouge info to final table...')
    keys = list(props_all.keys())
    # Add in MUSE informaiton
    for key in keys: 
        props_all[key] = join(props_all[key], muscat_table, 'muscat_id')
        props_all[key].sort('id')

    return(props_all)

def save_fits(props_all, output_dir):
    keys = list(props_all.keys())
    for key in keys: 
        filename = '%s/01_props_all_%s.fits' %(output_dir,key)
        print('[INFO] [save_fits] Saved to %s' %filename)
        props_all[key].write(filename, overwrite=True)
