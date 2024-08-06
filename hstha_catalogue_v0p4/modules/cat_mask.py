import sys
sys.path.append('../modules/')
from cat_imports import *

def get_threshmask(data, rms, offset=0, thresh=0):
    mask = data > (rms * thresh) + offset
    return mask


def get_circmask(radius=5, h=11, w=11):
    center = int(h/2)
    X, Y = np.ogrid[:h, :w]
    dist = np.sqrt((X - center)**2 + (Y-center)**2)
    mask = dist <= radius
    return mask *1


def get_prunemask(mask, thresh=10):

    mask_out = mask.copy()
    l, j = ndimage.label(mask_out)
    hist = ndimage.measurements.histogram(l, 0, j+1, j+1)
    os = ndimage.find_objects(l)

    for i in range(j):
        if hist[i+1]<thresh:
            mask_out[os[i]] = 0

    return(mask_out)


def get_regions(ra, dec, height, ellipse=False, width=0, angle=0):

    center_sky = SkyCoord(ra, dec, unit='deg', frame='fk5')
    if ellipse:
        region_sky = [EllipseSkyRegion(center, h, w, a) for center, h, w, a in zip(center_sky, height, width, angle)]
    else: 
        region_sky = [CircleSkyRegion(center, rad) for center, rad in zip(center_sky, height)]
    region_sky = Regions(region_sky)

    return region_sky


def get_ds9regions(props_all, outputfile='ds9'):

    region_sky = get_regions(props_all['ra_com'], props_all['dec_com'], props_all['radius_circ'])
    region_sky.write(outputfile+'_com.reg', format='ds9', overwrite=True)

    region_sky = get_regions(props_all['ra_com'], props_all['dec_com'], np.ones(len(props_all['dec_com']))*(0.1/3600.)*au.deg)
    region_sky.write(outputfile+'_compoint.reg', format='ds9', overwrite=True)

    region_sky = get_regions(props_all['ra_max'], props_all['dec_max'], props_all['radius_circ'])
    region_sky.write(outputfile+'_max.reg', format='ds9', overwrite=True)

    region_sky = get_regions(props_all['ra_max'], props_all['dec_max'], np.ones(len(props_all['dec_com']))*(0.1/3600.)*au.deg)
    region_sky.write(outputfile+'_maxpoint.reg', format='ds9', overwrite=True)

    region_sky = get_regions(props_all['ra_mom'], props_all['dec_mom'], props_all['major_sigma']*2, True, props_all['minor_sigma']*2, props_all['position_angle'])
    region_sky.write(outputfile+'_mom.reg', format='ds9', overwrite=True)


def get_hdumask(hdu_ha, hdu_mask, outputfile='tmp.fits'):
    
    # Set the shape of the reference H-alpha data and initialize its data array to -1
    shape_ha = hdu_ha.data.shape
    hdu_ha.data = np.ones(shape_ha) * -1
    wcs_ha = WCS(hdu_ha.header)

    # Loop over the HDUs provided in the 'indexmap_trunk_hdus_3sig' key
    for hdu_id in tqdm(range(len(hdu_mask)), desc='Masking regions'):

        # Count of non-negative values in the mask
        n_mask = np.sum(hdu_mask[hdu_id].data!=-1)  
        if n_mask == 0: 
            continue

        # Extract the current HDU mask and its WCS
        hdu_mask_ = hdu_mask[hdu_id]
        data_mask = (hdu_mask_.data != -1)  # Convert to a Boolean mask
        wcs_mask = WCS(hdu_mask_.header)

        sky = wcs_mask.pixel_to_world(0,0)
        # Convert the sky coordinate to pixel coordinate in the reference H-alpha HDU
        pix = wcs_ha.world_to_pixel(sky)
        pix = np.array(pix, dtype=np.float64)
        x = round(pix[1])
        y = round(pix[0])
        w, h = hdu_mask_.data.shape
        hdu_ha.data[x:x+w, y:y+h][data_mask] = hdu_mask_.data[data_mask]

    # hdu_ha.data = np.int16(hdu_ha.data)
    hdu_ha.data = hdu_ha.data.astype('int16')

    # Write the resulting mask to an output file
    hdu_ha.writeto(outputfile, overwrite=True)


def get_hducomplex(props_all_final, inputfile, outputfile):
    hdu = fits.open(inputfile)[0]
    for region_ID, complexity_score in tqdm(props_all_final['region_ID', 'complexity_score']):
        hdu.data[hdu.data==region_ID] = complexity_score
    hdu.writeto(outputfile, overwrite=True)


def get_maskeddata(muscat_hdu, hst_hdu, muscat_id): 
    
    hst_masked_hdu = hst_hdu.copy() #only ID in mask
    hst_masked_ones_hdu = hst_hdu.copy() #only ID in mask
    hst_maskedall_hdu = hst_hdu.copy() #all MUSE regions out of mask - for noise measure
    
    #get mask   
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

    hstha_sm_arr = []

    hstha_m_arr = []
    hstha_m_sm_arr = []
    hstha_m_ones_arr = []
    hstha_m_ones_sm_arr = []
    hstha_mall_arr = []
    hstha_mall_sm_arr = []
    muscat_id_arr = []
    mask_muse_arr = []

    print(f'[INFO] [get_maskedhdus] Getting HST maps masked by MUSE catalouge...')
    for i in tqdm(range(len(regions['ra'])), desc='Masking regions', position=0):

        # Load
        muscat_hdu = hdus['muscat_hdu'][i]
        hstha_hdu = hdus[hstha_hdu_name][i].copy()
        hstha_hdu_smooth = hdus[hstha_hdu_name][i].copy()
        muscat_regionID = muscat_regionIDs[i]

        # Smooth the data
        kernel = Gaussian2DKernel(x_stddev=0.5)
        hstha_hdu_smooth.data = convolve(hstha_hdu_smooth.data, kernel)

        # Mask the data
        output1 = get_maskeddata(muscat_hdu, hstha_hdu, muscat_regionID)
        hstha_m_hdu, hst_m_ones_hdu, hstha_mdall_hdu, _ = output1

        output2 = get_maskeddata(muscat_hdu, hstha_hdu_smooth, muscat_regionID)
        hstha_m_sm_hdu, hst_m_ones_sm_hdu, hstha_mall_sm_hdu, _ = output2

        # Make HDU
        mask_muse = fits.PrimaryHDU(~np.isnan(hstha_m_hdu.data) * 1, hstha_hdu_smooth.header)

        # Extract the results in the correct order
        hstha_sm_arr += [hstha_hdu_smooth]
        hstha_m_arr += [hstha_m_hdu]
        hstha_m_sm_arr += [hstha_m_sm_hdu]
        hstha_m_ones_arr += [hst_m_ones_hdu]
        hstha_m_ones_sm_arr += [hst_m_ones_sm_hdu]
        hstha_mall_arr += [hstha_mdall_hdu]
        hstha_mall_sm_arr += [hstha_mall_sm_hdu]

        muscat_id_arr += [muscat_regionID]
        mask_muse_arr += [mask_muse]

    # Assign the processed data to the corresponding keys in the hdus dictionary
    hdus['%s_smooth' % hstha_hdu_name] = hstha_sm_arr
    hdus['%s_masked' % hstha_hdu_name] = hstha_m_arr
    hdus['%s_smooth_masked' % hstha_hdu_name] = hstha_m_sm_arr
    hdus['%s_masked_ones' % hstha_hdu_name] = hstha_m_ones_arr
    hdus['%s_smooth_masked_ones' % hstha_hdu_name] = hstha_m_ones_sm_arr
    hdus['%s_maskedall' % hstha_hdu_name] = hstha_mall_arr
    hdus['%s_smooth_maskedall' % hstha_hdu_name] = hstha_mall_sm_arr
    hdus['musmask_hdu'] = mask_muse_arr

    return hdus