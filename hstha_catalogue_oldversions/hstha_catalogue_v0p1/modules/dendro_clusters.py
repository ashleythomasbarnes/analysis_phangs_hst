from astropy.table import Table, vstack, join
import numpy as np
import astropy.units as au
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm 
import os

import warnings
warnings.filterwarnings('ignore')

def get_clusters(cluster_table, hdus, props_all):

    # Get list... 
    cluster_table_masked_list = []

    # Create dummy table in case nothing is found...  
    cluster_table_dummy = Table(cluster_table, masked=True, copy=True) 
    for colname in cluster_table_dummy.colnames: 
        cluster_table_dummy[colname].mask = True
    cluster_table_dummy = cluster_table_dummy[0]

    for i in tqdm(range(len(props_all))):

        props_all_masked = props_all[i]

        indexmap_trunk_close_hdu = hdus['indexmap_trunk_close_hdu'][i].copy()
        indexmap_trunk_close_data = indexmap_trunk_close_hdu.data
        indexmap_trunk_close_mask = indexmap_trunk_close_data > -1

        ra_cl = cluster_table['PHANGS_RA'].quantity *au.deg
        dec_cl = cluster_table['PHANGS_DEC'].quantity *au.deg

        w = WCS(indexmap_trunk_close_hdu.header)
        c = SkyCoord(ra_cl,dec_cl,frame='icrs')
        x, y = w.world_to_pixel(c)

        lenx = indexmap_trunk_close_hdu.data.shape[0]
        leny = indexmap_trunk_close_hdu.data.shape[1]

        xmask = ~((x<5)|(x>lenx-5)|np.isnan(x)) 
        ymask = ~((y<5)|(y>leny-5)|np.isnan(y))
        xymask = xmask&ymask

        x = x[xymask]
        y = y[xymask]
        cluster_table_masked = cluster_table[xymask]

        x_int = np.array(np.round(x), dtype=int)
        y_int = np.array(np.round(y), dtype=int)
        inmask = indexmap_trunk_close_mask[y_int, x_int]
        x = x[inmask]
        y = y[inmask]
        x_int = x_int[inmask]
        y_int = y_int[inmask]
        cluster_table_masked = cluster_table_masked[inmask]

        # If nothing is found, set all to nan... 
        if len(cluster_table_masked) == 0: 

            cluster_table_dummy_tmp = Table(cluster_table_dummy, masked=True, copy=True) 
            cluster_table_dummy_tmp['region_ID'] = props_all_masked['region_ID']
            cluster_table_masked_list += [cluster_table_dummy_tmp]

            continue 

        # Get closest cluster to max... 
        c1 = SkyCoord(cluster_table_masked['PHANGS_RA'] *au.deg, cluster_table_masked['PHANGS_DEC'] *au.deg, frame='icrs')
        c2 = SkyCoord(props_all_masked['ra_max'] *au.deg, props_all_masked['dec_max'] *au.deg, frame='icrs')
        sep = c1.separation(c2)
        closemask = np.nanargmin(sep.to('arcsec'))

        x = x[closemask]
        y = y[closemask]
        x_int = x_int[closemask]
        y_int = y_int[closemask]
        cluster_table_masked = cluster_table_masked[closemask]

        cluster_table_masked = Table(cluster_table_masked)

        cluster_table_masked['x_cluster'] = x *au.pix
        cluster_table_masked['y_cluster'] = y *au.pix
        cluster_table_masked['region_ID'] = props_all_masked['region_ID']

        cluster_table_masked_list += [cluster_table_masked]

        cluster_table_masked_all = vstack(cluster_table_masked_list)
        props_all_cluster = join(props_all, cluster_table_masked_all, keys='region_ID')

    return(props_all_cluster)

def get_clusters_regions(props_all_cluster, dendro_dir, galaxy):

    ra, dec = props_all_cluster['PHANGS_RA'].data, props_all_cluster['PHANGS_DEC'].data

    os.system('rm %s/%s_clusters_regions.reg'  %(dendro_dir, galaxy))
    with open('%s/%s_clusters_regions.reg'  %(dendro_dir, galaxy), "w") as file:
        file.write("# Region file format: DS9 version 4.1\n")
        file.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        file.write("fk5\n")
        for i in tqdm(range(len(ra))):
            if ra[i] is not np.ma.masked:
                file.write("circle(%f,%f,0.15\")\n" %(ra[i],dec[i]))