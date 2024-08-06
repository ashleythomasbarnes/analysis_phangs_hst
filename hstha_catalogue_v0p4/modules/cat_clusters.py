import sys
sys.path.append('../modules/')
from cat_imports import *

def load_clusters(cluster_table_file, cluster_mask_file):

    print('[INFO] Loading %s' %cluster_table_file)
    print('[INFO] Loading %s' %cluster_mask_file)

    cluster_table = Table.read(cluster_table_file)

    def addunits(cluster_table):
        cluster_table['reg_dolflux_Age_MinChiSq'] = cluster_table['reg_dolflux_Age_MinChiSq'].quantity*au.Myr
        cluster_table['reg_dolflux_Mass_MinChiSq'] = cluster_table['reg_dolflux_Mass_MinChiSq'].quantity*au.Msun
        return(cluster_table)

    cluster_table = addunits(cluster_table)
    cluster_table = Table(cluster_table, masked=True, copy=False)

    cluster_mask = fits.open(cluster_mask_file)[0]
    
    return(cluster_table, cluster_mask)

def get_clusters_cutouts(cluster_mask, hdus, props_all):

    cluster_maskre_all = []

    # Loop through each region and extract cluster data
    for i in tqdm(range(len(props_all))):
    # for i in tqdm(range(10)):

        indexmap_trunk_close_hdu = hdus['indexmap_trunk_close_hdu'][i].copy()

        # Regrid the data from each associated catalog to match the indexmap trunk.
        cluster_maskre = reproject_interp(cluster_mask, indexmap_trunk_close_hdu.header, return_footprint=False, order='nearest-neighbor')
        cluster_maskre_all += [cluster_maskre]

    hdus['cluster_masks_hdu'] = cluster_maskre_all

    return hdus

def get_clusters(cluster_id, cluster_table):

    # Filter the table to get only rows that match the given IDs. 
    # `np.searchsorted` helps in efficiently finding indices where elements should be inserted to maintain order.
    cluster_table_masked = cluster_table[np.searchsorted(cluster_table['reg_id'], cluster_id)]

    # From the filtered table, get clusters with the minimum age value.
    min_age = np.min(cluster_table_masked['reg_dolflux_Age_MinChiSq'])
    cluster_table_masked = cluster_table_masked[cluster_table_masked['reg_dolflux_Age_MinChiSq'] == min_age]

    # If there are still multiple entries after the previous step, 
    # further filter the table to get only the cluster(s) with the maximum mass value.
    max_mass = np.max(cluster_table_masked['reg_dolflux_Mass_MinChiSq']) 
    cluster_table_masked = cluster_table_masked[cluster_table_masked['reg_dolflux_Mass_MinChiSq'] == max_mass]

    # If there is still more than one entry after all filtering steps, raise an alert.
    if len(cluster_table_masked) > 1:
        print('[WARNING] [get_clusters] STOP - please check, too many clusters!')
    
    cluster_table_masked.rename_column('reg_id', 'cluster_ID')

    return cluster_table_masked


def get_clusters_all(cluster_table, hdus, props_all):

    # Get list... 
    cluster_table_masked_list = []

    # Create dummy table in case nothing is found...  
    cluster_table_dummy = Table(cluster_table, masked=True, copy=True) 
    for colname in cluster_table_dummy.colnames: 
        cluster_table_dummy[colname].mask = True
    cluster_table_dummy = cluster_table_dummy[0]

    # Loop through each region and extract cluster data
    for i in tqdm(range(len(props_all))):
    # for i in tqdm(range(10)):

        props_all_masked = props_all[i]

        indexmap_trunk_close_hdu = hdus['indexmap_trunk_close_hdu'][i].copy()
        indexmap_trunk_close_data = indexmap_trunk_close_hdu.data
        indexmap_trunk_close_mask = indexmap_trunk_close_data > -1

        # Regrid the data from each associated catalog to match the indexmap trunk.
        # cluster_maskre = reproject_interp(cluster_mask, indexmap_trunk_close_hdu.header, return_footprint=False, order='nearest-neighbor')
        cluster_maskre = hdus['cluster_masks_hdu'][i].copy()

        # Mask out regions outside the trunk.
        cluster_maskre[~indexmap_trunk_close_mask] = np.nan
        cluster_maskre[cluster_maskre==0] = np.nan

        # Extract unique IDs for each associated catalog after masking.
        cluster_id = np.unique(cluster_maskre)[~np.isnan(np.unique(cluster_maskre))]

        # Depending on the available IDs, extract cluster information from the relevant catalog.
        if len(cluster_id) > 0:
            cluster_table_masked = get_clusters(cluster_id, cluster_table)
        
        # If nothing is found, set all to nan... 
        if len(cluster_id) == 0: 
            cluster_table_dummy_tmp = Table(cluster_table_dummy, masked=True, copy=True) 
            cluster_table_dummy_tmp['region_ID'] = props_all_masked['region_ID']
            cluster_table_masked_list += [cluster_table_dummy_tmp]
            continue 

        # Add the current region ID to the extracted table.
        cluster_table_masked['region_ID'] = props_all_masked['region_ID']

        cluster_table_masked_list += [cluster_table_masked]

    cluster_table_masked_all = vstack(cluster_table_masked_list)
    props_all_cluster = join(props_all, cluster_table_masked_all, keys='region_ID')
    
    return props_all_cluster