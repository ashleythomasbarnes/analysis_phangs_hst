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


def get_padding(hdu_catalog_mask_, hdu_cluster_mask_):

    # Get the shapes of the data arrays
    catalog_shape = hdu_catalog_mask_.data.shape
    cluster_shape = hdu_cluster_mask_.data.shape

    # Calculate the difference in shape
    rows_diff = cluster_shape[0] - catalog_shape[0]
    cols_diff = cluster_shape[1] - catalog_shape[1]

    # Create arrays of -1 to append
    rows_to_append = -1 * np.ones((rows_diff, catalog_shape[1]))
    cols_to_append = -1 * np.ones((cluster_shape[0], cols_diff))

    # Append rows and columns to match the shape
    hdu_catalog_mask_.data = np.vstack((hdu_catalog_mask_.data, rows_to_append))
    hdu_catalog_mask_.data = np.hstack((hdu_catalog_mask_.data, cols_to_append))

    print('[INFO] Padding added to catalog mask to match cluster mask shape... ')

    return(hdu_catalog_mask_)


def get_cluster_IDs(region_ID, i, hdus_catalog_mask_new, hdus_cluster_mask_new):

    hdu_catalog_mask_ = hdus_catalog_mask_new[i] # Get the mask of the region Ha
    hdu_cluster_mask_ = hdus_cluster_mask_new[i] # Get the mask of the region cluster 

    # Check if the shapes of the masks are the same - if not, add padding to the catalog mask
    # Was case for NGC2835s due to half map missing
    if hdu_catalog_mask_.shape != hdu_cluster_mask_.shape:
        hdu_catalog_mask_ = get_padding(hdu_catalog_mask_, hdu_cluster_mask_)

    mask = hdu_catalog_mask_.data == region_ID

    data_cluster_masked = hdu_cluster_mask_.data[mask]

    cluster_IDs = np.unique(data_cluster_masked)
    cluster_IDs = cluster_IDs[cluster_IDs != 0]

    return cluster_IDs


def get_no_clusters(region_ID, props_clusters):
    
    props_clusters_new = QTable(props_clusters[0], masked=True) # Make new table with same columns as props_clusters 
    for colname in props_clusters_new.colnames: 
        props_clusters_new[colname] = np.nan # Fill with nans
    
    props_clusters_new['reg_dolflux_Age_MinChiSq_ave'] = np.nan
    props_clusters_new['reg_dolflux_Mass_MinChiSq_sum'] = np.nan

    props_clusters_new.mask = np.ones(len(props_clusters_new.columns), dtype=bool) # Mask all values

    props_clusters_new['no_clusters'] = True
    props_clusters_new['one_clusters'] = False
    props_clusters_new['multiple_clusters'] = False
    props_clusters_new['region_ID'] = region_ID

    return props_clusters_new


def get_one_clusters(region_ID, props_clusters, cluster_IDs):

    mask = props_clusters['reg_id'] == cluster_IDs[0]
    props_clusters_new = props_clusters[mask]

    props_clusters_new['reg_dolflux_Age_MinChiSq_ave'] = np.nan
    props_clusters_new['reg_dolflux_Mass_MinChiSq_sum'] = np.nan

    props_clusters_new['no_clusters'] = False
    props_clusters_new['one_clusters'] = True
    props_clusters_new['multiple_clusters'] = False
    props_clusters_new['region_ID'] = region_ID

    return props_clusters_new


def get_multi_clusters(region_ID, props_clusters, cluster_IDs):

    ages = []
    masses = []
    ages_ave = []
    masses_sum = []

    for cluster_ID in cluster_IDs: 
        mask = props_clusters['reg_id'] == cluster_ID
        ages += [props_clusters['reg_dolflux_Age_MinChiSq'][mask][0]]
        masses += [props_clusters['reg_dolflux_Mass_MinChiSq'][mask][0]]

    argmin = np.argmin(ages) # Get minimum age cluster
    ages_ave = np.mean(ages) # Get mean mass of clusters
    masses_sum = np.sum(masses) # Get sum of masses of clusters
    ages_massweighted = np.sum(np.array(masses) * np.array(ages)) / np.sum(masses) # Calculate the mass-weighted age

    mask = props_clusters['reg_id'] == cluster_IDs[argmin]
    props_clusters_new = props_clusters[mask]

    props_clusters_new['reg_dolflux_Age_MinChiSq_ave'] = ages_ave
    props_clusters_new['reg_dolflux_Mass_MinChiSq_sum'] = masses_sum
    props_clusters_new['reg_dolflux_Age_MinChiSq_massweighted'] = ages_massweighted

    props_clusters_new['no_clusters'] = False
    props_clusters_new['one_clusters'] = False
    props_clusters_new['multiple_clusters'] = True
    props_clusters_new['region_ID'] = region_ID

    return props_clusters_new