import sys
sys.path.append('../modules/')
from cat_imports import *

def load_associations(association_table_file, association_mask_file):

    print('[INFO] Loading %s' %association_table_file)
    print('[INFO] Loading %s' %association_mask_file)

    association_table = Table.read(association_table_file)

    def addunits(association_table):
        association_table['reg_dolflux_Age_MinChiSq'] = association_table['reg_dolflux_Age_MinChiSq'].quantity*au.Myr
        association_table['reg_dolflux_Mass_MinChiSq'] = association_table['reg_dolflux_Mass_MinChiSq'].quantity*au.Msun
        return(association_table)

    association_table = addunits(association_table)
    association_table = Table(association_table, masked=True, copy=False)

    association_mask = fits.open(association_mask_file)[0]
    
    return(association_table, association_mask)

def get_associations_cutouts(association_mask, hdus, props_all):

    association_maskre_all = []

    # Loop through each region and extract association data
    for i in tqdm(range(len(props_all))):
    # for i in tqdm(range(10)):

        indexmap_trunk_close_hdu = hdus['indexmap_trunk_close_hdu'][i].copy()

        # Regrid the data from each associated catalog to match the indexmap trunk.
        association_maskre = reproject_interp(association_mask, indexmap_trunk_close_hdu.header, return_footprint=False, order='nearest-neighbor')
        association_maskre_all += [association_maskre]

    hdus['association_masks_hdu'] = association_maskre_all

    return hdus

def get_padding(hdu_catalog_mask_, hdu_association_mask_):

    # Get the shapes of the data arrays
    catalog_shape = hdu_catalog_mask_.data.shape
    association_shape = hdu_association_mask_.data.shape

    # Calculate the difference in shape
    rows_diff = association_shape[0] - catalog_shape[0]
    cols_diff = association_shape[1] - catalog_shape[1]

    # Create arrays of -1 to append
    rows_to_append = -1 * np.ones((rows_diff, catalog_shape[1]))
    cols_to_append = -1 * np.ones((association_shape[0], cols_diff))

    # Append rows and columns to match the shape
    hdu_catalog_mask_.data = np.vstack((hdu_catalog_mask_.data, rows_to_append))
    hdu_catalog_mask_.data = np.hstack((hdu_catalog_mask_.data, cols_to_append))

    print('[INFO] Padding added to catalog mask to match association mask shape... ')

    return(hdu_catalog_mask_)


def get_association_IDs(region_ID, i, hdus_catalog_mask_new, hdus_association_mask_new):

    hdu_catalog_mask_ = hdus_catalog_mask_new[i] # Get the mask of the region Ha
    hdu_association_mask_ = hdus_association_mask_new[i] # Get the mask of the region association 

    # Check if the shapes of the masks are the same - if not, add padding to the catalog mask
    # Was case for NGC2835s due to half map missing
    if hdu_catalog_mask_.shape != hdu_association_mask_.shape:
        hdu_catalog_mask_ = get_padding(hdu_catalog_mask_, hdu_association_mask_)

    mask = hdu_catalog_mask_.data == region_ID

    data_association_masked = hdu_association_mask_.data[mask]

    association_IDs = np.unique(data_association_masked)
    association_IDs = association_IDs[association_IDs != 0]

    return association_IDs


def get_no_associations(region_ID, props_associations):
    
    props_associations_new = QTable(props_associations[0], masked=True) # Make new table with same columns as props_associations 
    for colname in props_associations_new.colnames: 
        props_associations_new[colname] = np.nan # Fill with nans
    
    props_associations_new['reg_dolflux_Age_MinChiSq_ave'] = np.nan
    props_associations_new['reg_dolflux_Mass_MinChiSq_sum'] = np.nan

    props_associations_new.mask = np.ones(len(props_associations_new.columns), dtype=bool) # Mask all values

    props_associations_new['no_associations'] = True
    props_associations_new['one_associations'] = False
    props_associations_new['multiple_associations'] = False
    props_associations_new['region_ID'] = region_ID

    return props_associations_new


def get_one_associations(region_ID, props_associations, association_IDs):

    mask = props_associations['reg_id'] == association_IDs[0]
    props_associations_new = props_associations[mask]

    props_associations_new['reg_dolflux_Age_MinChiSq_ave'] = np.nan
    props_associations_new['reg_dolflux_Mass_MinChiSq_sum'] = np.nan

    props_associations_new['no_associations'] = False
    props_associations_new['one_associations'] = True
    props_associations_new['multiple_associations'] = False
    props_associations_new['region_ID'] = region_ID

    return props_associations_new


def get_multi_associations(region_ID, props_associations, association_IDs):

    ages = []
    masses = []
    ages_ave = []
    masses_sum = []

    for association_ID in association_IDs: 
        mask = props_associations['reg_id'] == association_ID
        ages += [props_associations['reg_dolflux_Age_MinChiSq'][mask][0]]
        masses += [props_associations['reg_dolflux_Mass_MinChiSq'][mask][0]]

    argmin = np.argmin(ages) # Get minimum age association
    ages_ave = np.mean(ages) # Get mean mass of associations
    masses_sum = np.sum(masses) # Get sum of masses of associations
    ages_massweighted = np.sum(np.array(masses) * np.array(ages)) / np.sum(masses) # Calculate the mass-weighted age

    mask = props_associations['reg_id'] == association_IDs[argmin]
    props_associations_new = props_associations[mask]

    props_associations_new['reg_dolflux_Age_MinChiSq_ave'] = ages_ave
    props_associations_new['reg_dolflux_Mass_MinChiSq_sum'] = masses_sum
    props_associations_new['reg_dolflux_Age_MinChiSq_massweighted'] = ages_massweighted

    props_associations_new['no_associations'] = False
    props_associations_new['one_associations'] = False
    props_associations_new['multiple_associations'] = True
    props_associations_new['region_ID'] = region_ID

    return props_associations_new

def get_all_associations(region_ID, association_IDs, props_associations, association_name):

    # If no association, fill with nans...
    if len(association_IDs) == 0: 
        props_associations_new = get_no_associations(region_ID, props_associations)

    # Only one association needs to be selected...
    if len(association_IDs) == 1:
        props_associations_new = get_one_associations(region_ID, props_associations, association_IDs)
        
    # If more than one association, select the one with the minimum age...
    elif len(association_IDs) > 1:  
        props_associations_new = get_multi_associations(region_ID, props_associations, association_IDs)


    colnames_old = ['reg_id', 
                'reg_ra', 
                'reg_dec', 
                'reg_area', 
                'reg_rad', 
                'reg_dolflux_Age_MinChiSq', 
                'reg_dolflux_Age_MinChiSq_err', 
                'reg_dolflux_Mass_MinChiSq', 
                'reg_dolflux_Mass_MinChiSq_err', 
                'reg_dolflux_Ebv_MinChiSq', 
                'reg_dolflux_Ebv_MinChiSq_err', 
                'reg_dolflux_Age_MinChiSq_ave', 
                'reg_dolflux_Mass_MinChiSq_sum', 
                'no_associations', 
                'one_associations', 
                'multiple_associations', 
                'region_ID']

    colnames_new = ['region_ID_association_%s' %association_name, 
                    'ra_association_%s' %association_name, 
                    'dec_association_%s' %association_name, 
                    'area_association_%s' %association_name, 
                    'rad_association_%s' %association_name, 
                    'age_association_%s' %association_name, 
                    'age_err_association_%s' %association_name, 
                    'mass_association_%s' %association_name, 
                    'mass_err_association_%s' %association_name, 
                    'ebv_association_%s' %association_name, 
                    'ebv_err_association_%s' %association_name, 
                    'age_ave_association_%s' %association_name, 
                    'mass_sum_association_%s' %association_name, 
                    'no_associations_association_%s' %association_name, 
                    'one_associations_association_%s' %association_name, 
                    'multiple_associations_association_%s' %association_name,
                    'region_ID']

    props_associations_new = props_associations_new[colnames_old]
    props_associations_new.rename_columns(colnames_old, colnames_new)

    return props_associations_new