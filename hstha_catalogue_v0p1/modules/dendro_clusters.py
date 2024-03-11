import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table, hstack, vstack, join
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as ac
import astropy.units as au
from glob import glob
from spectral_cube import SpectralCube
import scipy 
from reproject import reproject_interp
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm 
from astropy.io import fits
import matplotlib as mpl
import pyregion
import aplpy
import math
import os
import pickle

def get_asscats(dir_asscat='../../../data/cluster_catalogs/ngc628-vselect/'):

    # TO DO: Make more general
    asscat1_file = '%s/ws16pc/PHANGS_IR3_hst_wfc3_ngc628_v1p3_multi_assoc-vselect-ws16pc-main.fits' %dir_asscat #HST association catalogue
    asscat2_file = '%s/ws32pc/PHANGS_IR3_hst_wfc3_ngc628_v1p3_multi_assoc-vselect-ws32pc-main.fits' %dir_asscat #HST association catalogue
    asscat3_file = '%s/ws64pc/PHANGS_IR3_hst_wfc3_ngc628_v1p3_multi_assoc-vselect-ws64pc-main.fits' %dir_asscat #HST association catalogue

    print(f'[INFO] [get_asscats] Cluster catalouge file: {asscat1_file}')
    print(f'[INFO] [get_asscats] Cluster catalouge file: {asscat2_file}')
    print(f'[INFO] [get_asscats] Cluster catalouge file: {asscat3_file}')

    asscat1_table = Table.read(asscat1_file)
    asscat2_table = Table.read(asscat2_file)
    asscat3_table = Table.read(asscat3_file)

    def addunits(asscat_table):
        asscat_table['reg_dolflux_Age_MinChiSq'] = asscat_table['reg_dolflux_Age_MinChiSq'].quantity*au.Myr
        asscat_table['reg_dolflux_Mass_MinChiSq'] = asscat_table['reg_dolflux_Mass_MinChiSq'].quantity*au.Msun
        return(asscat_table)

    asscat1_table = addunits(asscat1_table)
    asscat2_table = addunits(asscat2_table)
    asscat3_table = addunits(asscat3_table)
    
    asscat1_table = Table(asscat1_table, masked=True, copy=False)
    asscat2_table = Table(asscat2_table, masked=True, copy=False)
    asscat3_table = Table(asscat3_table, masked=True, copy=False)
    
    return((asscat1_file, asscat2_file, asscat3_file), (asscat1_table, asscat2_table, asscat3_table))


def get_clusters(asscat_id, asscat_table, asscat_file):
    """
    Extracts cluster information based on given IDs and associated catalog table.
    
    Parameters:
    - asscat_id: Array of IDs for which cluster data needs to be fetched.
    - asscat_table: The table containing cluster data for various IDs.
    - asscat_file: The filename of the associated catalog.
    
    Returns:
    - asscat_table: A table with extracted cluster data based on the given IDs and certain conditions.
    """

    # Append the filename to the table so that we can trace back to the source file.
    asscat_table['asscat'] = asscat_file.split('/')[-1]  # Extract only the file name from the path

    # Filter the table to get only rows that match the given IDs. 
    # `np.searchsorted` helps in efficiently finding indices where elements should be inserted to maintain order.
    asscat_table = asscat_table[np.searchsorted(asscat_table['reg_id'], asscat_id)]

    # From the filtered table, get clusters with the minimum age value.
    min_age = np.min(asscat_table['reg_dolflux_Age_MinChiSq'])
    asscat_table = asscat_table[asscat_table['reg_dolflux_Age_MinChiSq'] == min_age]

    # If there are still multiple entries after the previous step, 
    # further filter the table to get only the cluster(s) with the maximum mass value.
    max_mass = np.max(asscat_table['reg_dolflux_Mass_MinChiSq']) 
    asscat_table = asscat_table[asscat_table['reg_dolflux_Mass_MinChiSq'] == max_mass]

    # If there is still more than one entry after all filtering steps, raise an alert.
    if len(asscat_table) > 1:
        print('[WARNING] [get_clusters] STOP - please check, too many clusters!')
        
    return asscat_table


def get_clusters_all(regions, asscat_files, asscat_tables, hdus):
    """
    Function to get all clusters for each region based on associated catalog files and tables.

    Parameters:
    - regions: The regions of interest.
    - asscat_files: List of associated catalog files for each region.
    - asscat_tables: List of tables corresponding to the associated catalog files.

    Returns:
    - clust_all_final: A table containing information on all the clusters across all regions.
    """

    # Decompose the associated catalog files and tables into individual variables.
    asscat1_file, asscat2_file, asscat3_file = asscat_files
    asscat1_table, asscat2_table, asscat3_table = asscat_tables

    asscat_hdus = (hdus['asscat1_hdus'], hdus['asscat2_hdus'], hdus['asscat3_hdus'])
    indexmap_trunk_hdus = hdus['indexmap_trunk_hdus_3sig']
    asscat1_hdus, asscat2_hdus, asscat3_hdus = asscat_hdus

    n = len(regions['ra'])  # Number of regions
    clust_nonefound = []  # List to store region IDs where no clusters are found

    # Try to delete the clust_all_final variable if it exists from previous runs.
    try: 
        del clust_all_final
    except: 
        pass
   
    # Loop through each region and extract cluster data
    for regionID in tqdm(range(n), desc='Get clusters', position=0): 

        # Extract data for each associated catalog for the current region.
        asscat1_hdu = asscat1_hdus[regionID]
        asscat2_hdu = asscat2_hdus[regionID]
        asscat3_hdu = asscat3_hdus[regionID]
        indexmap_trunk_hdu = indexmap_trunk_hdus[regionID]

        mask_trunk = indexmap_trunk_hdu.data != -1

        # Regrid the data from each associated catalog to match the indexmap trunk.
        asscat1_data = reproject_interp(asscat1_hdu, indexmap_trunk_hdu.header, return_footprint=False, order='nearest-neighbor')
        asscat2_data = reproject_interp(asscat2_hdu, indexmap_trunk_hdu.header, return_footprint=False, order='nearest-neighbor')
        asscat3_data = reproject_interp(asscat3_hdu, indexmap_trunk_hdu.header, return_footprint=False, order='nearest-neighbor')

        # Mask out regions outside the trunk.
        asscat1_data[~mask_trunk] = np.nan
        asscat2_data[~mask_trunk] = np.nan
        asscat3_data[~mask_trunk] = np.nan

        # Extract unique IDs for each associated catalog after masking.
        asscat1_id = np.unique(asscat1_data)[~np.isnan(np.unique(asscat1_data))]
        asscat2_id = np.unique(asscat2_data)[~np.isnan(np.unique(asscat2_data))]
        asscat3_id = np.unique(asscat3_data)[~np.isnan(np.unique(asscat3_data))]

        # Depending on the available IDs, extract cluster information from the relevant catalog.
        if len(asscat2_id) > 0:
            asscat_table = get_clusters(asscat2_id, asscat2_table, asscat2_file)
            asscat_table.rename_column('reg_id', 'asscat2_id')
        elif len(asscat1_id) > 0:
            asscat_table = get_clusters(asscat1_id, asscat1_table, asscat1_file)
            asscat_table.rename_column('reg_id', 'asscat1_id')
        elif len(asscat3_id) > 0:
            asscat_table = get_clusters(asscat3_id, asscat3_table, asscat3_file)
            asscat_table.rename_column('reg_id', 'asscat3_id')
        else: 
            clust_nonefound.append(regionID)
            print(f'[INFO] [regionID={regionID}] No clusters found.')
            continue

        # Add the current region ID to the extracted table.
        asscat_table['id'] = regionID

        # Try stacking the extracted cluster table to the final table.
        # If the final table does not exist, initialize it with the current extracted table.
        try:
            clust_all_final = vstack([clust_all_final, asscat_table])
        except: 
            print(f'[INFO] [regionID={regionID}] Initializing clust_all_final table...')
            clust_all_final = asscat_table

    # For regions where no clusters were found, mask all columns in the final table except the region ID.
    for regionID in tqdm(clust_nonefound, desc='Masking regions with no clusters', position=0): 
        clust_all_final.add_row(clust_all_final[0])
        clust_all_final[-1]['id'] = regionID
        for colname in clust_all_final[0].colnames: 
            if colname != 'id':
                clust_all_final[-1][colname] = np.ma.masked

    # Sort the final table based on region ID.
    clust_all_final.sort('id')
    
    return clust_all_final
    
def save_fits(props_all_final, props_all_file):
    filename = props_all_file.replace('01','02')
    print('[INFO] [save_fits] Saved to %s' %filename)
    props_all_final.write(filename, overwrite=True)