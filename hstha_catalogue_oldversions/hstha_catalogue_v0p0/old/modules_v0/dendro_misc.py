import pickle
import numpy as np
from astropy.table import Table, join

def save_pickle(a, filename):
    """
    Save an object to a pickle file.
    """
    with open(filename, 'wb') as handle:
        pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('[INFO] [save_pickle] Saved to %s' %filename)

def load_pickle(filename):
    """
    Load an object from a pickle file.
    """
    with open(filename, 'rb') as handle:
        b = pickle.load(handle)
    print('[INFO] [load_pickle] Load %s' %filename)
    return b

def find_nearest(array, value):
    """
    Find nearest value in array, and return id
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def convert_to_float32(hdu):
    """
    Conversion function to float32
    """
    hdu.data = hdu.data.astype('float32')
    return(hdu)


def get_galaxyprops(galaxy, sampletable_file):
    """
    Get galaxy properties from sample table
    """
    print(f'[INFO] [get_galaxyprops] Getting sample table properties for {galaxy}...')
    sampletable = Table.read(sampletable_file)
    mask = sampletable['name'] == galaxy
    sampletable = sampletable[mask]

    return(sampletable)

def get_museprops(galaxy, muscat_file):
    """
    Get properties from MUSE catalouge
    """
    muscat_table = Table.read(muscat_file)
    muscat_table = muscat_table[muscat_table['gal_name'] == galaxy.swapcase()]
    muscat_table.rename_column('region_ID', 'muscat_id')

    # Get MUSE REFITTED data for galaxy
    muscat_newfits_file = '/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/catalouge/v2/raw/Nebulae_catalogue_v2_refitNII_refitTe.fits' #MUSE catalogue with properties

    muscat_newfits_table = Table.read(muscat_newfits_file)
    muscat_newfits_table = muscat_newfits_table[muscat_newfits_table['GAL_NAME'] == galaxy.swapcase()]
    muscat_newfits_table.rename_column('REGION_ID', 'muscat_id')
    muscat_newfits_table = muscat_newfits_table['muscat_id', 'T_N2_REFIT']
    muscat_table = join(muscat_table, muscat_newfits_table, 'muscat_id')
    
    print(f'[INFO] [get_MuseProps] Getting MUSE catalouge properties for {galaxy}...')
    return(muscat_table)

# def get_museprops(galaxy, muscat_file):
#     """
#     Get properties from MUSE catalouge
#     """
#     muscat_table = Table.read(muscat_file)
#     muscat_table = muscat_table[muscat_table['gal_name'] == galaxy.swapcase()]
#     muscat_table.rename_column('region_ID', 'muscat_id')

#     # Get MUSE REFITTED data for galaxy
#     muscat_newfits_file = '/Users/abarnes/Dropbox/work/Projects/pressures/phangs/data/catalouge/v2/raw/Nebulae_catalogue_v2_refitNII_refitTe.fits' #MUSE catalogue with properties

#     muscat_newfits_table = Table.read(muscat_newfits_file)
#     muscat_newfits_table = muscat_newfits_table[muscat_newfits_table['GAL_NAME'] == galaxy.swapcase()]
#     muscat_newfits_table.rename_column('REGION_ID', 'muscat_id')
#     muscat_newfits_table = muscat_newfits_table['muscat_id', 'T_N2_REFIT']
#     muscat_table = join(muscat_table, muscat_newfits_table, 'muscat_id')
    
#     print(f'[INFO] [get_MuseProps] Getting MUSE catalouge properties for {galaxy}...')
#     return(muscat_table)
