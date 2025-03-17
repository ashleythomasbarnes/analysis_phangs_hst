import sys
sys.path.append('../modules/')
from cat_imports import *

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


def load_table(filename):
    """
    Get galaxy properties from sample table
    """
    
    print(f'[INFO] [load_table] %s' %filename)
    sampletable = QTable.read(filename)
    return(sampletable)


def get_galaxyprops(galaxy, sampletable_file):
    """
    Get galaxy properties from sample table
    """
    
    print(f'[INFO] [get_galaxyprops] Getting sample table properties for {galaxy}...')
    sampletable = QTable.read(sampletable_file)
    mask = sampletable['name'] == galaxy
    sampletable = sampletable[mask]
    return(sampletable)

def get_museprops(galaxy, muscat_file):
    """
    Get properties from MUSE catalouge
    """

    print(f'[INFO] [get_MuseProps] Getting MUSE catalouge properties for {galaxy}...')
    muscat_table = QTable.read(muscat_file)
    muscat_table = muscat_table[muscat_table['gal_name'] == galaxy.swapcase()]
    return(muscat_table)

def checkmakedir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return

def get_unpack(file):
    if os.path.isfile(file):
        print('File found: %s' %file)
    elif os.path.isfile(file+'.gz'): 
        file_untar = file+'.gz'
        print('Unpacking: %s' %file_untar)
        os.system('gunzip --keep %s' %file_untar)
    else: 
        print('File not found: %s' %file)

def clean_unpack(file): 
    if os.path.isfile(file):
        if file.endswith('-idmask.fits'):
            print('Cleaning: %s' %file)
            os.system('rm %s' %file)
    else: 
        print('File not found: %s' %file)