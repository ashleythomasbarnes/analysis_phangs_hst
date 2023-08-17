import pickle
import numpy as np

def save_pickle(a, filename):
    """
    Save an object to a pickle file.

    Args:
        a: Object to be saved.
        filename (str): Path to the pickle file.

    Returns:
        None

    """
    with open(filename, 'wb') as handle:
        pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('[INFO] [save_pickle] Saved to %s' %filename)

def load_pickle(filename):
    """
    Load an object from a pickle file.

    Args:
        filename (str): Path to the pickle file.

    Returns:
        object: Loaded object.

    """
    with open(filename, 'rb') as handle:
        b = pickle.load(handle)
    print('[INFO] [load_pickle] Load %s' %filename)
    return b

def find_nearest(array, value):
    """find nearest value in array, and return id"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx
