import warnings
warnings.filterwarnings('ignore')

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
from tqdm import tqdm
from astropy.io import fits
import matplotlib as mpl
import pyregion
import aplpy
import math
import os
import pickle
from regions import Regions, RectangleSkyRegion, CircleSkyRegion, EllipseSkyRegion
from astropy.nddata import Cutout2D
import gc
from scipy import ndimage
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.table import QTable
import cv2 # Import OpenCV


import os
import sys
import warnings
import numpy as np
from glob import glob 
from scipy import ndimage
# from tqdm.auto import tqdm
from tqdm import tqdm
from astropy import stats
import astropy.units as au
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import QTable, join, vstack, Column
from astrodendro import Dendrogram, pp_catalog
from reproject import reproject_interp

sys.path.append('./../')
from modules import cat_cutouts, cat_misc, cat_props, cat_mask, cat_associations

warnings.filterwarnings('ignore')

# import warnings
# warnings.filterwarnings('ignore')

# from astropy.table import Table, hstack, vstack, join
# import numpy as np
# import matplotlib.pyplot as plt
# import astropy.constants as ac
# import astropy.units as au
# from glob import glob
# from spectral_cube import SpectralCube
# import scipy 
# from reproject import reproject_interp
# from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
# from tqdm.auto import tqdm 
# from astropy.io import fits
# import matplotlib as mpl
# import pyregion
# import aplpy
# import math
# import os
# import pickle

# import warnings
# warnings.filterwarnings('ignore')

# from regions import Regions, RectangleSkyRegion
# import numpy as np
# import astropy.units as au
# from astropy.io import fits
# from tqdm import tqdm 
# from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
# from astropy.nddata import Cutout2D
# import gc

# from astropy.wcs import WCS
# from astropy.io import fits
# from tqdm import tqdm
# from astropy.coordinates import SkyCoord
# from regions import CircleSkyRegion, EllipseSkyRegion, Regions
# import numpy as np
# from scipy import ndimage
# from reproject import reproject_interp
# import astropy.units as au
# from astropy.convolution import Gaussian2DKernel
# from astropy.convolution import convolve


# import pickle
# import numpy as np
# from astropy.table import Table
# import os 

# import numpy as np
# import astropy.constants as ac
# import astropy.units as au
# import matplotlib.pyplot as plt
# from tqdm import tqdm
