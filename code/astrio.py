#import astropy
#from astropy import coordinates
from astropy import wcs
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import hstack, vstack, Table, Column
from astropy.io import ascii
from astropy.visualization import (MinMaxInterval, ZScaleInterval, PercentileInterval,
                                   SqrtStretch, LinearStretch, LogStretch, HistEqStretch,
                                   ImageNormalize)
#from astropy.units import Quantity
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus
from astroquery.skyview import SkyView

# from astropy.visualization import (MinMaxInterval, ZScaleInterval, PercentileInterval, )
# jupyter matplotlib backends
# inconsistent? windows popping up?
# calling "notebook" reverts to "nbagg"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize
import pandas as pd
from tkinter import *

import hdbscan

from getCI import *

from datetime import datetime

# TODO: needs documentation update
# TODO: needs implementation
def file_to_table(self, fname, tname=None):

    pass
    # code to pull a file into an astropy table
    #returns the table if tname=None
    #appends the table to tname if != None

# TODO: needs documentation update
# TODO: needs implementation
def table_to_file(self, table, fname):

    pass
    # code to push an astropy table into a file
    # returns nothing

# TODO: needs documentation update
def get_time(self, specification=None):

    time = datetime.now()

    if specification is not None:
        try:
            return time.strftime(specification)
        except TypeError as e:
            print e.info
            print e.message

    return time.strftime("%Y-%m-%d-%H:%M:%S")