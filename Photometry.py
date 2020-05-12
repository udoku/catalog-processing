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

class Photometry:

    def __init__(self, full_table, full_column=False):

        self.full_table = full_table

        self.none = None





    def apply_cut(self, cut):

        """
        Used to remove rows of self.full_table that don't meet a certain cut.

        Arguments:
            cut [string]: a string representing a cut, for example:
                "J-H > 0.7"
                "J-K > H-K"

        Returns:
            cut_true [astropy table]: table of sources meeting the cut
            cut_false [astropy table]: table of sources not meeting the cut

        """
        """

        J-H vs H-K

        J-K vs K-W1

        J-K vs K-W2

        J-K vs K-W3

        W1-W2 vs W2-W3

        W1-W2 vs W3-W4

        G-K vs K-W2

        BP-RP vs J-K

        """

        # table of color values. Allows shorthand filter names to be used and then replaced by column names
        colors =  {"G": """self.full_table["phot_g_mean_mag"]""",
                   "H": """self.full_table["h_m"]""",
                   "J": """self.full_table["j_m"]""",
                   "K": """self.full_table["k_m"]""",
                   "W1": """self.full_table["w1mpro"]""",
                   "W2": """self.full_table["w2mpro"]""",
                   "W3": """self.full_table["w3mpro"]""",
                   "W4": """self.full_table["w4mpro"]""",
                   "BP-RP": """self.full_table["bp_rp"]"""}

        # replaces all instances of a color keyword with the full_table column for that color
        for color in colors.keys():
            cut = cut.replace(color, colors[color])

        # evaluates section of the table that passes the cut and section that does not
        cut_true = self.full_table[pd.eval(cut)]
        cut_false = self.full_table[~pd.eval(cut)]

        # return cut_true and cut_false as a tuple
        return((cut_true, cut_false))



    def cut_list(self, iterable, cut, collist=None):

        """

        Returns tuple of two lists the first corresponding to the values where the cut is true and the other to the values where the cut is false.

        """

        cut_true =[]
        cut_false = []

        if collist == None:
            for val in iterable:
                if pd.eval(cut):
                    cut_true.append(val)
                else:
                    cut_false.append(val)
        else:
            for val, col_val in zip(iterable, col):
                if pd.eval(cut):
                    cut_true.append(val)
                else:
                    cut_false.append(val)



    def variability_cut(self):

        """
        Returns a copy of the full_table with systems that satisfy the variability cut.
        """

        cut_table = self.full_table.copy()

        mask_shape = []

        variability = [np.sqrt( g_n_obs ) / mean_flux_over_error for g_n_obs, mean_flux_over_error in zip(self.full_table["phot_g_n_obs"], self.full_table["phot_g_mean_flux_over_error"])]

        for g, v in zip(self.full_table["phot_g_mean_mag"], variability):

            if (v > 0.05 and g < 18):
                mask_shape.append(False)
            else:
                mask_shape.append(True)

        cut_table["gaia_designation"].mask = mask_shape

        cut_table.remove_rows(np.where(mask_shape))

        return(cut_table)


    def gaia_color_mag_cut(self):

        """
        Returns a copy of the full_table with systems not satisfying the following photometric cuts:
        
        """

        M_g = [g + 5 - np.log10(1000/p) -10 for g,p in zip(self.full_table["phot_g_mean_mag"], self.full_table["parallax"])]

        mask_shape = []

        for col, m_g in zip(self.full_table["bp_rp"], M_g):

            if (m_g < 2.46*col + 2.76) and (0.3 < col < 1.8):

                mask_shape.append(False)

            elif (m_g < 2.8*col + 2.16) and (1.8 < col):

                mask_shape.append(False)

            elif (m_g > 2.14*col - 0.57) and (0.5 < col < 1.2):

                mask_shape.append(False)

            elif (m_g > 1.11*col + 0.66) and (1.2 < col < 3):

                mask_shape.append(False)

            else:

                mask_shape.append(True)

        cut_table = self.full_table.copy()

        cut_table["gaia_designation"].mask = mask_shape

        cut_table.remove_rows(np.where(mask_shape))

        return(cut_table)


    def keep_reddest(self, join_method="and"):
        """
        Returns a list of sources that are redder than the reddest photospheres. Default criteria applied:
        W1-W2 > 0.15
        K-W2 > 0.3
        K-W1 > 0.2
        J-W2 > 1.3
        J-W1 > 1.1
        H-K > 0.4
        J-K > 0.8
        J-H > 0.7

        Inputs:
            join_method [string, optional, default: "and"]: string designating how to combine the criteria. Typical values would be "and" and "or".

        Outputs:
            cut_true, cut_false [tuple of astropy tables]: tables which satisfy and fail to satisfy the specified criteria, respectively.
        """

        join_argument = " " + join_method + " "

        self.cut_string = "" 

        master = Tk()
        master.title="criteria menu"

        def var_states():
            if var_w1w2.get():
                val = val_w1w2.get()
                print("W1-W2 > " + val)
                self.cut_string += "W1-W2 > " + val + join_argument
            if var_kw2.get():
                val = val_kw2.get()
                print("K-W2 > " + val)
                self.cut_string += "K-W2 > " + val + join_argument
            if var_kw1.get():
                val = val_kw1.get()
                print("K-W1 > " + val)
                self.cut_string += "K-W1 > " + val + join_argument
            if var_jw2.get():
                val = val_jw2.get()
                print("J-W2 > " + val)
                self.cut_string += "J-W2 > " + val + join_argument
            if var_jw1.get():
                val = val_jw1.get()
                print("J-W1 > " + val)
                self.cut_string += "J-W1 > " + val + join_argument
            if var_hk.get():
                val = val_hk.get()
                print("H-K > " + val)
                self.cut_string += "H-K > " + val + join_argument
            if var_jk.get():
                val = val_jk.get()
                print("J-K > " + val)
                self.cut_string += "J-K > " + val + join_argument
            if var_jh.get():
                val = val_jh.get()
                print("J-H > " + val)
                self.cut_string += "J-H > " + val + join_argument
            if var_custom.get():
                val = val_custom.get()
                print(val)
                self.cut_string += val
        

        Label(master, text="Criteria to use:").grid(row=0, sticky=W)

        var_w1w2 = BooleanVar()
        var_w1w2.set(True)
        Checkbutton(master, text="W1-W2 > ", variable=var_w1w2).grid(row=1, column=1, sticky=W)
        val_w1w2 = StringVar(master, value="0.15")
        Entry(master, textvariable=val_w1w2).grid(row=1, column=2, sticky=W)

        var_kw2 = BooleanVar()
        var_kw2.set(True)
        Checkbutton(master, text="K-W2 > ", variable=var_kw2).grid(row=2, column=1, sticky=W)
        val_kw2 = StringVar(master, value="0.3")
        Entry(master, textvariable=val_kw2).grid(row=2, column=2, sticky=W)

        var_kw1 = BooleanVar()
        var_kw1.set(True)
        Checkbutton(master, text="K-W1 > ", variable=var_kw1).grid(row=3, column=1, sticky=W)
        val_kw1 = StringVar(master, value="0.2")
        Entry(master, textvariable=val_kw1).grid(row=3, column=2, sticky=W)

        var_jw2 = BooleanVar()
        var_jw2.set(True)
        Checkbutton(master, text="J-W2 > ", variable=var_jw2).grid(row=4, column=1, sticky=W)
        val_jw2 = StringVar(master, value="1.3")
        Entry(master, textvariable=val_jw2).grid(row=4, column=2, sticky=W)

        var_jw1 = BooleanVar()
        var_jw1.set(True)
        Checkbutton(master, text="J-W1 > ", variable=var_jw1).grid(row=5, column=1, sticky=W)
        val_jw1 = StringVar(master, value="1.1")
        Entry(master, textvariable=val_jw1).grid(row=5, column=2, sticky=W)

        var_hk = BooleanVar()
        var_hk.set(True)
        Checkbutton(master, text="H-K > ", variable=var_hk).grid(row=6, column=1, sticky=W)
        val_hk = StringVar(master, value="0.4")
        Entry(master, textvariable=val_hk).grid(row=6, column=2, sticky=W)

        var_jk = BooleanVar()
        var_jk.set(True)
        Checkbutton(master, text="J-K > ", variable=var_jk).grid(row=7, column=1, sticky=W)
        val_jk = StringVar(master, value="0.8")
        Entry(master, textvariable=val_jk).grid(row=7, column=2, sticky=W)

        var_jh = BooleanVar()
        var_jh.set(True)
        Checkbutton(master, text="J-H > ", variable=var_jh).grid(row=8, column=1, sticky=W)
        val_jh = StringVar(master, value="0.7")
        Entry(master, textvariable=val_jh).grid(row=8, column=2, sticky=W)

        var_custom = BooleanVar()
        var_custom.set(False)
        Checkbutton(master, text="Custom: ", variable = var_custom).grid(row=9, column=1, sticky=W)
        val_custom = StringVar(master, value="")
        Entry(master, textvariable=val_custom).grid(row=9, column=2, sticky=W)
        

        #Button(master, text='Cancel', command=master.destroy).grid(row=9, sticky=W, pady=4)
        Button(master, text='Enter', command=lambda:[master.destroy(),var_states()]).grid(row=10, sticky=W, pady=4)

        master.mainloop()

        self.cut_string = self.cut_string.strip().strip('and').strip()

        print(self.cut_string)
        return self.apply_cut(self.cut_string)


    def keep_likely_disks(self, join_method="and"):
        """
        Returns a list of sources that are so red, they cannot be explained by a reddened photosphere (and are thus likely to have disks). Default criteria applied:
        W3-W4 > 1
        W2-W3 > 0.8 and W2-W3 < 4
        W1-W2 > 0.35
        K-W4 > 1.2 and K-W4 < 10
        K-W3 > 0.9 and K-W3 < 7
        K-W2 > 0.6
        K-W1 > 0.3
        H-K > 1.5 and K-W2 > 2.8

        Inputs:
            join_method [string, optional, default: "and"]: string designating how to combine the criteria. Typical values would be "and" and "or".

        Outputs:
            cut_true, cut_false [tuple of astropy tables]: tables which satisfy and fail to satisfy the specified criteria, respectively.
        """

        join_argument = " " + join_method + " "

        self.cut_string = "" 

        master = Tk()
        master.title="criteria menu"

        def var_states():
            if var_w1w2.get():
                val = val_w1w2.get()
                print("W1-W2 > " + val)
                self.cut_string += "W1-W2 > " + val + join_argument
            if var_w2w3g.get():
                val = val_w2w3g.get()
                print("W2-W3 > " + val)
                self.cut_string += "W2-W3 > " + val + join_argument
            if var_w2w3l.get():
                val = val_w2w3l.get()
                print("W2-W3 < " + val)
                self.cut_string += "W2-W3 < " + val + join_argument
            if var_w3w4.get():
                val = val_w3w4.get()
                print("W3-W4 > " + val)
                self.cut_string += "W3-W4 > " + val + join_argument
            if var_kw1.get():
                val = val_kw1.get()
                print("K-W1 > " + val)
                self.cut_string += "K-W1 > " + val + join_argument
            if var_kw2.get():
                val = val_kw2.get()
                print("K-W2 > " + val)
                self.cut_string += "K-W2 > " + val + join_argument
            if var_kw3g.get():
                val = val_kw3g.get()
                print("K-W3 > " + val)
                self.cut_string += "K-W3 > " + val + join_argument
            if var_kw3l.get():
                val = val_kw3l.get()
                print("K-W3 < " + val)
                self.cut_string += "K-W3 < " + val + join_argument
            if var_kw4g.get():
                val = val_kw4g.get()
                print("K-W4 > " + val)
                self.cut_string += "K-W4 > " + val + join_argument
            if var_kw4l.get():
                val = val_kw4l.get()
                print("K-W4 < " + val)
                self.cut_string += "K-W4 < " + val + join_argument
            if var_hk.get():
                val = val_hk.get()
                print("H-K > " + val)
                self.cut_string += "H-K > " + val + join_argument
            if var_custom.get():
                val = val_custom.get()
                print(val)
                self.cut_string += val
        

        Label(master, text="Criteria to use:").grid(row=0, sticky=W)

        var_w1w2 = BooleanVar()
        var_w1w2.set(True)
        Checkbutton(master, text="W1-W2 > ", variable=var_w1w2).grid(row=1, column=1, sticky=W)
        val_w1w2 = StringVar(master, value="0.35")
        Entry(master, textvariable=val_w1w2).grid(row=1, column=2, sticky=W)

        var_w2w3g = BooleanVar()
        var_w2w3g.set(True)
        Checkbutton(master, text="W2-W3 > ", variable=var_w2w3g).grid(row=2, column=1, sticky=W)
        val_w2w3g = StringVar(master, value="0.8")
        Entry(master, textvariable=val_w2w3g).grid(row=2, column=2, sticky=W)

        var_w2w3l = BooleanVar()
        var_w2w3l.set(True)
        Checkbutton(master, text="W2-W3 < ", variable=var_w2w3l).grid(row=3, column=1, sticky=W)
        val_w2w3l = StringVar(master, value="4.0")
        Entry(master, textvariable=val_w2w3l).grid(row=3, column=2, sticky=W)

        var_w3w4 = BooleanVar()
        var_w3w4.set(True)
        Checkbutton(master, text="W3-W4 > ", variable=var_w3w4).grid(row=4, column=1, sticky=W)
        val_w3w4 = StringVar(master, value="1.0")
        Entry(master, textvariable=val_w3w4).grid(row=4, column=2, sticky=W)

        var_kw1 = BooleanVar()
        var_kw1.set(True)
        Checkbutton(master, text="K-W1 > ", variable=var_kw1).grid(row=5, column=1, sticky=W)
        val_kw1 = StringVar(master, value="0.3")
        Entry(master, textvariable=val_kw1).grid(row=5, column=2, sticky=W)

        var_kw2 = BooleanVar()
        var_kw2.set(True)
        Checkbutton(master, text="K-W2 > ", variable=var_kw2).grid(row=6, column=1, sticky=W)
        val_kw2 = StringVar(master, value="2.8")
        Entry(master, textvariable=val_kw2).grid(row=6, column=2, sticky=W)

        var_kw3g = BooleanVar()
        var_kw3g.set(True)
        Checkbutton(master, text="K-W3 > ", variable=var_kw3g).grid(row=7, column=1, sticky=W)
        val_kw3g = StringVar(master, value="0.9")
        Entry(master, textvariable=val_kw3g).grid(row=7, column=2, sticky=W)

        var_kw3l = BooleanVar()
        var_kw3l.set(True)
        Checkbutton(master, text="K-W3 < ", variable=var_kw3l).grid(row=8, column=1, sticky=W)
        val_kw3l = StringVar(master, value="7.0")
        Entry(master, textvariable=val_kw3l).grid(row=8, column=2, sticky=W)

        var_kw4g = BooleanVar()
        var_kw4g.set(True)
        Checkbutton(master, text="K-W4 > ", variable=var_kw4g).grid(row=9, column=1, sticky=W)
        val_kw4g = StringVar(master, value="1.2")
        Entry(master, textvariable=val_kw4g).grid(row=9, column=2, sticky=W)

        var_kw4l = BooleanVar()
        var_kw4l.set(True)
        Checkbutton(master, text="K-W4 < ", variable=var_kw4l).grid(row=10, column=1, sticky=W)
        val_kw4l = StringVar(master, value="10.0")
        Entry(master, textvariable=val_kw4l).grid(row=10, column=2, sticky=W)

        var_hk = BooleanVar()
        var_hk.set(True)
        Checkbutton(master, text="H-K > ", variable=var_hk).grid(row=11, column=1, sticky=W)
        val_hk = StringVar(master, value="1.5")
        Entry(master, textvariable=val_hk).grid(row=11, column=2, sticky=W)

        var_custom = BooleanVar()
        var_custom.set(False)
        Checkbutton(master, text="Custom: ", variable = var_custom).grid(row=12, column=1, sticky=W)
        val_custom = StringVar(master, value="")
        Entry(master, textvariable=val_custom).grid(row=12, column=2, sticky=W)
        

        #Button(master, text='Cancel', command=master.destroy).grid(row=9, sticky=W, pady=4)
        Button(master, text='Enter', command=lambda:[master.destroy(),var_states()]).grid(row=13, sticky=W, pady=4)

        master.mainloop()

        self.cut_string = self.cut_string.strip().strip('and').strip()

        print(self.cut_string)
        return self.apply_cut(self.cut_string)


    def remove_star_forming_galaxies(self):
        """
        Returns a copy of full_table without star forming galaxies that may mimic young stars. Criteria used:
        IF W1 > 14, keep only if
            W2-W3 < 2.3 AND W1-W2 > 1  AND W1-W2 > 0.46*(W2-W3) - 0.78 AND W1 < 1.8*(W1-W3) + 4.1 

        Outputs:
            cut_true, cut_false [tuple of astropy tables]: tables which satisfy and fail to satisfy the specified criteria, respectively.
        """

        cut = "(W1 > 14 and W2-W3 < 2.3 AND W1-W2 > 1  AND W1-W2 > 0.46*(W2-W3) - 0.78 AND W1 < 1.8*(W1-W3) + 4.1) or W1 < 14"

        return self.apply_cut(cut)


    def subtract_cols(self, colname1, colname2):

        """

        Returns a list consisting of the difference between the values of column colname1 and colname2.

        """

        return([i-j for i,j in zip(self.full_table[colname1], self.full_table[colname2])])