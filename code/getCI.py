#getCI.py
# Returns shortest confidence interval
# NOTE: Must avoid rounding issues

import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import pdb

def getcdf(x,px):
    cdf=cumtrapz(px,x=x,initial=0)
    return cdf

def getCI(x,px,ci):
    mostlike=x[np.argmax(px)]
    cdf=getcdf(x,px)
    min_index=0
    max_index=np.min(np.where((cdf[-1]-cdf)<ci)[0])
    # Compute the width of each low/high pair
    low_index=np.zeros(max_index-min_index)*np.nan
    high_index=np.zeros(max_index-min_index)*np.nan
    all_widths=np.zeros(max_index-min_index)
    for ii in range(min_index,max_index):
        low_index[ii]=ii
        high_index[ii]=np.min(np.where(cdf-cdf[ii]>ci)[0])
        all_widths[ii]=x[int(high_index[ii])]-x[int(low_index[ii])]
    # Determine all pairs with width within dx spacing (error tolerance)
    width_err_tol=x[1]-x[0] # tolerance in width error (dx)
    best_index=np.argmin(all_widths)
    best_indices=np.where(np.abs(all_widths-np.min(all_widths))<width_err_tol)[0]
    # Choose the best pair as the median of all such pairs
    best_index=np.int(np.median(best_indices))
    low=x[int(low_index[int(best_index)])]
    high=x[int(high_index[int(best_index)])]
    return (mostlike,low,high,cdf)

def get_high(cdf_interpf,low,ci):
    myfun = lambda x: np.abs(cdf_interpf(x)-cdf_interpf(low)-ci)
    best_high=minimize(myfun,low)
    return best_high