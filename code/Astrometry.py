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

from Photometry import *

import hdbscan

from getCI import *

from code import path, logger, out

def get_bin_mid(bins):
    """
    Given a list of bin delimiters, returns a list of bin midpoints.

    Arguments:
        bins [list of floats]: list/array of bin delimit values

    Return:
        list consisting of the bin midpoints
    """

    bins_middle_values = []

    for i,j in zip(bins[:-1],bins[1:]):

        bins_middle_values.append((i+j)/2)

    return(bins_middle_values)


def gaussian(x, a, b, c):
    """
    Defines the Gaussian function.

    Arguments:
        x [float]: value at which to evaluate the Gaussian
        a [float]: scaling factor for the function
        b [float]: mean for the distribution
        c [float]: standard deviation for the distribution

    Return:
        a * exp(-(x-b)^2/(2c^2))
    """

    return(a*np.exp(-(x-b)**2/(2*c**2)))


def remove_masked(table, colname):
    """
    Takes a CatalogTable and a column name and returns the unmasked elements.

    !!! not yet reviewed for functionality; original function from 2019 !!!

    Arguments:
        table [CatalogTable]: Table to access
        header [string]: column name

    Return:
        List of data values in the specified column that are not masked.
    """

    data_list = []

    # this will fail if "header" is not a valid column name
    try:
        for i in table[colname]: # does this actually work?

            if type(i)== np.float64:

                data_list.append(i)
    except:
        out("Specified column is not a member of the table:")
        out(colname)

    return(data_list)


def make_hist(data_list, bin_method="sqrt"):
    """
    Makes a histogram of data_list. 

    Available bin methods include:
     ‘auto’: Maximum of the ‘sturges’ and ‘fd’ estimators. Good all around.
     ‘fd’: Robust estimator that uses data variability and data size.
     ‘doane’: Improved Sturges’, works better with non-normal datasets.
     ‘rice’: Does not take variability into account, only data size.
     ‘sqrt’: Square root (of data size) estimator.

    Arguments:
        data_list [list of floats]: data to make the histogram
        bin_method [int or str]: method for creating bins. If an int, specifies the number of bins; otherwise uses the string methods above.

    Return:
        (histogram, bin midpoints)
    """

    data_list = data_list[~data_list.mask]

    histogram = np.histogram(data_list.data, bins=bin_method)

    midpoints_of_bins = get_bin_mid(histogram[1])

    return (histogram, midpoints_of_bins)


def fit_gaussian(table, colname):
    """
    Fits a gaussian to a list of data points by automatically sorting them into
    bins, and fitting a gaussian using scipy's curve_fit.

    Arguments:
        table [CatalogTable]; table to access for data
        colname [string]: column of the table to access

    Return: 
        3-tuple (a,b,c) consisting of the coefficients of the
        gaussian: a*e^{-(x-b)^2/(2c^2)}.

    """

    histogram = make_hist(table.table[colname])

    fit = scipy.optimize.curve_fit(gaussian, histogram[1], histogram[0][0])

    return((fit[0][0], fit[0][1], fit[0][2]))


# TODO: investigate options on HDBSCAN
# (e.g. can we run several analyses in parallel?)
def identify_clusters(table, columns, expected_clusters=None, verbose=False):
    """
    Uses HDBSCAN to identify candidate clusters by specified parameters.

    Arguments:
        table [CatalogTable]: table to access
        columns [list of strings]: list of columns to use for cluster detection
        expected_clusters [int]: Parameter for HDBSCAN; number of clusters expected. If not provided, no default is passed.
        verbose [bool]: verbose flag option, set to True for more diagnostics

    Return:
        CatalogTable of cluster candidates, and list of cluster membership

    """
    # HDBSCAN arguments can include:
    # min_samples
    # min_cluster_size
    # cluster_selection_epsilon
    # allow_single_cluster
    # alpha

    #data = []
    #candidates_table = table[:]

    out("Organizing clustering data...")

    if verbose:
        out("Columns in supplied data table:")
        out(table.table.colnames)
        out("Columns to be used for clustering analysis:")
        out(columns)

    out("Initializing data structures...")
    data = []
    candidates_table = table.table[:]
    count = 0

    out("Selecting data...")
    for i in range(len(table.table[columns[0]])):
        datai=[]
        include_flag = True
        for c in columns:
            if not (table.table[c][i] != None):
                include_flag = False
        if include_flag:
            for c in columns:
                datai.append(table.table[c][i])
            data.append(datai)
        else:
            candidates_table.remove_row(i-count)
            count += 1

    if verbose:
        out("Columns in dataset: ")
        out(len(data[0]))
        out("Entries in dataset: ")
        out(len(data))

    out("Calculating clusters...")
    #clusterer = hdbscan.HDBSCAN(min_cluster_size=(len(data_A) / 10))
    #clusterer = hdbscan.HDBSCAN(allow_single_cluster=False, min_cluster_size=(len(data_A)/20))
    if expected_clusters:
        out("Iterating minimum cluster size...")
        
        iteration=1
        while(1):
            out(str(iteration) + "%")
            min_cluster_size = abs(int(iteration * len(data) / 100))+ 1 # 1% of the data size
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=10)

            clusterer.fit(data)

            membership = clusterer.labels_

            if max(membership) < expected_clusters:
                break

            iteration += 1
    else:
    
        clusterer = hdbscan.HDBSCAN(min_cluster_size = 50, min_samples=10)
        clusterer.fit(data)
        membership = clusterer.labels_

    out("Detected " + str(max(membership) + 1) + " clusters")
    out("Clustering calculation complete.")#\n Detected " + str(max(membership)) + " clusters in a population of " + str(len(data_A)) " objects.")
    return CatalogTable(table.catalogs, candidates_table), membership


# TODO: make cluster filters more robust/read in from file?
def process_clusters(candidates_table, membership, columns):
    """
    Given candidate clusters, applies filter(s) to remove spurious clusters.

    Currently implemented filters:
        mean pmra and pmdec both more than 1 std dev away from zero

    Arguments:
        candidates_table [CatalogTable]: table of cluster candidates
        membership [list of ints]: list of candidate membership
        columns [list of strings]: list of columns, currently unused

    Return:
        CatalogTable of sources that passed the filters, with their cluster membership
    """

    clusters_data = []
    col_membership = Column(name="cluster_membership", data=membership)
    candidates_table.table.add_column(col_membership)

    for i in range(max(membership)):

        cut_string = "table.table['cluster_membership'] == " + str(i)
        members, _ = apply_cut(candidates_table, cut_string)

        if abs(np.mean(members.table['pmra'])) - np.std(members.table['pmra']) > 0 and abs(np.mean(members.table['pmdec'])) - np.std(members.table['pmdec']) > 0:

            #members_trimmed = find_cluster_members(table, members, columns)
            clusters_data.append(members)

    return CatalogTable(candidates_table.catalogs, clusters_data)

# TODO: needs documentation updates
# takes a candidate cluster and pares down the members
def find_cluster_members(cluster, columns):
    """
    Given a candidate cluster, eliminates the sources that decrease the likelihood of the cluster the most.

    Arguments:
        cluster [CatalogTable]: Candidate cluster
        columns [list of strings]: 

    Return:

    """

    loglikelihoods = []
    cluster_size = []
    cluster_members = []
    means = []
    stdevs = []
    for c in columns:
        means.append(np.mean(cluster[c]))
        stdevs.append(np.std(cluster[c]))
    while(len(cluster) > 0):
        cluster_size.append(len(cluster))
        cluster_members.append(cluster[:])
        #out(len(cluster))
        loglikelihoods_i = []
        for i in range(len(cluster)):
            likelihood_i = 0
            #out(len(columns))
            for j in range(len(columns)):
                #out("column: " + str(columns[j]))
                #out("mean: " + str(means[j]))
                #out("stdev: " + str(stdevs[j]))
                likelihood_i += ((cluster[columns[j]][i] - means[j])/stdevs[j])**2
            loglikelihoods_i.append(likelihood_i)
            #out("likelihood: " + str(likelihood_i))
        loglikelihoods.append(np.mean(loglikelihoods_i))
        cluster.remove_row(loglikelihoods_i.index(max(loglikelihoods_i)))
    #out(zip(loglikelihoods[:-1],loglikelihoods[1:]))
    diffs = [y-x for x, y in zip(loglikelihoods[:-1], loglikelihoods[1:])]
    #diffs2 = [y-x for x, y in zip(diffs[:-1], diffs[1:])]
    #diffsmedian = np.median(diffs)

    #index = len(diffs) - 1
    #while(1):
    #    if abs(diffs[index]) > 2 * abs(diffsmedian):
    #        break
    #    index -= 1

    #fig = plt.figure()
    #fig.clf()

    # create subfigure to plot lists
    #ax = fig.add_subplot(1,1,1)
    #ax.scatter(cluster_size, loglikelihoods)
    #plt.vlines(cluster_size[index],min(loglikelihoods),max(loglikelihoods))
    #plt.show()

    index = len(diffs) - 1
    while(1):
        if abs(diffs[index]) > abs(np.mean(diffs)):
            break
        index -= 1

    fig = plt.figure()
    fig.clf()

    # create subfigure to plot lists
    ax = fig.add_subplot(1,1,1)
    ax.scatter(cluster_size, loglikelihoods)
    plt.vlines(cluster_size[index],min(loglikelihoods),max(loglikelihoods))
    plt.show()
    
    fig = plt.figure()
    fig.clf()

    # create subfigure to plot lists
    ax = fig.add_subplot(1,1,1)
    ax.scatter(cluster_size[1:], diffs)
    #plt.ylim(-0.01,0)
    plt.show()

    return cluster_members[index]


# TODO: shift to one function instead, and call with the desired parameter.
def cluster_stats(cluster, param):
    """
    Given a candidate cluster and a parameter, finds summary stats and plots that parameter for the cluster.

    Arguments:
        cluster [CatalogTable]: candidate cluster
        param [string]: parameter to analyze

    Return:
        [none]
    """

    param_hist = []
    
    try:
        param_error = param + "_error"
        param_data = cluster.table[param]
        param_error_data = cluster.table[param_error]
    except KeyError as e:
        out("Passed parameter does not exist in cluster table.")
        out(e.message)

    minval = min(param_data)
    maxval = max(param_data)
    param_vals = np.linspace(minval, maxval, ,num=5*len(param_data))

    for param in param_vals:
        loglikelihood = 0
        for index in range(len(param_data)):
            loglikelihood += 0.5 * ((param - param_data[index])/param_error_data[index])**2 - param_error_data[index]
        param_hist.append(loglikelihood)

    param_hist /= min(param_hist)
    param_actual = []
    for value in param_hist:
        param_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(param_actual, param_vals)
    param_actual /= normalize

    most, low, high, _ = getCI(param_vals, param_actual, 0.68)

    out("best " + param + ": " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel(param)
    ax.set_ylabel('posterior probability')
    ax.plot(param_vals, param_actual)
    plt.show()
    
    """
    #pmra_best = 0
    #pmdec_hist = []
    #pmdec_best = 0
    #parallax_hist = []
    #parallax_best = 0
    #ra_hist = []
    #ra_best = 0
    #dec_hist = []
    #dec_best = 0

    pmra_vals = np.linspace(min(cluster['pmra']),max(cluster['pmra']),num=5*len(cluster['pmra']))
    pmdec_vals = np.linspace(min(cluster['pmdec']),max(cluster['pmdec']),num=5*len(cluster['pmdec']))
    parallax_vals = np.linspace(min(cluster['parallax']),max(cluster['parallax']),num=5*len(cluster['parallax']))
    ra_vals = np.linspace(min(cluster['gaia_ra']),max(cluster['gaia_ra']),num=5*len(cluster['gaia_ra']))
    dec_vals = np.linspace(min(cluster['gaia_dec']),max(cluster['gaia_dec']),num=5*len(cluster['gaia_dec']))

    for pmra in pmra_vals:
        loglikelihood = 0
        for index in range(len(cluster['pmra'])):
            loglikelihood += 0.5 * ((pmra - cluster['pmra'][index])/cluster['pmra_error'][index])**2 - cluster['pmra_error'][index]
        pmra_hist.append(loglikelihood)

    pmra_hist /= min(pmra_hist)
    pmra_actual = []
    for value in pmra_hist:
        pmra_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(pmra_actual, pmra_vals)
    pmra_actual /= normalize

    most, low, high, _ = getCI(pmra_vals, pmra_actual, 0.68)

    out("best pmra: " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel('pmra')
    ax.set_ylabel('posterior probability')
    ax.plot(pmra_vals, pmra_actual)
    plt.show()

    for pmdec in pmdec_vals:
        loglikelihood = 0
        for index in range(len(cluster['pmdec'])):
            loglikelihood += 0.5 * ((pmdec - cluster['pmdec'][index])/cluster['pmdec_error'][index])**2
        pmdec_hist.append(loglikelihood)
        if loglikelihood < min(pmdec_hist):
            pmra_best = pmdec

    pmdec_hist /= min(pmdec_hist)
    pmdec_actual = []
    for value in pmdec_hist:
        pmdec_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(pmdec_actual, pmdec_vals)
    pmdec_actual /= normalize

    most, low, high, _ = getCI(pmdec_vals, pmdec_actual, 0.68)

    out("best pmdec: " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel('pmdec')
    ax.set_ylabel('posterior probability')
    ax.plot(pmdec_vals, pmdec_actual)
    plt.show()

    for parallax in parallax_vals:
        loglikelihood = 0
        for index in range(len(cluster['parallax'])):
            loglikelihood += 0.5 * ((parallax - cluster['parallax'][index])/cluster['parallax_error'][index])**2 - cluster['parallax_error'][index]
        parallax_hist.append(loglikelihood)

    parallax_hist /= min(parallax_hist)
    parallax_actual = []
    for value in parallax_hist:
        parallax_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(parallax_actual, parallax_vals)
    parallax_actual /= normalize

    most, low, high, _ = getCI(parallax_vals, parallax_actual, 0.68)

    out("best parallax: " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel('parallax')
    ax.set_ylabel('posterior probability')
    ax.plot(parallax_vals, parallax_actual)
    plt.show()

    for ra in ra_vals:
        loglikelihood = 0
        for index in range(len(cluster['gaia_ra'])):
            loglikelihood += 0.5 * ((ra - cluster['gaia_ra'][index])/cluster['ra_error'][index])**2 - cluster['ra_error'][index]
        ra_hist.append(loglikelihood)

    ra_hist /= min(ra_hist)
    ra_actual = []
    for value in ra_hist:
        ra_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(ra_actual, ra_vals)
    ra_actual /= normalize

    most, low, high, _ = getCI(ra_vals, ra_actual, 0.68)

    out("best ra: " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel('ra')
    ax.set_ylabel('posterior probability')
    ax.plot(ra_vals, ra_actual)
    plt.show()

    for dec in dec_vals:
        loglikelihood = 0
        for index in range(len(cluster['gaia_dec'])):
            loglikelihood += 0.5 * ((dec - cluster['gaia_dec'][index])/cluster['dec_error'][index])**2 - cluster['dec_error'][index]
        dec_hist.append(loglikelihood)

    dec_hist /= min(dec_hist)
    dec_actual = []
    for value in pmra_hist:
        dec_actual.append(np.exp(-1.0 * value))

    normalize = np.trapz(dec_actual, dec_vals)
    dec_actual /= normalize

    most, low, high, _ = getCI(dec_vals, dec_actual, 0.68)

    out("best dec: " + str(most))
    out("68 percent interval: " + str(low) + " to " + str(high))

    fig = plt.figure()
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    # apply axis labels
    ax.set_xlabel('dec')
    ax.set_ylabel('posterior probability')
    ax.plot(dec_vals, dec_actual)
    plt.show()
    """