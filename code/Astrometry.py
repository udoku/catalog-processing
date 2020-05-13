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

# TODO: why is this a class?
class Astrometry:

    """

    This class defines methods for astrometrical analysis.

    Instantiated with a full_table.

    """

    def __init__(self, full_table):

        self.full_table = full_table

    # TODO: needs documentation updates
    def get_bin_mid(self, bins):

        """

        Input: list/array of bin delimit values

        Output: list consisting of the middle values of these

        """

        bins_middle_values = []

        for i,j in zip(bins[:-1],bins[1:]):

            bins_middle_values.append((i+j)/2)

        return(bins_middle_values)

    # Defines the gaussian
    # TODO: needs documentation updates
    def gaussian(self, x, a, b, c):

        return(a*np.exp(-(x-b)**2/(2*c**2)))


    def remove_masked(self, ap_column, header):

        """

        Takes an astroquery result and a column header (str) and returns a list not containing the masked elements.

        """

        data_list = []

        for i in ap_column[header]: # does this actually work?

            if type(i)== np.float64:

                data_list.append(i)

        return(data_list)

    # TODO: needs documentation update
    def make_hist(self, data_list):

        """

        Makes a histogram using np.histogram. Return a 2-tuple consisting of the output of np.histogram and the midpoints of the bins.

        """

        data_list = data_list[~data_list.mask]

        histogram = np.histogram(data_list.data, bins="auto")

        midpoints_of_bins = self.get_bin_mid(histogram[1])

        return((histogram, midpoints_of_bins))

    # TODO: needs documentation updates
    # TODO: is this used anywhere?
    def fit_gaussian(self, colname, full_column=False):

        """

        Fits a gaussian to a list of data points by automatically sorting them into

        bins, and fitting a gaussian using scipy's curve_fit.


        Input: List of data points


        Output: 3-tuple (a,b,c) consisting of the coefficients of the

                gaussian: a*e^{-(x-b)^2/(2c^2)}.


        Dependency: plt, numpy, scipy

        """


        if full_column:

            histogram = self.make_hist(colname)

        else:

            histogram = self.make_hist(self.full_table[colname])

        fit = scipy.optimize.curve_fit(self.gaussian, histogram[1], histogram[0][0])

        return((fit[0][0], fit[0][1], fit[0][2]))

    # TODO: needs documentation update
    # TODO: investigate options on HDBSCAN
    # (e.g. can we run several analyses in parallel?)
    def identify_clusters(self, table, columns, expected_clusters=None, verbose=False):

        # can include:
        # min_samples
        # min_cluster_size
        # cluster_selection_epsilon
        # allow_single_cluster
        # alpha

        #data = []
        #candidates_table = table[:]

        print("Organizing clustering data...")

        if verbose:
            print("Columns in supplied data table:")
            print(table.colnames)
            print("Columns to be used for clustering analysis:")
            print(columns)

        print("Initializing data structures...")
        data = []
        candidates_table = table[:]
        count = 0

        print("Selecting data...")
        for i in range(len(table[columns[0]])):
            datai=[]
            include_flag = True
            for c in columns:
                if not (table[c][i] != None):
                    include_flag = False
            if include_flag:
                for c in columns:
                    datai.append(table[c][i])
                data.append(datai)
            else:
                candidates_table.remove_row(i-count)
                count += 1

        if verbose:
            print("Columns in dataset: ")
            print(len(data[0]))
            print("Entries in dataset: ")
            print(len(data))

        print("Calculating clusters...")
        #clusterer = hdbscan.HDBSCAN(min_cluster_size=(len(data_A) / 10))
        #clusterer = hdbscan.HDBSCAN(allow_single_cluster=False, min_cluster_size=(len(data_A)/20))
        if expected_clusters:
            print("Iterating minimum cluster size...")
            
            iteration=1
            while(1):
                print(str(iteration) + "%")
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

        print("Detected " + str(max(membership) + 1) + " clusters")
        print("Clustering calculation complete.")#\n Detected " + str(max(membership)) + " clusters in a population of " + str(len(data_A)) " objects.")
        return candidates_table, membership

    # TODO: needs documentation update
    # TODO: make cluster filters more robust/read in from file?
    def process_clusters(self, candidates_table, membership, columns):
        clusters_data = []
        col_membership = Column(name="cluster_membership", data=membership)
        candidates_table.add_column(col_membership)
        memberphot = Photometry(candidates_table)

        for i in range(max(membership)):
            cut_string = "self.full_table['cluster_membership'] == " + str(i)
            members, _ = memberphot.apply_cut(cut_string)
            if abs(np.mean(members['pmra'])) - np.std(members['pmra']) > 0 and abs(np.mean(members['pmdec'])) - np.std(members['pmdec']) > 0:
                #members_trimmed = self.find_cluster_members(members, columns)
                clusters_data.append(members)

        return clusters_data

    # TODO: needs documentation update
    # takes a candidate cluster and pares down the members
    def find_cluster_members(self, cluster, columns):
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
            #print(len(cluster))
            loglikelihoods_i = []
            for i in range(len(cluster)):
                likelihood_i = 0
                #print(len(columns))
                for j in range(len(columns)):
                    #print("column: " + str(columns[j]))
                    #print("mean: " + str(means[j]))
                    #print("stdev: " + str(stdevs[j]))
                    likelihood_i += ((cluster[columns[j]][i] - means[j])/stdevs[j])**2
                loglikelihoods_i.append(likelihood_i)
                #print("likelihood: " + str(likelihood_i))
            loglikelihoods.append(np.mean(loglikelihoods_i))
            cluster.remove_row(loglikelihoods_i.index(max(loglikelihoods_i)))
        #print(zip(loglikelihoods[:-1],loglikelihoods[1:]))
        diffs = [y-x for x, y in zip(loglikelihoods[:-1], loglikelihoods[1:])]
        diffs2 = [y-x for x, y in zip(diffs[:-1], diffs[1:])]
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

    # TODO: needs documentation updates
    # TODO: shift to one function instead, and call with the desired parameter.
    def cluster_stats(self, cluster):
        pmra_hist = []
        pmra_best = 0
        pmdec_hist = []
        pmdec_best = 0
        parallax_hist = []
        parallax_best = 0
        ra_hist = []
        ra_best = 0
        dec_hist = []
        dec_best = 0

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

        print("best pmra: " + str(most))
        print("68 percent interval: " + str(low) + " to " + str(high))

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

        print("best pmdec: " + str(most))
        print("68 percent interval: " + str(low) + " to " + str(high))

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

        print("best parallax: " + str(most))
        print("68 percent interval: " + str(low) + " to " + str(high))

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

        print("best ra: " + str(most))
        print("68 percent interval: " + str(low) + " to " + str(high))

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

        print("best dec: " + str(most))
        print("68 percent interval: " + str(low) + " to " + str(high))

        fig = plt.figure()
        fig.clf()
        ax = fig.add_subplot(1,1,1)
        # apply axis labels
        ax.set_xlabel('dec')
        ax.set_ylabel('posterior probability')
        ax.plot(dec_vals, dec_actual)
        plt.show()