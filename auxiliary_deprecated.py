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


"""

Here is an outline of the final notebook



Class catalogueProcessing:

    Instantiated with:

        HAe/Be coords/list

        physical or angular radius



    > Methods for querying



    > Methods for xmatching



    > Method for generating one table with Gaia, 2mass, and allwise data (full_table)



Class astrometry:

    Instantiated with:

        full_table



    > Method to generate gaussian fit

    > Method to generate 2d gaussian fit

    > Method to do DBSCAN



Class plotting

    Instantiated with full_table



    > Method for plotting astrometry

        > includes the gaussians

    > Method for plotting photometry



"""

class catalogProcessing:

    """

    Instantiated with

        > a string representing HAeBe coordinates in degrees.

        > radius: number, either physical or angular radius (depending on if parallax is given).

        > parallax (optional): number, in milliarcseconds. If given, radius given is assumed to be physical, and will

          be recalculated into an angular one.

    """

    def __init__(self, ra, dec, radius, parallax=None):

        #self.skycoord = SkyCoord(HAeBe_coord, unit=(u.hourangle, u.deg))
        self.skycoord = SkyCoord(ra, dec)

        self.radius = radius

        self.parallax = parallax



        if self.parallax == None:

            pass

        else:

            self.radius = self.get_radius()


        print("Querying Gaia...")
        self.gaia = self.gaia_query()
        print("Done querying Gaia!\n")
        print("Querying 2MASS...")
        #self.tmass = self.ir_query("fp_psc")
        print("Done querying 2MASS!\n")
        print("Querying AllWISE...")
        #self.allwise = self.ir_query("allwise_p3as_psd")
        print("Done querying AllWISE!\n")
        # other catalogs to be added



    def get_radius(self):

        """
        Input: parallax in millliarcseconds, the desired physical search radius in parsec.
        Output: angular search radius in degrees.
        """

        # since parallax is given in milliarcseconds, we have to multiply the reciprocal of the parallax by 1000 to get distance in parsecs
        # in general, distance (pc) ~ 1 / parallax (arcsec)
        dist = 1000/(self.parallax)

        angular_radius_as = ((self.radius)*2.06265*10**5)/dist # magic number: number of AU in a parsec (why??)

        return(angular_radius_as/3600)





    def gaia_query(self):

        """

        Performs a query on gaia DR2.

        Input: coordinates given as a 2-tuple (305.11765944213,41.36428841871) for instance, and a radius (float) given in degrees.



        Output: astropy table with query result.



        Note: "Exception: 500" tends to happen when the coordinates used in search_string is not formatted correctly.

        """







        #try:

        #    iterator = iter(radius)

        #except TypeError:





        ra = self.skycoord.ra.degree

        dec = self.skycoord.dec.degree

        #print((ra, dec, self.radius))

        #print("Done 1")

        #print(coords.ra, coords.dec)

        search_string = "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{0},{1},{2}))=1;".format(ra, dec, self.radius)





        #table = Gaia.load_tables(only_names=True)

        #print("Done 2")

        #print(search_string)

        # Angles are in degrees.

        # "gaierror: [Errno 11001] getaddrinfo failed" is typically raised when there is no internet connection.

        #try:

        print("Creating query...")
        job = Gaia.launch_job_async(search_string, dump_to_file=False)

        #except gaierror as e:

        #    if str(e) != "[Errno 11001] getaddrinfo failed":

        #        raise

        #    else:

        #        print("This error is typically raised when there is no internet connection.")

        #        raise
        print("Retrieving results...")
        g = job.get_results()

        #print("Done 3")

        #table_list  .append(g)



        return(g)





    def ir_query(self, ircat_name, view_adql=False):

        """

        Performs a TAP+ query to irsa.ipac and returns the results.



        ircat_name: string, name of the catalogue according to the irsa.ipac TAP specs

        center_coords: coordinates in SkyCoord format.

        radius: search radius in degrees

        view_adql: False by default. If True,prints the adql of the query generated.



        -----------------------

        Examples:

        "allwise_p3as_psd": AllWISE Source Catalog

        "fp_psc": 2MASS Point Source Catalog

        "glimpse_s07": GLIMPSE I Spring 07 Catalog (Spitzer)

        "cosmos_phot": COSMOS Photometry Catalog

        "iraspsc": IRAS Point Source Catalog

        -----------------------



        Main reason for this function is to have basis for easily incorporating TAP queries in a general query function later.

        """







        if view_adql is True:

            print("SELECT * FROM {} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))=1 ".format(ircat_name, str(self.skycoord.ra.degree), str(self.skycoord.dec.degree), self.radius))



        ipac = TapPlus(url="https://irsa.ipac.caltech.edu/TAP/")


        print("Creating query...")
        job = ipac.launch_job_async("SELECT * FROM {} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))=1 ".format(ircat_name, str(self.skycoord.ra.degree), str(self.skycoord.dec.degree), self.radius))


        print("Retrieving results...")
        return(job.get_results())





    def xmatch(self, cat1, cat2, colnames1, colnames2, rad=-1, show_dynamic_radii=True):

        """

        valid surveynames: "2mass", "allwise", "gaia"



        Crossmatch cat1 with cat2 at a maximum xmatch radius of rad.

        Uses SkyCoord.match_to_catalog_sky() method which implements a nearest neighbour algorithm.



        cat1, cat2: full query results (astropy tables resulting from querying catalogs)



        colnames1, colnames2: list of the column names of ra dec and respective errors in cat1, cat2. Either [ra, dec] or [ra, dec, err_ra, err_dec] for



        Returns: list of lists [cat1_idx, new_idx, new_d2d] or if show_dynamic_radii=True [cat1_idx, new_idx, new_d2d, xmatch_rad]

            cat1_idx

            cat1_idx: index of the element in cat2 that is nearest to the corresponding element in cat1

                              (cat1_idx[0] = 69 means that the first element in cat1 is matched towards the 69th element in cat2)



        """





        surveys_coord_colnames = {"2mass":["ra", "dec", "err_maj", "err_min"],

                                  "allwise":["ra", "dec", "sigra", "sigdec"],

                                  "gaia":["ra", "dec", "ra_error", "dec_error"]}



        ra1 = cat1[colnames1[0]]

        dec1 = cat1[colnames1[1]]

        ra2 = cat2[colnames2[0]]

        dec2 = cat2[colnames2[1]]

        #print(ra1)

        #print(dec1)



        if rad < 0:

            ra_err1 = cat1[colnames1[2]]

            dec_err1 = cat1[colnames1[3]]



            ra_err2 = cat2[colnames2[2]]

            dec_err2 = cat2[colnames2[3]]



        skycoord1 = SkyCoord(ra1, dec1)

        skycoord2 = SkyCoord(ra2, dec2)





        idx, d2d, d3d = skycoord1.match_to_catalog_sky(skycoord2)



        #print(idx, d2d)



        #Imposes constraint on the xmatch radius.

        if rad < 0:

            constraint = "dynamic"

            #print("True")

        else:

            constraint = d2d < rad*u.arcsec







        #print(idx, len(idx))

        """

        Iterates through the index (index[0] = 69 means that the first element in cat1 is matched towards the 69th element in cat2),

        constraint, and

        d2d (2D angular xmatch distance.)

        """



        new_idx = idx

        new_d2d = d2d



        if constraint == "dynamic":



            xmatch_rad = []

            cat1_idx = []



            i = 0

            j = 0

            for indpair in zip(idx, d2d):

                #print(i)

                #Done: change the RA error to RAerror*cos(DEC) [not cos(DECerror)]



                #error_radius1 = min(ra_err1[j]*np.cos(dec1[j]), dec_err1[j])

                #error_radius2 = min(ra_err2[j]*np.cos(dec2[j]), dec_err2[j])

                error_radius1 = max(ra_err1[j]*np.cos(dec1[j]), dec_err1[j])

                error_radius2 = max(ra_err2[j]*np.cos(dec2[j]), dec_err2[j])

                #match_constraint = np.sqrt(error_radius_min**2 + error_radius_max**2)

                match_constraint = error_radius1 + error_radius2

                #print(type(match_constraint))





                if show_dynamic_radii is True:

                    xmatch_rad.append(match_constraint)

                else:

                    pass

                dist = indpair[1].arcsec

                #dist = coordinates.Angle(dist*u.arcsec)



                #print(type(dist), type(coordinates.Angle(match_constraint*u.arcsec)))

                #const = coordinates.Angle(match_constraint*u.arcsec)

                #print(dist, const)



                if dist > match_constraint:

                    #print(dist, match_constraint)

                    #print(i[1].arcsec, coordinates.Angle(match_constraint*u.arcsec))

                    #print(new_idx[j])

                    new_idx = np.delete(new_idx, j)

                    new_d2d = np.delete(new_d2d, j)

                    j-=1

                else:

                    cat1_idx.append(i)

                    #print(dist, match_constraint)

                #print(j)



                j+=1

                i+=1









        else:

            j = 0 # what does j do?

            for i in zip(idx, constraint, d2d):

                #print(i[1,d2d])



                # Checks if the
                # if the... what?

                if i[1] == False: #constraint = false

                    # if the constraint is false, then we delete the entry

                    #print(idx)

                    new_idx = np.delete(idx, j)

                    new_d2d = np.delete(d2d, j)

                    #print(i[2]) # print d2d[i]

                    # Removes a step from the iterator in case an element was removed.
                    # (why?)

                    j-=1

                j+=1



        if show_dynamic_radii is True:

            return([cat1_idx, new_idx, new_d2d, xmatch_rad])

        else:

            return([cat1_idx, new_idx, new_d2d])



    def empty_combined_table(self, table_list):

        """

        Takes a list of astropy tables and returns an empty table with the columns of the old tables (hstacked).

        Used as an auxiliary to generate_full_table()

        """

        combined_table = hstack(table_list)

        combined_table.remove_rows(slice(0, len(combined_table)))

        return(combined_table)



    def merge_tables(self, table1, table2, xmatch_result):

        """

        Returns a new astropy table that consists of table1 and table2 merged horizontally (columns of table1 coming first).



        table1: astropy table

        table2: astropy table

        xmatch_result: the output from an xmatch (using catalogProcessing.xmatch()) of table1 with table2.

        """



        combined_table = self.empty_combined_table([table1, table2])



        assert len(xmatch_result[0]) == len(xmatch_result[1])



        for row_index in range(len(xmatch_result[0])):

            table1_row = [i for i in table1[xmatch_result[0][row_index]]]

            table2_row = [j for j in table2[xmatch_result[1][row_index]]]



            combined_row = table1_row + table2_row

            combined_table.add_row(combined_row, mask=[type(i) is np.ma.core.MaskedConstant for i in combined_row])



        return(combined_table)



    def table_difference(self, table1, table2, id_1, id_2, not_in_table1, gaia_cat, tmass_cat, allwise_cat):

        """

        Returns an astropy table (with the dimenstions of table2) that consists of the objects in table 1 but not in table 2, where the surveys composing table1 is fully contained in table2.

        table1: astropy table

        table2: astropy table

        not_in_table1: list of names of the surveys (gaia, 2mass, allwise) not in table1 but in table2. The purpose of this parameter is to adjust the

                      mask on the full_table returned.

        id_1, id_2: The name by which to identify similar objects (such as gaia designation)

        """



        diff_table = table2



        gaia_mask_shape = [-1 for i in range(len(gaia_cat.colnames))]

        tmass_mask_shape = [-1 for i in range(len(tmass_cat.colnames))]

        allwise_mask_shape = [-1 for i in range(len(allwise_cat.colnames))]



        for row in table1:

            if row[id_1] in table2[id_2]:

                pass

            else:

                row_as_list = [i for i in row]



                full_row = []



                if "gaia" in not_in_table1:

                    full_row = full_row + gaia_mask_shape

                else:

                    full_row = full_row + row_as_list



                if "2mass" in not_in_table1:

                    full_row = full_row + tmass_mask_shape

                elif not "gaia" in not_in_table1:

                    pass

                else:

                    full_row = full_row + row_as_list



                if "allwise" in not_in_table1:

                    full_row = full_row + allwise_mask_shape

                elif not "2mass" in not_in_table1:

                    pass

                else:

                    full_row = full_row + row_as_list



                mask = [i is -1 for i in full_row]

                diff_table.add_row(full_row, mask=mask)



        return(diff_table)







    def generate_full_table(self, gen_small_table = False):

        """

        WORK IN PROGRESS:



        Plan:

        Gaia X tmass -> res1



        res1 X wise -> res2



        append res2



        append res1\res2



        append gaia not xmatched at all



        tmass X wise -> res3

            Here should be no systems in gaia - only tmass\res1 and wise\res2! Check!



        append res3



        append tmass unmatched <<<



        append wise unmatched



        """





        """

        Small list:

        Gaia

            designation, ra, dec, plx, pmra, pmdec, phot_g_mean_mag, bp_rp, bp_g, g_rp

        2MASS

            designation, ra, dec, j_m, h_m, k_m

        AllWISE

            designation, ra, dec, w1mpro, w2mpro, w3mpro, w4mpro

        """

        print("\nGenerating full table...")

        print("Extracting specified data columns...")

        if gen_small_table:

            gaia_cat = self.gaia["designation", "ra", "dec", "ra_error", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag", "bp_rp", "bp_g", "g_rp", "phot_g_n_obs", "phot_g_mean_flux_over_error", "phot_g_mean_flux_error", "phot_g_mean_flux", "phot_bp_mean_mag", "phot_rp_mean_mag"]

            tmass_cat = self.tmass["designation", "ra", "dec","err_maj", "err_min", "j_m", "h_m", "k_m", "j_msig", "h_msig", "k_msig", "k_snr"]

            allwise_cat = self.allwise["designation", "ra", "dec", "sigra", "sigdec", "w1mpro", "w2mpro", "w3mpro", "w4mpro", "w1sigmpro", "w2sigmpro", "w3sigmpro", "w4sigmpro", "w4snr"]

        else:

            gaia_cat = self.gaia

            tmass_cat = self.tmass

            allwise_cat = self.allwise

        print("Crossmatching Gaia with 2MASS...")

        res1 = self.xmatch(gaia_cat, tmass_cat, ["ra", "dec", "ra_error", "dec_error"], ["ra", "dec", "err_maj", "err_min"])



        gaia_X_tmass_table = self.merge_tables(gaia_cat, tmass_cat, res1)

        print("Crossmatching with AllWISE...")

        res2 = self.xmatch(gaia_X_tmass_table, allwise_cat, ["ra_1", "dec_1", "ra_error", "dec_error"], ["ra", "dec", "sigra", "sigdec"])



        gaia_X_tmass_X_allwise = self.merge_tables(gaia_X_tmass_table, allwise_cat, res2)



        full_table = gaia_X_tmass_X_allwise

        print("Adding table differences...")

        print("(Gaia x TMASS) vs (Gaia x TMASS x AllWISE)")

        diff1 = self.table_difference(gaia_X_tmass_table, gaia_X_tmass_X_allwise, "designation_1", "designation", ["allwise"], gaia_cat, tmass_cat, allwise_cat)

        full_table = vstack(full_table, diff1)



        diff2 = self.table_difference(gaia_cat, full_table, "designation", "designation_1", ["2mass", "allwise"], gaia_cat, tmass_cat, allwise_cat)



        full_table = vstack(full_table, diff2)



        res3 = self.xmatch(tmass_cat, allwise_cat, ["ra", "dec", "err_maj", "err_min"], ["ra", "dec", "sigra", "sigdec"])



        tmass_X_allwise = self.merge_tables(tmass_cat, allwise_cat, res3)



        diff3 = self.table_difference(tmass_X_allwise, full_table, "designation_1", "designation_2", ["gaia"], gaia_cat, tmass_cat, allwise_cat)



        full_table = vstack(full_table, diff3)



        diff4 = self.table_difference(tmass_cat, full_table, "designation", "designation_2", ["gaia", "allwise"], gaia_cat, tmass_cat, allwise_cat)



        full_table = vstack(full_table, diff4)



        diff5 = self.table_difference(allwise_cat, full_table, "designation", "designation", ["gaia", "2mass"], gaia_cat, tmass_cat, allwise_cat)



        full_table  = vstack(full_table, diff5)

        print("Renaming columns...")

        full_table.rename_column("designation_1", "gaia_designation")

        full_table.rename_column("designation_2", "2mass_designation")

        full_table.rename_column("designation", "allwise_designation")



        full_table.rename_column("ra_1", "gaia_ra")

        full_table.rename_column("ra_2", "2mass_ra")

        full_table.rename_column("ra", "allwise_ra")



        full_table.rename_column("dec_1", "gaia_dec")

        full_table.rename_column("dec_2", "2mass_dec")

        full_table.rename_column("dec", "allwise_dec")

        v = np.sqrt(full_table['phot_g_n_obs']/full_table['phot_g_mean_flux_over_error'])
        v.name="variability"
        v.unit="(n_obs / mag)^0.5"
        full_table.add_column(v)

        d = 1000 / full_table['parallax']
        d.name = "radial_distance"
        d.unit="pc"
        full_table.add_column(d)


        # summary statistics for data
        print("Full table generated. Summary:")
        print("\nGaia:")
        print("Total rows: "+ str(len(full_table['gaia_designation'])))
        print("G: " + str(len(full_table['phot_g_mean_mag'].mask.nonzero()[0])))
        print("BP: " + str(len(full_table['phot_bp_mean_mag'].mask.nonzero()[0])))
        print("RP: " + str(len(full_table['phot_rp_mean_mag'].mask.nonzero()[0])))
        print("PLX: " + str(len(full_table['parallax'].mask.nonzero()[0])))
        print("PMRA: " + str(len(full_table['pmra'].mask.nonzero()[0])))
        print("PMDEC: " + str(len(full_table['pmdec'].mask.nonzero()[0])))

        print("\n2MASS:")
        print("Total rows: " +str(len(full_table['2mass_designation'])))
        print("J: " + str(len(full_table['j_m'].mask.nonzero()[0])))
        print("H: " + str(len(full_table['h_m'].mask.nonzero()[0])))
        print("K: " + str(len(full_table['k_m'].mask.nonzero()[0])))

        print("\nallwise:")
        print("Total rows: " + str(len(full_table['allwise_designation'])))
        print("W1: " + str(len(full_table['w1mpro'].mask.nonzero()[0])))
        print("W2: " + str(len(full_table['w2mpro'].mask.nonzero()[0])))
        print("W3: " + str(len(full_table['w3mpro'].mask.nonzero()[0])))
        print("W4: " + str(len(full_table['w4mpro'].mask.nonzero()[0])))

        return(full_table)



    def extract_ft_part(self, cat_itr, full_table):

        """

        returns a copy of the full table without the parts that doesn't have data from the surveys in caty_itr.

        """





        select_dict = {"gaia": "gaia_designation",

                       "2mass": "2mass_designation",

                       "allwise": "allwise_designation"}

        # Creates an empty copy of full_table

        reduced_table = full_table.copy()

        reduced_table.remove_rows(slice(0, len(reduced_table)))





        for row_ind in range(len(full_table)):

            add_row = []

            for i in cat_itr:

                if not(np.ma.is_masked(full_table[row_ind][select_dict[i]])):

                    add_row.append(True)

                else:

                    add_row.append(False)



            if all(add_row):

                mask = []

                for val in full_table[row_ind]:

                    mask.append(np.ma.is_masked(val))



                reduced_table.add_row(vals=full_table[row_ind],mask=mask)



        return(reduced_table)

    def identify_clusters(self, table, columns, expected_clusters=None):

        # can include:
        # min_samples
        # min_cluster_size
        # cluster_selection_epsilon
        # allow_single_cluster
        # alpha

        #data = []
        #candidates_table = table[:]

        print("Organizing clustering data...")

        print(table.colnames)

        '''cut_string = ""
        for c in columns:
            cut_string += "self.full_table['" + c + "'] > -1000 and "
        cut_string = cut_string.strip().strip('and').strip()

        print(cut_string)

        phot_sp = Photometry(table)
        candidates_table, _ = phot_sp.apply_cut(cut_string)

        data = Table()
        for c in columns:
            data.add_column(candidates_table[c])
        '''

        '''for c in columns:
            data_c = []

            for i in range(len(table[c])):
                if candidates_table[i] != None and table[c][i] != None:
                    data_c.append(table[c][i])
                else:
                    candidates_table[i] = None

            data = np.vstack((data, data_c))

        data = data.t
        '''

        '''
        for i in range(len(table)):
            include_flag = True
            for c in columns:
                if table[c][i] == None:
                    include_flag = False
            if include_flag:
                candidates_table.add_row(table[i])
        
        candidates_table.colnames = table.colnames

        data = Table()
        for c in columns:
            data.add_column(candidates_table[c])
        '''
        '''column_A = columns[0]
        column_B = columns[1]
        
        data_A = []
        data_B = []
        candidates_table = table[:]
        count = 0
        count_included = 0
        for i in range(len(table[columns[0]])):
            if table[column_A][i] != None and table[column_B][i] != None:
                count_included += 1
                data_A.append(table[column_A][i])
                data_B.append(table[column_B][i])
            else:
                candidates_table.remove_row(i-count)
                count += 1

        print(count_included)

        data = np.vstack((data_A, data_B)).T'''
        

        data = []
        candidates_table = table[:]
        count = 0
        print(columns)
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

        print(len(data))
        print(len(data[0]))

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





class Astrometry:

    """

    This class defines methods for astrometrical analysis.

    Instantiated with a full_table.

    """

    def __init__(self, full_table):

        self.full_table = full_table





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

    def gaussian(self, x, a, b, c):

        return(a*np.exp(-(x-b)**2/(2*c**2)))





    def remove_masked(self, ap_column, header):

        """

        Takes an astroquery result and a column header (str) and returns a list not containing the masked elements.

        """

        data_list = []

        for i in ap_column[header]:

            if type(i)== np.float64:

                data_list.append(i)



        return(data_list)





    def make_hist(self, data_list):

        """

        Makes a histogram using np.histogram. Return a 2-tuple consisting of the output

        of the output of np.histogram and the midpoints of the bins.

        """



        #print(data_list.data)



        data_list = data_list[~data_list.mask]



        histogram = np.histogram(data_list.data, bins="auto")



        midpoints_of_bins = self.get_bin_mid(histogram[1])

        return((histogram, midpoints_of_bins))





    def fit_gaussian(self, colname, full_column=False):

        """

        Fits a gaussian to a list of data points by automatically sorting them into

        bins, and fitting a gaussian using scipy's curve_fit.



        Input: List of data points



        Output: 3-tuple (a,b,c) consisting of the coefficients of the

                gaussian: a*e^{-(x-b)^2/(2c^2)}.



        Dependency: plt, numpy, scipy

        """

        #number_of_bins = bin_number(data_list)



        #hist_values = bin_data(data_list, number_of_bins)





        if full_column:

            histogram = self.make_hist(colname)

        else:

            histogram = self.make_hist(self.full_table[colname])

        fit = scipy.optimize.curve_fit(self.gaussian, histogram[1], histogram[0][0])



        return((fit[0][0], fit[0][1], fit[0][2]))









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





class Visualization:

    """

    This provides functionality for plotting and visualization.



    Methods in this class:



    > To plot any two columns of the full_table against each other



    > To plot a histogram and gaaussian of any column of the full_table



    > Overlay coordinates from all three surveys (gaia, 2mass, allwise) over images of the sky (DSS, AllWISE, 2MASS, etc.)



    > Plot any cut equation as a line in the diagrams.



    > Define certain <standard plots> that are allways displayed when the main is run with new HAeBe coordinates.



    """



    def __init__(self, full_table, ra, dec):

        self.full_table = full_table

        self.ra = ra

        self.dec = dec

        self.phot = Photometry(full_table)





    def format_coords(self, coord_list):

        #TODO: Make this extract the coordinates from the full_list and use that

        return(SkyCoord(coord_list, unit='deg'))





    def get_surveys(self, survey_list, field_side_length):

        """

        Center: skycoord coordinate+

        survey_list: list or iterable consisting of strings in the form of SkyView surveynames.

        Full description of all the surveys availible can be found here:

        https://skyview.gsfc.nasa.gov/current/cgi/survey.plhttps://skyview.gsfc.nasa.gov/current/cgi/survey.pl



        Some surveynames are:

        "2MASS-J"

        "2MASS-H"

        "2MASS-K"

        "WISE 3.4"

        "WISE 4.6"

        "WISE 12"

        "WISE 22"

        "DSS"

        "SDSSg"

        "SDSSi"

        "SDSSr"

        "SDSSu"

        "SDSSz"











        Returns a list of astropy.fits.HDUList objects. All elements of this

        list has an attribute "".data" that can be passed to plt.imshow().

        """



        print("fetching surveys:")
        print(survey_list)
        print("with radius:")
        print(field_side_length)

        center=self.ra + self.dec

        return(SkyView.get_images(position=center, survey=survey_list, radius=(field_side_length)))





    def plot_surveys(self, img_list, formatted_coord_list, survey_list, plot_coords=True):



        #print(img_list)

        wcs_list = []

        for i in img_list:

            #print(i[0])

            wcs_list.append(wcs.WCS(i[0].header))



        fig = plt.figure(figsize=(20,20))

        fig.clf()

        for i in range(len(img_list)):

            ind = i+1

            #ax = fig.add_subplot(len(img_list), 3, ind, projection=wcs_list[0])
            ax = fig.add_subplot(len(img_list), 2, ind, projection=wcs_list[0])

            ax.set_xlabel('RA')

            ax.set_ylabel('Dec')

            ax.set_title(survey_list[i])

            #ax.imshow(img_list[i][0].data, cmap='gray')
            norm = ImageNormalize(img_list[i][0].data, interval=ZScaleInterval(), stretch=LinearStretch())
            ax.imshow(img_list[i][0].data, origin='lower', norm=norm, cmap='gray_r')


            #print((formattedcoord_list[i].ra, formattedcoord_list[i].dec))

            if plot_coords:

                for cat in formatted_coord_list:

                    ax.plot(cat.ra, cat.dec, 'o', transform=ax.get_transform('icrs'), mec='r', mfc='none')





    def plot_images(self, survey_list, cats_to_plot=None):

        """

        This method plotts images of the sky based on the surveys in survey_list.

        If coord_overlay is True, plots sources from self.full_table on top.

        If cats_to_plot != None, must be a list consisting of "gaia", "2mass", "allwise"

        Will overlay only sources with information from catalogs in cats_to_plot.



        survey_list: list/iterable of strings representing the surveys to pull images from

        """



        #formatted_center = self.format_coords(self.HAeBe_coord)
        formatted_center = SkyCoord(self.ra, self.dec)




        # coord_list = []

        # for row in self.full_table["ra", "dec"]:

        #     coord_list.append((row[0], row[1]))



        try:

            coord_list = []

            for i in range(len(cats_to_plot)):

                coord_list.append([])

                if cats_to_plot[i] == "gaia":
                    count_gaia = 0
                    for row in self.full_table["gaia_ra", "gaia_dec"]:
                        count_gaia += 1
                        if row[0] != 0:
                                coord_list[i].append((row[0], row[1]))
                    print("Added gaia coords!")
                    #print(str(count_gaia))

                if cats_to_plot[i] == "2mass":
                    count_2mass = 0
                    for row in self.full_table["2mass_ra", "2mass_dec"]:
                        count_2mass += 1
                        if row[0] != 0:
                                coord_list[i].append((row[0], row[1]))
                    print("Added 2mass coords!")

                if cats_to_plot[i] == "allwise":
                    count_allwise = 0
                    for row in self.full_table["allwise_ra", "allwise_dec"]:
                        count_allwise += 1
                        if row[0] != 0:
                                coord_list[i].append((row[0], row[1]))
                    print("Added allwise coords!")
                    #print(str(count_allwise))



            ra_max_vals = []

            ra_min_vals = []

            dec_max_vals = []

            dec_min_vals = []

            for cat in coord_list:

                ra_max_vals.append(max([pair[0] for pair in cat]))

                ra_min_vals.append(max([pair[0] for pair in cat]))

                dec_max_vals.append(max([pair[1] for pair in cat]))

                dec_min_vals.append(max([pair[1] for pair in cat]))



            ra_min = min(ra_min_vals)

            ra_max = max(ra_max_vals)

            dec_min = min(dec_min_vals)

            dec_max = max(dec_max_vals)



            field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

            ra_min = min(formatted_coord_list[0].ra)

            ra_max = max(formatted_coord_list[0].ra)

            dec_min = min(formatted_coord_list[0].dec)

            dec_max = max(formatted_coord_list[0].dec)

            #print(ra_min)
            #print(ra_max)
            #print(ra_max - ra_min)



            field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            #formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

            img_list = self.get_surveys(survey_list, field_side_length)

            self.plot_surveys(img_list, formatted_coord_list, survey_list, plot_coords=True)


        except TypeError:

            print("TypeError occurred!")

            '''ra_min = min(self.full_table["gaia_ra"])

            ra_max = max(self.full_table["gaia_ra"])

            dec_min = min(self.full_table["gaia_dec"])

            dec_max = max(self.full_table["gaia_dec"])

            #print(ra_min, ra_max)

            field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            #print(field_side_length)

            img_list = self.get_surveys(survey_list, field_side_length*u.deg)

            self.plot_surveys(img_list, [], survey_list, plot_coords=False)'''



        '''#formatted_coord_list = self.format_coords(coord_list)
        formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

        #print(type(formatted_coord_list))
        #print(formatted_coord_list[0])

        ra_min = min(formatted_coord_list[0].ra)

        ra_max = max(formatted_coord_list[0].ra)

        dec_min = min(formatted_coord_list[0].dec)

        dec_max = max(formatted_coord_list[0].dec)

        #print(ra_min)
        #print(ra_max)
        #print(ra_max - ra_min)



        field_side_length = max((ra_max - ra_min, dec_max - dec_min))

        #print(field_side_length)


        img_list = self.get_surveys(survey_list, field_side_length)

        print("Number of surveys:")
        print(len(survey_list))
        print("Number of images:")
        print(len(img_list))

        self.plot_surveys(img_list, formatted_coord_list, survey_list)'''



    def plot_hist(self, colname, xlable, cut=None):

        plt.figure()

        if cut != None:



            phot = Photometry(self.full_table)



            column = phot.apply_cut(cut)[0][colname]

            #print("column:")

            #column.pprint(max_lines=-1)

            astrometry = Astrometry(self.full_table)

            try:

                fit = astrometry.fit_gaussian(column, full_column=True)

            #except RuntimeError as rte:
            except:

                hist = astrometry.make_hist(column)



                x = np.linspace(min(column),max(column),1000)



                plt.xlabel(xlable,fontsize=14)

                plt.ylabel('# count',fontsize=14)



                plt.hist(column, bins=hist[0][1])

            else:

                hist = astrometry.make_hist(column)



                x = np.linspace(min(column),max(column),1000)







                plt.xlabel(xlable,fontsize=14)

                plt.ylabel('# count',fontsize=14)



                plt.hist(column, bins=hist[0][1])

                plt.plot(x, astrometry.gaussian(x, fit[0], fit[1], fit[2]))

            finally:

                plt.show()



        else:

            astrometry = Astrometry(self.full_table)



            fit = astrometry.fit_gaussian(colname)

            #self.full_table[colname].pprint(max_lines=-1)



            hist = astrometry.make_hist(self.full_table[colname])



            #print(fit)



            x = np.linspace(min(self.full_table[colname]),max(self.full_table[colname]),1000)



            #plt.title('',fontsize=16)

            plt.xlabel(xlable,fontsize=14)

            plt.ylabel('# count',fontsize=14)



            plt.hist(self.full_table[colname], bins=hist[0][1])

            plt.plot(x, astrometry.gaussian(x, fit[0], fit[1], fit[2]))

            plt.show()


    def plot_error(self, parameter):

        try:
            good_data = self.full_table["gaia_ra"] > 0
            gmag = self.full_table['phot_g_mean_mag']
            param_data = self.full_table[parameter]
            param_error = self.full_table[parameter + "_error"]
        except:
            print("Given parameter does not appear in the imported data. Check parameter name for consistency with full_table columns.")

        fractional = np.abs(100.* param_error / param_data)

        matplotlib.rcParams['figure.figsize'] = (18,6)
        plt.subplots(nrows=1, ncols=3, sharey=True)
        plt.subplot(1,3,1)
        plt.title('cluster')
        plt.gca().minorticks_on()
        plt.errorbar(gmag[good_data],param_data[good_data],yerr=param_error[good_data], fmt='.',ms=2)
        plt.scatter(gmag[good_data],param_data[good_data])
        plt.xlabel('G [mag]')
        plt.ylabel(parameter + ' ' + str(param_data.unit))
        plt.subplot(1,3,2)
        plt.scatter(param_error[good_data],param_data[good_data])
        plt.xlabel(parameter + '_error ' + str(param_data.unit))
        plt.ylabel(parameter + ' ' + str(param_data.unit))
        plt.subplot(1,3,3)
        plt.scatter(np.log10(fractional[good_data]),param_data[good_data])
        plt.axvline(np.log10(1), linestyle="solid",color="black")
        plt.axvline(np.log10(10), linestyle="dotted",color="blue")
        plt.axvline(np.log10(20), linestyle="dotted",color="blue")
        plt.axvline(np.log10(30), linestyle="dotted",color="blue")
        plt.xlabel('log fractional error [%]')
        plt.ylabel(parameter + ' ' + str(param_data.unit))
        plt.show()


    def plot_list(self, colname_list):

        #print(self.full_table[colname1], self.full_table[colname1])



        fig = plt.figure(figsize=(20,20))

        fig.clf()



        #fig, axs = plt.subplots((len(colname_list)))



        for i in range(len(colname_list)):

            ind = i+1

            ax = fig.add_subplot(len(colname_list), 2, ind)



            ax.set_xlabel(colname_list[i][0] + " [{}]".format(self.full_table[colname_list[i][0]].unit), fontsize=14)

            ax.set_ylabel(colname_list[i][1] + " [{}]".format(self.full_table[colname_list[i][1]].unit), fontsize=14)

            #print("hi")

            ax.scatter(self.full_table[colname_list[i][0]], self.full_table[colname_list[i][1]], s=8, marker="o")

        plt.show()



    def plot(self, colname1, colname2, xlim=None, ylim=None, squared=False, invert_y=False, invert_x=False):

        """

        plots colname1 on the x-axis and colname2 on the y-axis.

        colname1/2: tuples on two strings, (colname, label)



        """

        plt.figure()

        if type(colname1) is str:

            x = self.full_table[colname1]

            plt.xlabel(colname1 + " [{}]".format(self.full_table[colname1].unit), fontsize=14)

        elif type(colname1[0]) is str:

            #print("check")

            x = self.full_table[colname1[0]]

            plt.xlabel(colname1[1] + " [{}]".format(self.full_table[colname1[0]].unit), fontsize=14)

        else:

            x = colname1[0]

            plt.xlabel(colname1[1], fontsize=14)



        if type(colname2) is str:

            y = self.full_table[colname2]

            plt.ylabel(colname2 + " [{}]".format(self.full_table[colname2].unit), fontsize=14)

        elif type(colname2[0]) is str:

            y = self.full_table[colname2[0]]

            plt.ylabel(colname2[1] + " [{}]".format(self.full_table[colname2[0]].unit), fontsize=14)

        else:

            y = colname2[0]

            plt.ylabel(colname2[1], fontsize=14)



        if xlim != None:

            plt.xlim(xlim)

        if ylim != None:

            plt.ylim(ylim)





        if invert_y:

            plt.gca().invert_yaxis()

        if invert_x:

            plt.gca().invert_xaxis()

        #print(type(x), type(y))



        plt.scatter(x,y, s=8, marker="o", c="k")

        if squared:

            plt.gcf().set_size_inches(6,6)

        plt.show()





    def double_plot(self, colname_list, lims=None, squared=False):

        """

        colanme_list:

        [(colname11, colname12), (colname21, colname22)]

        [(colname11, colname12, xlabel1, ylabel1),(colname21, colname22, xlabel2, ylabel2)]

        or

        [(column11, column12, xlabel1, ylabel1),(col21, col22, xlabel2, ylabel2)]



        lims:

        [(xlim1, ylim1), (xlim1, ylim2)]

        """



        fig = plt.figure(figsize=(10,10))

        fig.clf()



        #if squared:

        #        plt.gca().set_aspect('equal', adjustable='box')





        for i, col_tup in enumerate(colname_list):



            ind = i+1

            if squared:

                ax = fig.add_subplot(1, 2, ind, aspect="equal")

            else:

                ax = fig.add_subplot(1, 2, ind)





            if type(col_tup[0]) is str:

                x = self.full_table[col_tup[0]]

                ax.set_xlabel(col_tup[0] + " [{}]".format(self.full_table[col_tup[0]].unit), fontsize=14)

            #elif type(col_tup[2]) is str:

            #    #print("check")

            #    x = self.full_table[col_tup[0]]

            #    ax.set_xlabel(col_tup[2] + " [{}]".format(self.full_table[col_tup[0]].unit), fontsize=14)

            else:

                x = col_tup[0]

                ax.set_xlabel(col_tup[2], fontsize=14)



            if type(col_tup[1]) is str:

                y = self.full_table[col_tup[1]]

                plt.ylabel(col_tup[1] + " [{}]".format(self.full_table[col_tup[1]].unit), fontsize=14)

            #elif type(col_tup[3]) is str:

            #    y = self.full_table[col_tup[1]]

            #    plt.ylabel(colname2[1] + " [{}]".format(self.full_table[colname2[0]].unit), fontsize=14)

            else:

                y = col_tup[1]

                plt.ylabel(col_tup[3], fontsize=14)



            #if lims != None:

            #    ax.set_xlim(lims[i][0])

            #    ax.set_ylim(lims[i][1]

            #print(type(x), type(y))



            ax.scatter(x,y, s=8, marker="o")

        plt.show()



    def plot_double_zoom(self, col1, col2):
        """
        Not yet implemented
        """

        pass





    def plot_removed(self, plot_list, xlabel, ylabel, xlim=None, ylim=None, invert_x=False, invert_y=False, squared=False):

        """
        Plots each list in the tuples in plot_list against each other in different colors. Plots up to six different lists of tuples.

        Arguments:
            plot_list [list of tuples]: list of tuples of two lists that represents x and y of different datasets to be plotted together.
            xlabel [string]: string to label the x axis of the plot
            ylabel [string]: string to label the y axis of the plot
            xlim [tuple]: (lower x limit, upper x limit) for plot. Default: None
            ylim [tuple]: (lower y limit, upper y limit) for plot. Default: None
            invert_x [boolean]: whether or not to invert the x axis. Default: False
            invert_y [boolean]: whether or not to invert the y axis. Default: False
            squared [boolean]: whether or not to force generation of a square plot. Default: False

        Returns:
            [none]
        """

        # create pyplot object
        fig = plt.figure()
        fig.clf()

        # create subfigure to plot lists
        ax = fig.add_subplot(1,1,1)

        # apply axis labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # create ordered list of colors
        colors = ["k", "r", "b", "g", "m", "y"]

        # scatter plot lists on the subfigure object
        for i, tup in enumerate(plot_list):
            ax.scatter(tup[0], tup[1], c=colors[i], s=8, marker="o")

        # set optional parameters
        if xlim:
            ax.set_xlim(xlim)

        if ylim:
            ax.set_ylim(ylim)

        if invert_y:
            plt.gca().invert_yaxis()

        if invert_x:
            plt.gca().invert_xaxis()

        if squared:
            plt.gcf().set_size_inches(6,6)

        # display plot
        plt.show()


    def cut_and_plot(self, cut, col1, col2, phot=None, xlim=None, ylim=None, squared=False, invert_y=False, invert_x=False):

        """
        Applies a cut to the full_table and then plots the sources in different colors based on if they meet the cut or not.

        col1/2: (colname, label)

        phot: instatiation of Photometry (((I think this is ugly, TODO: Fix)))
        """

        (cut_true, cut_false) = self.phot.apply_cut(cut)
        
        xlabel = str(col1[1]) + " [{}]".format(self.full_table[col1[0]].unit)
        ylabel = str(col2[1]) + " [{}]".format(self.full_table[col2[0]].unit)

        self.plot_removed([(cut_true[col1[0]], cut_true[col2[0]]), (cut_false[col1[0]], cut_false[col2[0]])],
                xlabel, ylabel, xlim, ylim, squared=squared, invert_y = invert_y, invert_x = invert_x)


    def plot_diagnostics(self):
        """
        generates hard-coded plots to display useful diagnostic information
        """

        
        
        # define variables
        M_g = [g + 5 - 5*np.log10(1000/p) for g,p in zip(self.full_table["phot_g_mean_mag"], self.full_table["parallax"])]
        g_k = self.phot.subtract_cols("phot_g_mean_mag", "k_m")
        j_k = self.phot.subtract_cols("j_m", "k_m")
        w1_w2 = self.phot.subtract_cols("w1mpro", "w2mpro")
        j_h = self.phot.subtract_cols("j_m", "h_m")
        h_k = self.phot.subtract_cols("h_m", "k_m")
        g_h = self.phot.subtract_cols("phot_g_mean_mag", "h_m")
        k_w1 = self.phot.subtract_cols("k_m", "w1mpro")
        k_w2 = self.phot.subtract_cols("k_m", "w2mpro")
        k_w3 = self.phot.subtract_cols("k_m", "w3mpro")

        try:
            cut_true = self.full_table[pd.eval("""(-1.81 -2*1.31 < self.phot.full_table["pmra"]) & (self.phot.full_table["pmra"] < -1.81 +2*1.31) & (-2.67 -2*1.44 < self.phot.full_table["pmdec"]) & (self.phot.full_table["pmdec"] < -2.67 +2*1.44)""")]
            cut_false = self.full_table[~pd.eval("""(-1.81 -2*1.31 < self.phot.full_table["pmra"]) & (self.phot.full_table["pmra"] < -1.81 +2*1.31) & (-2.67 -2*1.44 < self.phot.full_table["pmdec"]) & (self.phot.full_table["pmdec"] < -2.67 +2*1.44)""")]

            cut_mg_true = [g + 5 - 5*np.log10(1000/p) for g,p in zip(cut_true["phot_g_mean_mag"], cut_true["parallax"])]
            cut_bprp_true = cut_true["bp_rp"]
            cut_mg_false = [g + 5 - 5*np.log10(1000/p) for g,p in zip(cut_false["phot_g_mean_mag"], cut_false["parallax"])]
            cut_bprp_false = cut_false["bp_rp"]

            cut_1s = """(self.full_table["pmra"] < -1.81 +1.31) & (self.full_table["pmra"] > -1.81 -1.31) & (self.full_table["pmdec"] < -2.67 +1.44) & (self.full_table["pmdec"] > -2.67 -1.44)"""

            cut_2s = """(self.full_table["pmra"] < -1.81 +2*1.31) & (self.full_table["pmra"] > -1.81 -2*1.31) & (self.full_table["pmdec"] < -2.67 +2*1.44) & (self.full_table["pmdec"] > -2.67 -2*1.44)"""

            cut_outliers = """(self.full_table["pmra"] < 15) & (self.full_table["pmra"] > -15) & (self.full_table["pmdec"] < 15) & (self.full_table["pmdec"] > -15)"""

            cut_plx_1s = """(self.full_table["parallax"]< 1.04+0.23) & (self.full_table["parallax"]> 1.04-0.23)"""



            # generate and display plots with documentation

            print("2mass RA vs 2mass Dec. J-H > 0.7 shown in black, J-H < 0.7 shown in red.")
            self.cut_and_plot("J-H > 0.7", ("2mass_ra", "2mass RA"), ("2mass_dec", "2mass Dec"), squared=True, invert_x=True)

            print("Gaia RA vs Gaia Dec.")
            self.plot(("gaia_ra", "Gaia RA"), ("gaia_dec", "Gaia Dec"), squared=True, invert_x=True )

            print("Gaia PM RA vs Gaia PM Dec.")
            self.plot(("pmra", "pm RA (Gaia)"), ("pmdec", "pm Dec (Gaia)"), xlim=(-30,30), ylim=(-30,30), squared=True)

            print("Gaia PM RA vs Gaia PM Dec (closer detail).")
            self.plot(("pmra", "pm RA (Gaia)"), ("pmdec", "pm Dec (Gaia)"), xlim=(-10,10), ylim=(-10,10), squared=True)

            print("Parallax vs Gaia Dec.")
            self.plot(("parallax", "Parallax"), ("gaia_dec", "Dec (Gaia)"), xlim=(0,5))

            #print("what does double_plot do?")
            #self.double_plot([("pmra", "pmdec"), ("pmra", "pmdec")], [(-20,20), (-10, 10)])

            print("BP/RP vs G Mean Magnitude.")
            self.plot("bp_rp", "phot_g_mean_mag")

            print("BP/RP vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot("bp_rp", (M_g, "$M_G [mag]$"))

            print("G-K Magnitude vs G Mean Magnitude.")
            self.plot((g_k, "G-K [mag]"), "phot_g_mean_mag")

            print("G-K Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot((g_k, "G-K [mag]"), (M_g, "$M_G [mag]$"))

            print("J-K Magnitude vs J-M Magnitude.")
            self.plot((j_k, "J-K [mag]"), "j_m")

            print("W1-W2 Magnitude vs W1 Magnitude.")
            self.plot((w1_w2, "W1-W2 [mag]"), "w1mpro")

            print("J-H Magnitude vs H-K Magnitude.")
            self.plot((j_h, "J-H [mag]"), (h_k, "H-K [mag]"))

            print("G-H Magnitude vs H-K Magnitude.")
            self.plot((g_h, "G-H [mag]"), (h_k, "H-K [mag]"))

            print("K-W2 Magnitude vs G-K Magnitude.")
            self.plot((k_w2, "K - W2 [mag]"), (g_k, "G - K [mag]")) 

            print("BP/RP Magnitude vs G Mean Magnitude.")
            self.plot("bp_rp", "phot_g_mean_mag", invert_y=True)

            print("2mass RA vs 2mass Dec. K-W1 > 0.2 shown in black, K-w1 < 0.2 shown in red.")
            self.cut_and_plot("K-W1 > 0.2", ("2mass_ra", "RA"), ("2mass_dec", "Dec"), squared=True)

            print("PMRA vs PMDec. -3.12 < PMRA < -0.5 AND -4.17 < PMDec < -1.23 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_1s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-10,10), ylim=(-10,10), squared=True)

            print("PMRA vs PMDec - closer detail. -3.12 < PMRA < -0.5 AND -4.17 < PMDec < -1.23 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_1s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-6,1), ylim=(-6,1), squared=True)

            print("PMRA vs PMDec. -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_2s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-10,10), ylim=(-10,10), squared=True)

            print("PMRA vs PMDec - closer detail. -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_2s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-6,1), ylim=(-6,1), squared=True)

            print("Parallax vs PMRA. 0.81 < Parallax < 1.27 shown in black, otherwise shown in red.")
            self.cut_and_plot(cut_plx_1s, ("parallax", "Parallax"),("pmra", "pm RA"), xlim=(0,5), squared=True)

            print("PMRA histogram. Data shown satisfy -15 < PMRA < 15 AND -15 < PMDec < 15. Units in mas/yr.")
            self.plot_hist("pmra", "pm RA", cut=cut_outliers)

            print("PMDec histogram. Data shown satisfy -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21. Units in mas/yr.")
            self.plot_hist("pmdec", "pm Dec", cut=cut_2s)

            print("BP-RP Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot_removed([(cut_bprp_true, cut_mg_true), (cut_bprp_false, cut_mg_false)], "BP - RP", "$M_G$", invert_y=True)

            print("BP-RP Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax )). More detail.")
            self.plot_removed([(cut_bprp_true, cut_mg_true), (cut_bprp_false, cut_mg_false)], "BP - RP", "$M_G$", xlim=(0.1,3.5), ylim=(-1, 15), invert_y=True)
        except:
            print("An error occurred while plotting diagnostics. This usually occurs because no sources passed a photometric cut.")

    def plot_clusters(self, table, membership, column_A, column_B, xlim= None, ylim=None, invert_x=False, invert_y=True, squared=True):
        # create pyplot object
        fig = plt.figure()
        fig.clf()

        # create subfigure to plot lists
        ax = fig.add_subplot(1,1,1)

        # apply axis labels
        ax.set_xlabel(str(column_A))
        ax.set_ylabel(str(column_B))

        A_none = []
        B_none = []
        for j in range(len(membership)):
            if membership[j] == -1:
                A_none.append(table[column_A][j])
                B_none.append(table[column_B][j])
        ax.scatter(A_none, B_none)

        for i in range(max(membership)+1):
            print("plotting cluster" + str(i))
            A_i = []
            B_i = []
            for j in range(len(membership)):
                if membership[j] == i:
                    A_i.append(table[column_A][j])
                    B_i.append(table[column_B][j])
            ax.scatter(A_i, B_i)
            name = "cluster_" + str(i) + "_" + column_A + column_B + ".png"
            figsp = plt.figure()
            axsp = figsp.add_subplot(1,1,1)
            axsp.set_xlabel(str(column_A))
            axsp.set_ylabel(str(column_B))
            axsp.scatter(A_i, B_i)
            figsp.savefig(name)
            figsp.clf()

        if xlim:
            ax.set_xlim(xlim)

        if ylim:
            ax.set_ylim(ylim)

        if invert_y:
            plt.gca().invert_yaxis()

        if invert_x:
            plt.gca().invert_xaxis()

        if squared:
            plt.gcf().set_size_inches(6,6)

        plt.show()