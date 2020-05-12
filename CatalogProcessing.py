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

from CatalogTable import *

class catalogProcessing:

    """

    Instantiated with

        > a string representing HAeBe coordinates in degrees.

        > radius: number, either physical or angular radius (depending on if parallax is given).

        > parallax (optional): number, in milliarcseconds. If given, radius given is assumed to be physical, and will

          be recalculated into an angular one.

    """

    def __init__(self, ra, dec, radius, parallax=None, io=None):

        #self.skycoord = SkyCoord(HAeBe_coord, unit=(u.hourangle, u.deg))
        self.skycoord = SkyCoord(ra, dec)

        self.radius = radius

        self.parallax = parallax

        self.io = io


        if self.parallax == None:

            pass

        else:

            self.radius = self.get_radius()


        print("Querying Gaia...")
        self.gaia = self.gaia_query()
        print("Done querying Gaia!\n")
        print("Querying 2MASS...")
        self.tmass = self.ir_query("fp_psc")
        print("Done querying 2MASS!\n")
        print("Querying AllWISE...")
        self.allwise = self.ir_query("allwise_p3as_psd")
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


    def query_catalogs(self):
        """
        Performs queries of catalogs across the specified radius.
        For large queries, splits the query into evenly-sized pieces such that no query is larger than 1 deg

        Maybe implementation requires that we query over squares of sky?
        """

        pass


    def gaia_query(self):

        """
        Performs a query on gaia DR2 using self.skycoord and self.radius.
        Note: "Exception: 500" tends to happen when the coordinates used in search_string is not formatted correctly.

        Arguments:
            [none]

        Returns:
            Astropy table with query result.
        """

        catalog = ["gaia"]

        ra = self.skycoord.ra.degree

        dec = self.skycoord.dec.degree


        #TODO: include code to write diagnostic info to log file incl. ra, dec, radius, search string

        search_string = "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{0},{1},{2}))=1;".format(ra, dec, self.radius)

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
        query_results = job.get_results()

        return(CatalogTable(catalog,query_results))

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

        catalog = [irquery]

        search_string = "SELECT * FROM {} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))=1 ".format(ircat_name, str(self.skycoord.ra.degree), str(self.skycoord.dec.degree), self.radius)

        if view_adql:
            print(search_string)

        ipac = TapPlus(url="https://irsa.ipac.caltech.edu/TAP/")

        print("Creating query...")
        job = ipac.launch_job_async(search_string)


        print("Retrieving results...")
        return(CatalogTable(catalog,job.get_results()))

    # TODO: needs documentation updates
    # TODO: needs better error handling
    def xmatch(self, catalog1, catalog2, rad=-1, show_dynamic_radii=True):

        """
        valid surveynames: "2mass", "allwise", "gaia"


        Crossmatch cat1 with cat2 at a maximum xmatch radius of rad.

        Uses SkyCoord.match_to_catalog_sky() method which implements a nearest neighbour algorithm.


        cat1, cat2: full query results (astropy tables resulting from querying catalogs)
        colnames1, colnames2: list of the column names of ra dec and respective errors in cat1, cat2. Either [ra, dec] or [ra, dec, err_ra, err_dec] for


        Returns: list of lists [cat1_idx, new_idx, new_d2d] or if show_dynamic_radii=True [cat1_idx, new_idx, new_d2d, xmatch_rad]

            cat1_idx: index of the element in cat2 that is nearest to the corresponding element in cat1

                              (cat1_idx[0] = 69 means that the first element in cat1 is matched towards the 69th element in cat2)

        """

        surveys_coord_colnames = {"2mass":["2mass_ra", "2mass_dec", "err_maj", "err_min"],

                                  "allwise":["allwise_ra", "allwise_dec", "sigra", "sigdec"],

                                  "gaia":["gaia_ra", "gaia_dec", "ra_error", "dec_error"]}


        colnames1 = surveys_coord_colnames[catalog1[0][0]]
        cat1 = catalog1[1]

        colnames2 = surveys_coord_colnames[catalog2[0][0]]
        cat2 = catalog2[1]

        ra1 = cat1[colnames1[0]]
        dec1 = cat1[colnames1[1]]

        ra2 = cat2[colnames2[0]]
        dec2 = cat2[colnames2[1]]


        if rad < 0:

            ra_err1 = cat1[colnames1[2]]
            dec_err1 = cat1[colnames1[3]]

            ra_err2 = cat2[colnames2[2]]
            dec_err2 = cat2[colnames2[3]]



        skycoord1 = SkyCoord(ra1, dec1)

        skycoord2 = SkyCoord(ra2, dec2)


        idx, d2d, d3d = skycoord1.match_to_catalog_sky(skycoord2)


        #Imposes constraint on the xmatch radius.
        if rad < 0:

            constraint = "dynamic"

        else:

            constraint = d2d < rad*u.arcsec


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

                # compute max acceptable error radius
                # cosine because we're on a sphere
                error_radius1 = max(ra_err1[j]*np.cos(dec1[j]), dec_err1[j])
                error_radius2 = max(ra_err2[j]*np.cos(dec2[j]), dec_err2[j])

                # require that the sources are closer to each other than the sum of their errors
                match_constraint = error_radius1 + error_radius2

                # TODO: what is this for?
                if show_dynamic_radii is True:

                    xmatch_rad.append(match_constraint)

                else:

                    pass

                # indpair[1] is 2d distance (projected onto the sky) between the source in catalog 1
                # and the matched source in catalog 2
                dist = indpair[1].arcsec

                #const = coordinates.Angle(match_constraint*u.arcsec)

                # TODO: why do we know match_constraint is in arcsec?
                if dist > match_constraint:
                    # if this is true then remove the crossmatch entry

                    new_idx = np.delete(new_idx, j)
                    new_d2d = np.delete(new_d2d, j)

                    # subtract 1 from the index under consideration because we just deleted the entry
                    # we were looking at, and we'll increment j at the end of the loop - this ensures
                    # we stay at the same index.
                    j-=1

                else:
                    
                    cat1_idx.append(i)

                j+=1
                i+=1


        else:

            j = 0 # tracks the index currently under consideration

            for entry in zip(idx, constraint, d2d):


                if entry[1] == False: #constraint = false

                    # if the constraint is false, then we delete the entry

                    new_idx = np.delete(idx, j)
                    new_d2d = np.delete(d2d, j)

                    # Removes a step from the iterator in case an element was removed.
                    # Because we'll increment the iterator at the end of the loop
                    j-=1

                j+=1


        # cat1_idx is a list of indices into catalog1
        # new_idx is a list of indices into catalog2
        # new_d2d is a list of angular distances between catalog1[cat1_idx] and catalog2[new_idx]

        if show_dynamic_radii is True:

            return([cat1_idx, new_idx, new_d2d, xmatch_rad])

        else:

            return([cat1_idx, new_idx, new_d2d])

    def empty_combined_table(self, table_list):

        """
        Takes a list of astropy tables and returns an empty table with the columns of the old tables (hstacked).
        Used as an auxiliary to generate_full_table()

        Arguments:
            table_list [list of astropy tables]: list of tables to hstack and return list of columns

        Returns:
            First row of hstacked table_list
        """

        combined_table = hstack(table_list)

        combined_table.remove_rows(slice(0, len(combined_table)))

        return(combined_table)

    def merge_tables(self, cat_table_1, cat_table_2):

        """
        Returns a new astropy table that consists of table1 and table2 merged horizontally (columns of table1 coming first).

        Arguments:
            table1 [CatalogTable]: First table to merge
            table2 [CatalogTable]: Second table to merge

        Returns:
            Astropy table consisting of table1 and table2 merged horizontally, with corresponding objects in the same rows.
        """

        # xmatch the two tables
        xmatch_result = self.xmatch(tablecat1,tablecat2)

        table1 = cat_table_1.table
        table2 = cat_table_2.table

        cats1 = cat_table_1.catalogs
        cats2 = cat_table_2.catalogs

        # add catalogs for table 2 to catalogs for table 1, without repeating
        for cat in cats2:
            if cat not in cats1:
                cats1.append(cat)

        # Create table with the columns of table1 and table2
        combined_table = self.empty_combined_table([table1, table2])

        for row_index in range(len(xmatch_result[0])):

            table1_row = [i for i in table1[xmatch_result[0][row_index]]]
            table2_row = [j for j in table2[xmatch_result[1][row_index]]]

            combined_row = table1_row + table2_row

            combined_table.add_row(combined_row, mask=[type(i) is np.ma.core.MaskedConstant for i in combined_row])

        return(CatalogTable(cats1, combined_table))

    # TODO: needs documentation updates
    # TODO: needs updates to support more catalogs
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

    # TODO: needs documentation updates
    # TODO: needs updates to support more catalogs
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
            "designation", "ra", "dec", "ra_error", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag", "bp_rp", "bp_g", "g_rp", "phot_g_n_obs", "phot_g_mean_flux_over_error", "phot_g_mean_flux_error", "phot_g_mean_flux", "phot_bp_mean_mag", "phot_rp_mean_mag"

        2MASS
            "designation", "ra", "dec","err_maj", "err_min", "j_m", "h_m", "k_m", "j_cmsig", "h_cmsig", "k_cmsig", "k_snr"

        AllWISE
            "designation", "ra", "dec", "sigra", "sigdec", "w1mpro", "w2mpro", "w3mpro", "w4mpro", "w1sigmpro", "w2sigmpro", "w3sigmpro", "w4sigmpro", "w4snr"
        """

        print("\nGenerating full table...")

        print("Extracting specified data columns...")

        if gen_small_table:

            gaia_cat = self.gaia["designation", "ra", "dec", "ra_error", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag", "bp_rp", "bp_g", "g_rp", "phot_g_n_obs", "phot_g_mean_flux_over_error", "phot_g_mean_flux_error", "phot_g_mean_flux", "phot_bp_mean_mag", "phot_rp_mean_mag"]

            tmass_cat = self.tmass["designation", "ra", "dec","err_maj", "err_min", "j_m", "h_m", "k_m", "j_cmsig", "h_cmsig", "k_cmsig", "k_snr"]

            allwise_cat = self.allwise["designation", "ra", "dec", "sigra", "sigdec", "w1mpro", "w2mpro", "w3mpro", "w4mpro", "w1sigmpro", "w2sigmpro", "w3sigmpro", "w4sigmpro", "w4snr"]

        else:

            gaia_cat = self.gaia

            tmass_cat = self.tmass

            allwise_cat = self.allwise

        # rename columns to avoid name collisions later
        print("Renaming columns...")

        gaia_cat.rename_column("designation", "gaia_designation")
        gaia_cat.rename_column("ra", "gaia_ra")
        gaia_cat.rename_column("dec", "gaia_dec")

        tmass_cat.rename_column("designation", "2mass_designation")
        tmass_cat.rename_column("ra", "2mass_ra")
        tmass_cat.rename_column("dec", "2mass_dec")

        allwise_cat.rename_column("designation", "allwise_designation")
        allwise_cat.rename_column("ra", "allwise_ra")
        allwise_cat.rename_column("dec", "allwise_dec")

        print("Crossmatching Gaia with 2MASS...")
        gaia_X_tmass = self.merge_tables((("gaia",), gaia_cat), (("2mass",), tmass_cat))


        print("Crossmatching Gaia with AllWISE...")
        gaia_X_allwise = self.merge_tables((("gaia",), gaia_cat), (("allwise",), allwise_cat))


        print("Crossmatching 2MASS with AllWISE...")
        tmass_X_allwise = self.merge_tables((("2mass",), tmass_cat), (("allwise",), allwise_cat))


        gaia_X_tmass_X_allwise = self.merge_tables(gaia_X_tmass, (("allwise",), allwise_cat))
        
        full_table = gaia_X_tmass_X_allwise

        print("Adding objects that do not appear in all catalogs...")

        # objects in gaia_X_tmass_table that are not in gaia_X_tmass_X_allwise, i.e. objects in Gaia and Tmass but not Allwise
        diff1 = self.table_difference(gaia_X_tmass, gaia_X_tmass_X_allwise, "gaia_designation", "gaia_designation", ["allwise"], gaia_cat, tmass_cat, allwise_cat)
        full_table = vstack(full_table, diff1)


        # objects in gaia_cat that are not in full_table, i.e. objects that are in gaia but NOT in (all three catalogs, or both gaia and tmass)
        diff2 = self.table_difference(gaia_cat, full_table, "gaia_designation", "gaia_designation", ["2mass", "allwise"], gaia_cat, tmass_cat, allwise_cat)
        full_table = vstack(full_table, diff2)


        # objects in tmass_X_allwise but not in gaia
        diff3 = self.table_difference(tmass_X_allwise, full_table, "2mass_designation", "2mass_designation", ["gaia"], gaia_cat, tmass_cat, allwise_cat)
        full_table = vstack(full_table, diff3)


        # objects in tmass but not in gaia or allwise
        diff4 = self.table_difference(tmass_cat, full_table, "2mass_designation", "2mass_designation", ["gaia", "allwise"], gaia_cat, tmass_cat, allwise_cat)
        full_table = vstack(full_table, diff4)


        # objects in allwise but not in (gaia OR tmass_X_allwise OR all three OR gaia_X_tmass OR tmass)
        diff5 = self.table_difference(allwise_cat, full_table, "allwise_designation", "allwise_designation", ["gaia", "2mass"], gaia_cat, tmass_cat, allwise_cat)
        full_table  = vstack(full_table, diff5)

        print("Computing user-defined columns...")
        print("Variability...")
        v = np.sqrt(full_table['phot_g_n_obs']/full_table['phot_g_mean_flux_over_error'])
        v.name="variability"
        v.unit="(n_obs / mag)^0.5"
        full_table.add_column(v)

        print("Radial distance...")
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

    # TODO: needs updates to support extraction of rows with data in particular columns
    def extract_ft_part(self, cat_itr, table):
        """
        Returns a copy of the table without the parts that doesn't have data from the surveys in cat_itr.

        Arguments:
            cat_itr [list of strings]: list of catalog names to check for
            table [astropy table]: table to return the relevant portion of

        Returns:
            A copy of table including only the rows that contain data from all surveys in cat_itr
        """

        select_dict = {"gaia": "gaia_designation",

                       "2mass": "2mass_designation",

                       "allwise": "allwise_designation"}

        # Creates an empty copy of full_table

        # edited 2020-05-10
        #reduced_table = full_table.copy()

        #reduced_table.remove_rows(slice(0, len(reduced_table)))
        reduced_table = empty_combined_table(full_table)


        for row_ind in range(len(full_table)):

            add_row = []

            for cat in cat_itr:

                if not(np.ma.is_masked(full_table[row_ind][select_dict[cat]])):

                    add_row.append(True)

                else:

                    add_row.append(False)



            if all(add_row):

                mask = []

                for val in full_table[row_ind]:

                    mask.append(np.ma.is_masked(val))



                reduced_table.add_row(vals=full_table[row_ind],mask=mask)

        return(reduced_table)