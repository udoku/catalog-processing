import numpy as np
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import hstack, vstack, Table, Column
from astropy.io import ascii

import os

from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus
#from astroquery.vizier import Vizier
#from astroquery.skyview import SkyView

from . import path, dpath, cpath, tpath, logger, out, info_out

from CatalogTable import *

class catalogProcessing:
    """
    Class for querying astronomical databases.

    Contains methods for accessing data from Gaia and infrared surveys,
    as well as for manipulating the resulting data.

    Instantiated with

    @ra: Right ascension of the target coordinates
    @dec: Declination of the target coordinates
    @radius: Radius within which to search
    @parallax [optional]: Target parallax.

    If @parallax is provided, then @radius is assumed to be in pc, and is
    converted to degrees. If @parallax is not provided, @radius is assumed
    to be in degrees.
    """

    def __init__(self, ra, dec, radius, parallax=None, load=None):
        """
        Initialize a CatalogProcessing object.

        Performs catalog queries of Gaia, 2Mass, and WISE.

        Arguments:
            ra: Right Ascension of the target
            dec: Declination of the target
            radius: radius of the target (deg if @parallax is not supplied, pc if it is)
            parallax: Parallax of target, in mas
        """

        if load is not None:

            info_out("Loading data from directory: " + load)
            try:
                self.skycoord, self.radius, self.parallax = self.load_target(load)
                self.gaia, self.tmass, self.allwise = self.load_tables(load)

            except:
                info_out("The specified directory does not exist:")
                data_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                info_out(data_dir + "/" + load)

        else:

            self.skycoord = SkyCoord(ra, dec)
            self.radius = radius
            self.parallax = parallax

            if self.parallax != None:

                self.radius = self.get_radius()

            out("Initializing CatalogProcessing:")
            info_out("ra={}".format(ra, dec, self.radius, parallax))
            info_out("dec={}".format(dec))
            info_out("radius={} deg".format(self.radius))
            info_out("parallax={} mas".format(parallax))

            # Query catalogs and save the data to file
            info_out("Querying Gaia...")
            self.gaia = self.gaia_query()
            self.gaia.catalogs = ["gaia"]
            out("Done querying Gaia!\n")
            info_out("Querying 2MASS...")
            self.tmass = self.ir_query("tmass")
            out("Done querying 2MASS!\n")
            info_out("Querying AllWISE...")
            self.allwise = self.ir_query("allwise")
            out("Done querying AllWISE!\n")
            # other catalogs to be added

    def get_radius(self):

        """
        Converts a parallax (mas) and a physical search radius (pc) to an angular search radius.

        Arguments:
            [none]

        Returns:
            Angular search radius in degrees.
        """

        # since parallax is given in milliarcseconds, we have to multiply the reciprocal of the parallax by 1000 to get distance in parsecs
        # in general, distance (pc) ~ 1 / parallax (arcsec)
        dist = 1000/(self.parallax)

        radius_rad = self.radius / dist

        return(radius_rad * 57.2958) # magic number: degrees per radian

    
    def load_target(self, dir):
        """
        Loads ra, dec, radius, parallax data from a previously created directory

        Arguments:
            dir [string]: directory in which to look

        Returns:
            ra, dec, radius, parallax
        """

        data_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 

        try:
            with open(data_dir + "/" + dir + "/README", 'r') as f:
                out(dir + ": " + f.readline())
                ra = f.readline().lstrip('ra=').rstrip('\n')
                out(ra)
                dec = f.readline().lstrip('dec=').rstrip('\n')
                out(dec)
                radius = float(f.readline().lstrip('radius=').rstrip(' deg\n'))
                out(radius)
                parallax = float(f.readline().lstrip('parllax=').rstrip(' mas\n'))
                out(parallax)

        except:
            out("Specified directory does not contain the required data:")
            out(dir)
            return

        return (SkyCoord(ra, dec, unit=(u.hourangle, u.deg)), radius, parallax)


    def load_tables(self, dir):
        """
        Loads gaia, tmass, and allwise tables from a previously created directory

        Arguments:
            dir [string]: directory in which to look

        Returns:
            (gaia table, tmass table, allwise table)
            All as CatalogTables
        """
        
        data_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        try:
            gaia_name = data_dir + "/" + dir + "/tables/gaia_query.dat"
            out(gaia_name)
            gaia_table = Table.read(gaia_name, format='ascii')
            gaia = CatalogTable(['gaia'], gaia_table)
            out("read Gaia table")

            gaia_fname = tpath + "/gaia_query.dat"
            gaia_table.write(gaia_fname, format='ascii.ecsv')
            out("wrote Gaia table")

            tmass_name = data_dir + "/" + dir + "/tables/tmass_query.dat"
            tmass_table = Table.read(tmass_name, format='ascii')
            tmass = CatalogTable(['tmass'], tmass_table)
            out("read Tmass table")

            tmass_fname = tpath + "/tmass_query.dat"
            tmass_table.write(tmass_fname, format='ascii.ecsv')
            out("wrote Tmass table")

            allwise_name = data_dir + "/" + dir + "/tables/allwise_query.dat"
            allwise_table = Table.read(allwise_name, format='ascii')
            allwise = CatalogTable(['allwise'], allwise_table)
            out("read Allwise table")

            allwise_fname = tpath + "/allwise_query.dat"
            allwise_table.write(allwise_fname, format='ascii.ecsv')
            out("wrote Allwise table")
        
        except:
            out("The specified directory does not contain the required data [2]:")
            out(dir)
            return

        return (gaia, tmass, allwise)


    # TODO: integrate with irquery
    def gaia_query(self):

        """
        Performs a query on gaia DR2 using self.skycoord and self.radius.
        Note: "Exception: 500" tends to happen when the coordinates used in search_string is not formatted correctly.

        Arguments:
            [none]

        Returns:
            CatalogTable with query result.
        """

        catalog = ["gaia"]

        ra = self.skycoord.ra.degree

        dec = self.skycoord.dec.degree


        #TODO: include code to write diagnostic info to log file incl. ra, dec, radius, search string

        search_string = "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{0},{1},{2}))=1;".format(ra, dec, self.radius)

        #try:

        out("Creating query...")
        out(search_string)
        job = Gaia.launch_job_async(search_string, dump_to_file=False)

        #except gaierror as e:

        #    if str(e) != "[Errno 11001] getaddrinfo failed":

        #        raise

        #    else:

        #        out("This error is typically raised when there is no internet connection.")

        #        raise

        out("Retrieving results...")
        query_results = job.get_results()

        out("Results retrieved.")
        info_out(str(len(query_results['designation'])) + " sources detected.")

        # write Gaia query results to file
        fname = tpath + "/gaia_query.dat"
        query_results.write(fname, format='ascii.ecsv')

        return(CatalogTable(catalog,query_results))

    # TODO: Needs documentation update/standardization
    def ir_query(self, ircat_name, view_adql=True):

        """
        Performs a TAP+ query to irsa.ipac and returns the results.

        Main reason for this function is to have basis for easily incorporating TAP queries in a general query function later.
        -----------------------

        Example catalogs:

        "allwise_p3as_psd": AllWISE Source Catalog

        "fp_psc": 2MASS Point Source Catalog

        "glimpse_s07": GLIMPSE I Spring 07 Catalog (Spitzer)

        "cosmos_phot": COSMOS Photometry Catalog

        "iraspsc": IRAS Point Source Catalog

        -----------------------

        Arguments:
            ircat_name [string]: name of the catalogue in the irsa.ipac TAP spec
            view_adql [bool]: If True, prints the adql of the query generated.

        Returns:
            CatalogTable with query results.
        """

        catalogs = {"allwise": "allwise_p3as_psd", "tmass": "fp_psc", "glimpse": "glimpse_s07" }

        catalog = catalogs[ircat_name]

        search_string = "SELECT * FROM {} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{},{},{}))=1 ".format(catalog, str(self.skycoord.ra.degree), str(self.skycoord.dec.degree), self.radius)

        if view_adql:
            out(search_string)

        ipac = TapPlus(url="https://irsa.ipac.caltech.edu/TAP/")

        out("Creating query...")
        job = ipac.launch_job_async(search_string)

        out("Retrieving results...")
        query_results = job.get_results()

        out("Results retrieved.")
        info_out(str(len(query_results['designation'])) + " sources detected.")

        # Write query results to file.
        fname = tpath + "/" + ircat_name + "_query.dat"
        query_results.write(fname, format='ascii.ecsv')

        return(CatalogTable([ircat_name],query_results))

    # TODO: needs documentation updates
    # TODO: needs better error handling
    # TODO: read in survey definitions from file
    def xmatch(self, catalog1, catalog2, rad=-1, show_dynamic_radii=True):
        """
        Crossmatches two catalogs.

        Uses SkyCoord.match_to_catalog_sky() method which implements a nearest neighbour algorithm.

        valid surveynames: "2mass", "allwise", "gaia"

        Arguments:
            catalog1 [CatalogTable]: first catalog to crossmatch
            catalog2 [CatalogTable]: second catalog to crossmatch
            rad [float]: maximum radius to crossmatch. By default, specifies dynamic radius using the associated error.
            show_dynamic_radii [bool]: 

        


        cat1, cat2: full query results (astropy tables resulting from querying catalogs)
        colnames1, colnames2: list of the column names of ra dec and respective errors in cat1, cat2. Either [ra, dec] or [ra, dec, err_ra, err_dec] for


        Returns: 
            [cat1_idx, new_idx, new_d2d]
            or if show_dynamic_radii=True:
            [cat1_idx, new_idx, new_d2d, match radii]

            cat1_idx: index of the element in cat2 that is nearest to the corresponding element in cat1
            (cat1_idx[0] = 69 means that the first element in cat1 is matched towards the 69th element in cat2)

            new_idx: list of indices into catalog2

        """

        surveys_coord_colnames = {
            "2mass":["2mass_ra", "2mass_dec", "err_maj", "err_min"],
            "allwise":["allwise_ra", "allwise_dec", "sigra", "sigdec"],
            "gaia":["gaia_ra", "gaia_dec", "ra_error", "dec_error"]}


        colnames1 = surveys_coord_colnames[catalog1.catalogs[0]]
        cat1 = catalog1.table

        colnames2 = surveys_coord_colnames[catalog2.catalogs[0]]
        cat2 = catalog2.table

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

            constraint = d2d < rad#*u.arcsec


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

                # indpair[1] is 2d distance (projected onto the sky) between 
                # the source in catalog 1 and the matched source in catalog 2
                dist = indpair[1].arcsec

                #const = coordinates.Angle(match_constraint*u.arcsec)

                # TODO: why do we know match_constraint is in arcsec?
                if dist > match_constraint:
                    # if this is true then remove the crossmatch entry

                    new_idx = np.delete(new_idx, j)
                    new_d2d = np.delete(new_d2d, j)

                    # subtract 1 from the index under consideration because we 
                    # just deleted the entry we were looking at, and we'll 
                    # increment j at the end of the loop - this ensures
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
        Generates an empty table with the columns of all tables in table_list hstacked.

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
        xmatch_result = self.xmatch(cat_table_1,cat_table_2)

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
    # TODO: needs updates to support more/arbitrary catalogs
    def table_difference(self, table1, table2, id, not_in_table1, gaia_cat, tmass_cat, allwise_cat):
        """
        Returns an astropy table (with the dimenstions of table2) that consists of the objects in table 1 but not in table 2, where the surveys composing table1 is fully contained in table2.

        Arguments:
            table1 [CatalogTable]: first table to merge
            table2 [CatalogTable]: second table to merge
            id [string]: column by which to identify matching sources
            not_in_table1 [list of strings]: list of names of the surveys (gaia, 2mass, allwise) not in table1 but in table2
            # TODO: eliminate this parameter using the .catalogs data object of the CatalogTables

        """
        dtype = table2.table.dtype
        dtype2 = []
        for n in dtype.names:
            dtype2.append(dtype[n])

        diff_table = Table(names=table2.table.colnames, dtype=dtype2)

        t2cols = table2.table.colnames
        t1cols = table1.table.colnames

        for row in table1.table:

            if row[id] in table2.table[id]:

                pass

            else:

                full_row = []

                for col in t2cols:
                    if col in t1cols:
                        full_row.append(row[col])
                    else:
                        full_row.append(-1)


                mask = [i is -1 for i in full_row]

                diff_table.add_row(full_row, mask=mask)

        return(diff_table)

    # TODO: needs documentation updates
    # TODO: needs updates to support more catalogs
    def generate_full_table(self, gen_small_table = False, load=None):
        """
        Generates the full table from the gaia, tmass, and allwise queries.

        Xmatch Gaia and tmass
        Xmatch Gaia with allwise
        Xmatch tmass with allwise
        Xmatch (Gaia x tmass) with allwise

        We then append:
        sources in all three catalogs
        sources in two of the catalogs
        sources in just one catalog

        Arguments:
            gen_small_table [bool]: Whether to only use a subset of available columns. False by default.
            load [string]: Directory from which to load an already generated full_table. None by default.
        """

        out("\nGenerating full table...")

        out("Extracting specified data columns...")

        if gen_small_table:

            # Gaia columns in small_table
            gaia_cat_table = self.gaia.table["designation", "ra", "dec", "ra_error", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag", "bp_rp", "bp_g", "g_rp", "phot_g_n_obs", "phot_g_mean_flux_over_error", "phot_g_mean_flux_error", "phot_g_mean_flux", "phot_bp_mean_mag", "phot_rp_mean_mag"]
            gaia_cat = CatalogTable(["gaia"],gaia_cat_table)

            # 2Mass columns in small_table
            tmass_cat_table = self.tmass.table["designation", "ra", "dec","err_maj", "err_min", "j_m", "h_m", "k_m", "j_cmsig", "h_cmsig", "k_cmsig", "k_snr"]
            tmass_cat = CatalogTable(["2mass"],tmass_cat_table)

            # Allwise columns in small_table
            allwise_cat_table = self.allwise.table["designation", "ra", "dec", "sigra", "sigdec", "w1mpro", "w2mpro", "w3mpro", "w4mpro", "w1sigmpro", "w2sigmpro", "w3sigmpro", "w4sigmpro", "w4snr"]
            allwise_cat = CatalogTable(["allwise"],allwise_cat_table)

        else:

            gaia_cat = self.gaia

            tmass_cat = self.tmass

            allwise_cat = self.allwise

        # rename columns to avoid name collisions later
        out("Renaming columns...")

        gaia_cat.table.rename_column("designation", "gaia_designation")
        gaia_cat.table.rename_column("ra", "gaia_ra")
        gaia_cat.table.rename_column("dec", "gaia_dec")

        tmass_cat.table.rename_column("designation", "2mass_designation")
        tmass_cat.table.rename_column("ra", "2mass_ra")
        tmass_cat.table.rename_column("dec", "2mass_dec")

        allwise_cat.table.rename_column("designation", "allwise_designation")
        allwise_cat.table.rename_column("ra", "allwise_ra")
        allwise_cat.table.rename_column("dec", "allwise_dec")



        out("Crossmatching Gaia with 2MASS...")
        gaia_X_tmass = self.merge_tables(gaia_cat, tmass_cat)

        out("Crossmatching Gaia with AllWISE...")
        gaia_X_allwise = self.merge_tables(gaia_cat, allwise_cat)

        out("Crossmatching 2MASS with AllWISE...")
        tmass_X_allwise = self.merge_tables(tmass_cat, allwise_cat)

        out("Crossmatching all three catalogs...")
        gaia_X_tmass_X_allwise = self.merge_tables(gaia_X_tmass,  allwise_cat)
        
        info_out(str(len(gaia_X_tmass_X_allwise.table['gaia_designation'])) + " sources in all three catalogs.")

        full_table = gaia_X_tmass_X_allwise

        out("Adding objects that do not appear in all catalogs...")

        # sources in gaia_X_tmass that are not in gaia_X_tmass_X_allwise
        # i.e. sources in Gaia and Tmass but not Allwise
        diff1 = self.table_difference(gaia_X_tmass, full_table, "gaia_designation", ["allwise"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff1['gaia_designation'])) + " sources in Gaia and 2MASS but not AllWISE.")
        full_table.table = vstack([full_table.table, diff1])

        # sources in gaia_X_allwise that are not in gaia_X_tmass_X_allwise
        # i.e. sources in Gaia and Allwise but not Tmass.
        diff2 = self.table_difference(gaia_X_allwise, full_table, "gaia_designation", ["2mass"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff2['gaia_designation'])) + " sources in Gaia and AllWISE but not 2MASS.")
        full_table.table = vstack([full_table.table, diff2])

        # objects in tmass_X_allwise that are not in gaia_X_tmass_X_allwise, 
        # i.e. sources in tmass and allwise but not gaia.
        diff3 = self.table_difference(tmass_X_allwise, full_table, "allwise_designation", ["gaia"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff3['allwise_designation'])) + " sources in 2MASS and AllWISE but not Gaia.")
        full_table.table = vstack([full_table.table, diff3])

        # objects in gaia but not yet in full_table
        diff4 = self.table_difference(gaia_cat, full_table, "gaia_designation", ["2mass", "allwise"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff4['gaia_designation'])) + " sources in Gaia only.")
        full_table.table = vstack([full_table.table, diff4])

        # objects in tmass but not yet in full_table
        diff5 = self.table_difference(tmass_cat, full_table, "2mass_designation", ["gaia", "allwise"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff5['2mass_designation'])) + " sources in 2mass only.")
        full_table.table = vstack([full_table.table, diff5])

        # objects in allwise but not yet in full_table
        diff6 = self.table_difference(allwise_cat, full_table, "allwise_designation", ["gaia", "2mass"], gaia_cat, tmass_cat, allwise_cat)
        info_out(str(len(diff6['allwise_designation'])) + " sources in AllWISE only.")
        full_table.table  = vstack([full_table.table, diff6])


        out("Computing user-defined columns...")

        out("Variability...")
        v = np.sqrt(full_table.table['phot_g_n_obs']/full_table.table['phot_g_mean_flux_over_error'])
        v.name="variability"
        v.unit="(n_obs / mag)^0.5"
        full_table.table.add_column(v)

        out("Radial distance...")
        d = 1000 / full_table.table['parallax']
        d.name = "radial_distance"
        d.unit="pc"
        full_table.table.add_column(d)

        info_out("Full table generated.")

        fname = tpath + "/full_table.dat"
        full_table.table.write(fname, format='ascii.ecsv')

        return(full_table)

    # TODO: needs updates to support extraction of rows with data in particular columns
    def extract_ft_part(self, cat_itr, table):
        """
        Returns a copy of @table without the parts missing data from the surveys in cat_itr.

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
        reduced_table = self.empty_combined_table(table)


        for row_ind in range(len(table)):

            add_row = []

            for cat in cat_itr:

                if not(np.ma.is_masked(table[row_ind][select_dict[cat]])):

                    add_row.append(True)

                else:

                    add_row.append(False)



            if all(add_row):

                mask = []

                for val in table[row_ind]:

                    mask.append(np.ma.is_masked(val))



                reduced_table.add_row(vals=table[row_ind],mask=mask)

        return(reduced_table)


    def get_surveys(self, survey_list, field_side_length):
        """
        Retrieves images of the target from surveys

        Full description of all the surveys availible can be found here:

        https://skyview.gsfc.nasa.gov/current/cgi/survey.plhttps://skyview.gsfc.nasa.gov/current/cgi/survey.pl

        Some surveynames are:
        "2MASS-J", "2MASS-H", "2MASS-K"
        "WISE 3.4", "WISE 4.6", "WISE 12", "WISE 22"
        "DSS"
        "SDSSg", "SDSSi", "SDSSr", "SDSSu", "SDSSz"

        Arguments:
            survey_list [list of strings]: list of SkyView surveynames.

        Returns:
            list of astropy.fits.HDUList objects with attribute "".data" that can be passed to plt.imshow().
        """

        out("fetching surveys:")
        out(survey_list)
        out("with radius:")
        out(field_side_length)

        center=self.skycoord.ra.degree + self.skycoord.dec.degree

        return(SkyView.get_images(position=center, survey=survey_list, radius=(field_side_length)))