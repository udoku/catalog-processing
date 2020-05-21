from Photometry import *

from Astrometry import *

from CatalogProcessing import *

from code import path, logger, out

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
        list have an attribute "".data" that can be passed to plt.imshow().

        """

        out("fetching surveys:")
        out(survey_list)
        out("with radius:")
        out(field_side_length)

        center=self.ra + self.dec

        return(SkyView.get_images(position=center, survey=survey_list, radius=(field_side_length)))





    def plot_surveys(self, img_list, formatted_coord_list, survey_list, plot_coords=True):



        #out(img_list)

        wcs_list = []

        for i in img_list:

            #out(i[0])

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


            #out((formattedcoord_list[i].ra, formattedcoord_list[i].dec))

            if plot_coords:

                for cat in formatted_coord_list:

                    ax.plot(cat.ra, cat.dec, 'o', transform=ax.get_transform('icrs'), mec='r', mfc='none')





    def plot_images(self, table, survey_list, cats_to_plot=None):

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

        # for row in self.full_table.table["ra", "dec"]:

        #     coord_list.append((row[0], row[1]))

        out("Adding coordinates of catalog objects...")

        try:

            coord_list = []

            for cat in cats_to_plot:

                out("Adding objects from " + str(cat))

                cat_coord_list = []
                count = 0

                ra_str = str(cat) + "_ra"
                dec_str = str(cat) + "_dec"

                for row in table.table[ra_str, dec_str]:

                    if row[0] != 0:
                        count += 1
                        cat_coord_list.append((row[0], row[1]))
                
                coord_list.append(cat_coord_list)

                out("Added " + str(cat) + " coords: " + str(count) + " objects added.")


            #ra_min = min(min([coords[0] for coords in coord_list[cal]]) for cal in range(len(coord_list)))

            #ra_max = max(max([coords[0] for coords in coord_list[cal]]) for cal in range(len(coord_list)))

            #dec_min = min(min([coords[1] for coords in coord_list[cal]]) for cal in range(len(coord_list)))

            #dec_max = max(max([coords[1] for coords in coord_list[cal]]) for cal in range(len(coord_list)))




            #field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            out("Formatting coordinates...")

            formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

            out("Finding minimum and maximum RA and DEC...")

            ra_min = min(formatted_coord_list[0].ra)

            ra_max = max(formatted_coord_list[0].ra)

            dec_min = min(formatted_coord_list[0].dec)

            dec_max = max(formatted_coord_list[0].dec)

            out("Calculating field side length...")

            field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            #formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

            out("Retrieving images...")

            img_list = self.get_surveys(survey_list, field_side_length)

            self.plot_surveys(img_list, formatted_coord_list, survey_list, plot_coords=True)


        except TypeError as e:

            out("TypeError occurred!")

            out(e)

            '''ra_min = min(self.full_table.table["gaia_ra"])

            ra_max = max(self.full_table.table["gaia_ra"])

            dec_min = min(self.full_table.table["gaia_dec"])

            dec_max = max(self.full_table.table["gaia_dec"])

            #out(ra_min, ra_max)

            field_side_length = max((ra_max - ra_min, dec_max - dec_min))

            #out(field_side_length)

            img_list = self.get_surveys(survey_list, field_side_length*u.deg)

            self.plot_surveys(img_list, [], survey_list, plot_coords=False)'''



        '''#formatted_coord_list = self.format_coords(coord_list)
        formatted_coord_list = [self.format_coords(cl) for cl in coord_list]

        #out(type(formatted_coord_list))
        #out(formatted_coord_list[0])

        ra_min = min(formatted_coord_list[0].ra)

        ra_max = max(formatted_coord_list[0].ra)

        dec_min = min(formatted_coord_list[0].dec)

        dec_max = max(formatted_coord_list[0].dec)

        #out(ra_min)
        #out(ra_max)
        #out(ra_max - ra_min)



        field_side_length = max((ra_max - ra_min, dec_max - dec_min))

        #out(field_side_length)


        img_list = self.get_surveys(survey_list, field_side_length)

        out("Number of surveys:")
        out(len(survey_list))
        out("Number of images:")
        out(len(img_list))

        self.plot_surveys(img_list, formatted_coord_list, survey_list)'''

    def plot_hist(self, colname, xlable, cut=None):

        plt.figure()

        if cut != None:

            phot = Photometry(self.full_table)

            column = phot.apply_cut(cut)[0].table[colname]

            #out("column:")

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

            #self.full_table.table[colname].pout(max_lines=-1)



            hist = astrometry.make_hist(self.full_table.table[colname])



            #out(fit)



            x = np.linspace(min(self.full_table.table[colname]),max(self.full_table.table[colname]),1000)



            #plt.title('',fontsize=16)

            plt.xlabel(xlable,fontsize=14)

            plt.ylabel('# count',fontsize=14)



            plt.hist(self.full_table.table[colname], bins=hist[0][1])

            plt.plot(x, astrometry.gaussian(x, fit[0], fit[1], fit[2]))

            plt.show()

    def plot_error(self, parameter):

        try:
            good_data = self.full_table.table["gaia_ra"] > 0
            gmag = self.full_table.table['phot_g_mean_mag']
            param_data = self.full_table.table[parameter]
            param_error = self.full_table.table[parameter + "_error"]
        except:
            out("Given parameter does not appear in the imported data. Check parameter name for consistency with full_table columns.")

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

        #out(self.full_table.table[colname1], self.full_table.table[colname1])



        fig = plt.figure(figsize=(20,20))

        fig.clf()



        #fig, axs = plt.subplots((len(colname_list)))



        for i in range(len(colname_list)):

            ind = i+1

            ax = fig.add_subplot(len(colname_list), 2, ind)



            ax.set_xlabel(colname_list[i][0] + " [{}]".format(self.full_table.table[colname_list[i][0]].unit), fontsize=14)

            ax.set_ylabel(colname_list[i][1] + " [{}]".format(self.full_table.table[colname_list[i][1]].unit), fontsize=14)

            #out("hi")

            ax.scatter(self.full_table.table[colname_list[i][0]], self.full_table.table[colname_list[i][1]], s=8, marker="o")

        plt.show()

    def plot(self, colname1, colname2, xlim=None, ylim=None, squared=False, invert_y=False, invert_x=False):

        """
        plots colname1 on the x-axis and colname2 on the y-axis.

        colname1/2: tuples on two strings, (colname, label)
        """

        plt.figure()

        if type(colname1) is str:

            x = self.full_table.table[colname1]

            plt.xlabel(colname1 + " [{}]".format(self.full_table.table[colname1].unit), fontsize=14)

        elif type(colname1[0]) is str:

            #out("check")

            x = self.full_table.table[colname1[0]]

            plt.xlabel(colname1[1] + " [{}]".format(self.full_table.table[colname1[0]].unit), fontsize=14)

        else:

            x = colname1[0]

            plt.xlabel(colname1[1], fontsize=14)



        if type(colname2) is str:

            y = self.full_table.table[colname2]

            plt.ylabel(colname2 + " [{}]".format(self.full_table.table[colname2].unit), fontsize=14)

        elif type(colname2[0]) is str:

            y = self.full_table.table[colname2[0]]

            plt.ylabel(colname2[1] + " [{}]".format(self.full_table.table[colname2[0]].unit), fontsize=14)

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

        #out(type(x), type(y))



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

                x = self.full_table.table[col_tup[0]]

                ax.set_xlabel(col_tup[0] + " [{}]".format(self.full_table.table[col_tup[0]].unit), fontsize=14)

            #elif type(col_tup[2]) is str:

            #    #out("check")

            #    x = self.full_table.table[col_tup[0]]

            #    ax.set_xlabel(col_tup[2] + " [{}]".format(self.full_table.table[col_tup[0]].unit), fontsize=14)

            else:

                x = col_tup[0]

                ax.set_xlabel(col_tup[2], fontsize=14)



            if type(col_tup[1]) is str:

                y = self.full_table.table[col_tup[1]]

                plt.ylabel(col_tup[1] + " [{}]".format(self.full_table.table[col_tup[1]].unit), fontsize=14)

            #elif type(col_tup[3]) is str:

            #    y = self.full_table.table[col_tup[1]]

            #    plt.ylabel(colname2[1] + " [{}]".format(self.full_table.table[colname2[0]].unit), fontsize=14)

            else:

                y = col_tup[1]

                plt.ylabel(col_tup[3], fontsize=14)



            #if lims != None:

            #    ax.set_xlim(lims[i][0])

            #    ax.set_ylim(lims[i][1]

            #out(type(x), type(y))



            ax.scatter(x,y, s=8, marker="o")

        plt.show()

    def plot_double_zoom(self, col1, col2):
        """
        Not yet implemented
        """

        pass

    def plot_tuples(self, title, plot_list, xlabel, ylabel, xlim=None, ylim=None, invert_x=False, invert_y=False, squared=False):

        """
        Plots each list in the tuples in plot_list against each other in different colors. Plots up to six different lists of tuples.

        Arguments:
            title [string]: title of the plot
            plot_list [list of tuples]: list of tuples of two lists that represent x and y of different datasets to be plotted together.
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
        ax.set_title(title)

        # apply axis labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # create ordered list of colors
        colors = ["k", "r", "b", "g", "m", "y"]

        # scatter plot the lists on the subfigure object
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

        self.plot_tables('stand-in title', (cut_true, cut_false), col1, col2, xlim=xlim, ylim=ylim, squared=squared, invert_y=invert_y, invert_x=invert_x)
        
        #xlabel = str(col1[1]) + " [{}]".format(self.full_table.table[col1[0]].unit)
        #ylabel = str(col2[1]) + " [{}]".format(self.full_table.table[col2[0]].unit)

        #self.plot_removed([(cut_true[col1[0]], cut_true[col2[0]]), (cut_false[col1[0]], cut_false[col2[0]])],
        #        xlabel, ylabel, xlim, ylim, squared=squared, invert_y = invert_y, invert_x = invert_x)

    def plot_tables(self, title, tables, col1, col2, xlim=None, ylim=None, invert_x=False, invert_y=False, squared=False):

        """
        Plots col1 against col2 for each table in the list tables. Plots are superimposed and differentiated by color.

        Arguments:
            tables [list of astropy tables]: col1 will be plotted against col2 for each table in this list
            col1 [tuple of strings]: (column name, column description) to plot on the x-axis
            col2 [tuple of strings]: (column name, column description) to plot on the y-axis
            xlim [tuple]: (lower x limit, upper x limit) for plot. Default: None
            ylim [tuple]: (lower y limit, upper y limit) for plot. Default: None
            invert_x [boolean]: whether or not to invert the x axis. Default: False
            invert_y [boolean]: whether or not to invert the y axis. Default: False
            squared [boolean]: whether or not to force generation of a square plot. Default: False

        Returns:
            [none]
        """

        try:
        
            xlabel = str(col1[1]) + " [{}]".format(tables[0].table[col1[0]].unit)
            ylabel = str(col2[1]) + " [{}]".format(tables[0].table[col2[0]].unit)

        except KeyError as e:
            out("The specified colname does not match a column in the provided tables. " + 
                  "This might occur if you have provided a (column name, column description) tuple with an invalid column name, or if you have provided a single string instead of a tuple.")
            out(e)

        tuples = []

        for t in tables:
            try:
                tuples.append((t.table[col1[0]], t.table[col2[0]]))
            except KeyError as e:
                out("Table " + str(tables.index(t)) + " does not include column " + e)
        
        self.plot_tuples(title, tuples, xlabel, ylabel, xlim = xlim, ylim = ylim, invert_x = invert_x, invert_y = invert_y, squared=squared)

    def plot_diagnostics(self):
        """
        generates hard-coded plots to display useful diagnostic information
        """

        
        
        # define variables

        out("Computing photometric variables...")
        M_g = [g + 5 - 5*np.log10(1000/p) for g,p in zip(self.full_table.table["phot_g_mean_mag"], self.full_table.table["parallax"])]
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
            out("Applying photometric cuts...")
            cut_true = self.full_table.table[pd.eval("""(-1.81 -2*1.31 < self.phot.full_table["pmra"]) & (self.phot.full_table["pmra"] < -1.81 +2*1.31) & (-2.67 -2*1.44 < self.phot.full_table["pmdec"]) & (self.phot.full_table["pmdec"] < -2.67 +2*1.44)""")]
            cut_false = self.full_table.table[~pd.eval("""(-1.81 -2*1.31 < self.phot.full_table["pmra"]) & (self.phot.full_table["pmra"] < -1.81 +2*1.31) & (-2.67 -2*1.44 < self.phot.full_table["pmdec"]) & (self.phot.full_table["pmdec"] < -2.67 +2*1.44)""")]

            cut_mg_true = [g + 5 - 5*np.log10(1000/p) for g,p in zip(cut_true["phot_g_mean_mag"], cut_true["parallax"])]
            cut_bprp_true = cut_true["bp_rp"]
            cut_mg_false = [g + 5 - 5*np.log10(1000/p) for g,p in zip(cut_false["phot_g_mean_mag"], cut_false["parallax"])]
            cut_bprp_false = cut_false["bp_rp"]

            cut_1s = """(self.full_table.table["pmra"] < -1.81 +1.31) & (self.full_table.table["pmra"] > -1.81 -1.31) & (self.full_table.table["pmdec"] < -2.67 +1.44) & (self.full_table.table["pmdec"] > -2.67 -1.44)"""

            cut_2s = """(self.full_table.table["pmra"] < -1.81 +2*1.31) & (self.full_table.table["pmra"] > -1.81 -2*1.31) & (self.full_table.table["pmdec"] < -2.67 +2*1.44) & (self.full_table.table["pmdec"] > -2.67 -2*1.44)"""

            cut_outliers = """(self.full_table.table["pmra"] < 15) & (self.full_table.table["pmra"] > -15) & (self.full_table.table["pmdec"] < 15) & (self.full_table.table["pmdec"] > -15)"""

            cut_plx_1s = """(self.full_table.table["parallax"]< 1.04+0.23) & (self.full_table.table["parallax"]> 1.04-0.23)"""



            # generate and display plots with documentation

            out("2mass RA vs 2mass Dec. J-H > 0.7 shown in black, J-H < 0.7 shown in red.")
            self.cut_and_plot("J-H > 0.7", ("2mass_ra", "2mass RA"), ("2mass_dec", "2mass Dec"), squared=True, invert_x=True)

            out("Gaia RA vs Gaia Dec.")
            self.plot(("gaia_ra", "Gaia RA"), ("gaia_dec", "Gaia Dec"), squared=True, invert_x=True )

            out("Gaia PM RA vs Gaia PM Dec.")
            self.plot(("pmra", "pm RA (Gaia)"), ("pmdec", "pm Dec (Gaia)"), xlim=(-30,30), ylim=(-30,30), squared=True)

            out("Gaia PM RA vs Gaia PM Dec (closer detail).")
            self.plot(("pmra", "pm RA (Gaia)"), ("pmdec", "pm Dec (Gaia)"), xlim=(-10,10), ylim=(-10,10), squared=True)

            out("Parallax vs Gaia Dec.")
            self.plot(("parallax", "Parallax"), ("gaia_dec", "Dec (Gaia)"), xlim=(0,5))

            #out("what does double_plot do?")
            #self.double_plot([("pmra", "pmdec"), ("pmra", "pmdec")], [(-20,20), (-10, 10)])

            out("BP/RP vs G Mean Magnitude.")
            self.plot("bp_rp", "phot_g_mean_mag")

            out("BP/RP vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot("bp_rp", (M_g, "$M_G [mag]$"))

            out("G-K Magnitude vs G Mean Magnitude.")
            self.plot((g_k, "G-K [mag]"), "phot_g_mean_mag")

            out("G-K Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot((g_k, "G-K [mag]"), (M_g, "$M_G [mag]$"))

            out("J-K Magnitude vs J-M Magnitude.")
            self.plot((j_k, "J-K [mag]"), "j_m")

            out("W1-W2 Magnitude vs W1 Magnitude.")
            self.plot((w1_w2, "W1-W2 [mag]"), "w1mpro")

            out("J-H Magnitude vs H-K Magnitude.")
            self.plot((j_h, "J-H [mag]"), (h_k, "H-K [mag]"))

            out("G-H Magnitude vs H-K Magnitude.")
            self.plot((g_h, "G-H [mag]"), (h_k, "H-K [mag]"))

            out("K-W2 Magnitude vs G-K Magnitude.")
            self.plot((k_w2, "K - W2 [mag]"), (g_k, "G - K [mag]")) 

            out("BP/RP Magnitude vs G Mean Magnitude.")
            self.plot("bp_rp", "phot_g_mean_mag", invert_y=True)

            out("2mass RA vs 2mass Dec. K-W1 > 0.2 shown in black, K-w1 < 0.2 shown in red.")
            self.cut_and_plot("K-W1 > 0.2", ("2mass_ra", "RA"), ("2mass_dec", "Dec"), squared=True)

            out("PMRA vs PMDec. -3.12 < PMRA < -0.5 AND -4.17 < PMDec < -1.23 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_1s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-10,10), ylim=(-10,10), squared=True)

            out("PMRA vs PMDec - closer detail. -3.12 < PMRA < -0.5 AND -4.17 < PMDec < -1.23 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_1s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-6,1), ylim=(-6,1), squared=True)

            out("PMRA vs PMDec. -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_2s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-10,10), ylim=(-10,10), squared=True)

            out("PMRA vs PMDec - closer detail. -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21 shown in black, otherwise shown in red. Units in mas/yr.")
            self.cut_and_plot(cut_2s, ("pmra", "pm RA [$mas\ yr^{-1}$]"), ("pmdec", "pm Dec [$mas\ yr^{-1}$]"), xlim=(-6,1), ylim=(-6,1), squared=True)

            out("Parallax vs PMRA. 0.81 < Parallax < 1.27 shown in black, otherwise shown in red.")
            self.cut_and_plot(cut_plx_1s, ("parallax", "Parallax"),("pmra", "pm RA"), xlim=(0,5), squared=True)

            out("PMRA histogram. Data shown satisfy -15 < PMRA < 15 AND -15 < PMDec < 15. Units in mas/yr.")
            self.plot_hist("pmra", "pm RA", cut=cut_outliers)

            out("PMDec histogram. Data shown satisfy -4.43 < PMRA < 0.81 AND -5.55 < PMDec < 0.21. Units in mas/yr.")
            self.plot_hist("pmdec", "pm Dec", cut=cut_2s)

            out("BP-RP Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax ))")
            self.plot_removed([(cut_bprp_true, cut_mg_true), (cut_bprp_false, cut_mg_false)], "BP - RP", "$M_G$", invert_y=True)

            out("BP-RP Magnitude vs M_g = (G Mean Magnitude + 5 - 5 * log10( 1000 / parallax )). More detail.")
            self.plot_removed([(cut_bprp_true, cut_mg_true), (cut_bprp_false, cut_mg_false)], "BP - RP", "$M_G$", xlim=(0.1,3.5), ylim=(-1, 15), invert_y=True)
        except:
            out("An error occurred while plotting diagnostics. This usually occurs because no sources passed a photometric cut.")

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
            out("plotting cluster" + str(i))
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