import os
from datetime import datetime
import logging

current_time = datetime.now().strftime('%Y-%m-%d-%H-%M')
path = current_time + "/"
dpath = path + "diagnostics/"
cpath = path + "clustering/"
tpath = path + "tables/"

for p in (path, dpath, cpath, tpath):
    if not os.path.exists(p):
        os.mkdir(p)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

catalog_logfile = path + "logs.txt"
lh = logging.FileHandler(catalog_logfile)
lh.setLevel(logging.DEBUG)
lh.setFormatter(formatter)
logger.addHandler(lh)

catalog_readme = path + "README"

def out(message):
    print(message)
    logger.info(message)

def info_out(message):
    print(message)
    logger.info(message)

    readme = open(catalog_readme,'a')
    readme.write(message)
    readme.write("\n")
    readme.close()

info_out(path)

def list_to_str(list,separator):
    result = ""

    for e in list:
        result += str(e) + separator

    result = result.strip(separator)

    return result

gaia_cols_list = ["designation", "ra", "dec", "ra_error", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag", "bp_rp", "bp_g", "g_rp", "phot_g_n_obs", "phot_g_mean_flux_over_error", "phot_g_mean_flux_error", "phot_g_mean_flux", "phot_bp_mean_mag", "phot_rp_mean_mag"]
gaia_cols = list_to_str(gaia_cols_list,",")

tmass_cols_list = ["designation", "ra", "dec","err_maj", "err_min", "j_m", "h_m", "k_m", "j_cmsig", "h_cmsig", "k_cmsig", "k_snr"]
tmass_cols = list_to_str(tmass_cols_list,",")

allwise_cols_list = ["designation", "ra", "dec", "sigra", "sigdec", "w1mpro", "w2mpro", "w3mpro", "w4mpro", "w1sigmpro", "w2sigmpro", "w3sigmpro", "w4sigmpro", "w4snr"]
allwise_cols = list_to_str(allwise_cols_list,",")