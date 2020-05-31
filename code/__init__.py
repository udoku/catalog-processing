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