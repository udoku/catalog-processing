import os
from datetime import datetime
import logging

current_time = datetime.now().strftime('%Y-%m-%d-%H-%M')
path = current_time

os.mkdir(path)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

catalog_logfile = path + "/logs.txt"
fh = logging.FileHandler(catalog_logfile)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

def out(message):
    print(message)
    logger.info(message)