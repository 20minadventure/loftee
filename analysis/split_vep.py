import random
from datetime import datetime

import hail as hl
import dxpy
from analysis.utils.load_spark import hl_init, SC
from analysis.utils.dxpathlib import PathDx


def main():
    db_ref = 'wgs_mt'
    tmp_path = PathDx(database=db_ref)

    db_hail_tmp = 'hail_tmp'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_hail_tmp} LOCATION 'dnax://'")
    hail_tmp_path = PathDx(database='hail_tmp')

    log_path = f'/tmp/{datetime.now().strftime("%Y%m%d-%H%M")}-{random.randrange(16 ** 6):04x}.log'
    vep_config_path = PathDx('analysis/utils/vep-config.json')


    hl_init(tmp_dir=tmp_path.rstr, log=log_path)
