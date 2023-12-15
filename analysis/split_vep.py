import re
import sys
import random
import gzip
from math import ceil
from itertools import chain
from datetime import datetime

import hail as hl
import dxpy
import numpy as np
from analysis.utils.load_spark import hl_init, SC
from analysis.utils.dxpathlib import PathDx
from analysis.utils.variant_filtering import VCFFilter

import pkg_resources
file_path = pkg_resources.resource_filename('analysis', 'utils/vep-config.json')

def split_annotate(p, out, permit_shuffle=False, vep_config_path=PathDx(file_path)):
    mt = hl.import_vcf(
        p.rstr,
        force_bgz=True,
        array_elements_required=False,
    )
    mt = hl.split_multi_hts(mt, permit_shuffle=permit_shuffle)
    mt = hl.vep(mt, vep_config_path.rstr)

    mt.write(out.rstr, overwrite=True)


def reannotate_vep(mt_path, out_path, vep_config_path=PathDx(file_path)):
    mt = hl.read_matrix_table(mt_path.rstr)
    mt = mt.drop('vep', 'vep_proc_id')
    mt = hl.vep(mt, vep_config_path.rstr)

    mt.write(out_path.rstr, overwrite=True)


def main():
    db_ref = 'wes_mt_bl'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_ref} LOCATION 'dnax://'")
    
    tmp_path = PathDx(database=db_ref)

    db_hail_tmp = 'hail_tmp'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_hail_tmp} LOCATION 'dnax://'")
    hail_tmp_path = PathDx(database=db_hail_tmp)

    log_path = f'/tmp/{datetime.now().strftime("%Y%m%d-%H%M")}-{random.randrange(16 ** 6):04x}.log'

    hl_init(tmp_dir=hail_tmp_path.rstr, log=log_path)

    chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']
    blocks = [str(i) for i in range(100)]
    eids = None
    if len(sys.argv) > 1 and sys.argv[1] in chrs:
        chrs = [sys.argv[1]]
    if len(sys.argv) > 2:
        eids_path = PathDx('/mnt/project/') / sys.argv[2]
        with open(eids_path) as f:
            eids = [eid.rstrip() for eid in f.readlines()]
        

    wes_mt = PathDx(database='wes_mt')
    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    for p in b_vcf.listdir():
        m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
        if m:
            contig = m.group(1)
            block = m.group(2)
            chr_b_path = tmp_path / f'chr-{contig}-b{block}.mt'
            chr_b_path_out = wes_mt / f'chr-{contig}-b{block}.mt'
            if contig in chrs and block in blocks:
                print(chr_b_path_out, flush=True)
                reannotate_vep(chr_b_path, chr_b_path_out)
