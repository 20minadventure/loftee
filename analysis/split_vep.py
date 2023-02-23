import re
import sys
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

    if len(sys.argv) > 1:
        chrs = sys.argv[1:]
    else:
        chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']

    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    for p in b_vcf.listdir():
        m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
        if m:
            contig = m.group(1)
            block = m.group(2)
            chr_b_path = tmp_path / f'chr-{contig}-b{block}.mt'
            if contig in chrs and chr_b_path not in tmp_path.listdir():
                print(chr_b_path)
                try:
                    mt = hl.import_vcf(
                        p.rstr,
                        force_bgz=True,
                        array_elements_required=False,
                    )
                    mt = hl.split_multi_hts(mt)
                    mt = hl.vep(mt, vep_config_path.rstr)

                    mt.write(chr_b_path.rstr)
                except Exception as e:
                    print('ERROR: ', p)
                    print(e)
