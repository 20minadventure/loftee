import re
import sys
import random
from datetime import datetime

import hail as hl
import dxpy
from analysis.utils.load_spark import hl_init, SC
from analysis.utils.dxpathlib import PathDx
from analysis.utils.variant_filtering import VCFFilter

import pkg_resources
file_path = pkg_resources.resource_filename('analysis', 'utils/vep-config.json')

def split_annotate(p, out, permit_shuffle=False, vep_config_path=file_path):
    mt = hl.import_vcf(
        p.rstr,
        force_bgz=True,
        array_elements_required=False,
    )
    mt = hl.split_multi_hts(mt, permit_shuffle=permit_shuffle)
    mt = hl.vep(mt, vep_config_path)

    mt.write(out.rstr, overwrite=True)

def main():
    db_ref = 'wgs_mt'
    tmp_path = PathDx(database=db_ref)

    db_hail_tmp = 'hail_tmp'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_hail_tmp} LOCATION 'dnax://'")
    hail_tmp_path = PathDx(database='hail_tmp')

    log_path = f'/tmp/{datetime.now().strftime("%Y%m%d-%H%M")}-{random.randrange(16 ** 6):04x}.log'

    hl_init(tmp_dir=tmp_path.rstr, log=log_path)

    if len(sys.argv) > 1:
        chrs = sys.argv[1:]
    else:
        chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']

    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    out_mts = []
    for p in b_vcf.listdir():
        m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
        if m:
            contig = m.group(1)
            block = m.group(2)
            chr_b_path = tmp_path / f'chr-{contig}-b{block}.mt'
            if contig in chrs and block == '10':
                print(chr_b_path)
                if  chr_b_path in tmp_path.listdir():
                    try:
                        _mt = hl.read_matrix_table(chr_b_path.rstr)
                    except Exception as e:
                        print('Saved matrix table corrupted: rerunning with permit_shuffle=True')
                        split_annotate(p, chr_b_path, permit_shuffle=True)
                else:
                    try:
                        split_annotate(p, chr_b_path)
                    except Exception as e:
                        print('ERROR: ', p)
                        print(e)
                        print('Rerunning with permit_shuffle=True')
                        try:
                            split_annotate(p, chr_b_path, permit_shuffle=True)
                        except Exception as x:
                            print(x)
                            continue
                out_mts.append(chr_b_path)
    # out table
    mt_lof = hl.MatrixTable.union_rows(
        *(hl.read_matrix_table(b.rstr) for b in out_mts)
    )

    CANONICAL = 1
    mt_lof = mt_lof.explode_rows(
        mt_lof.vep.transcript_consequences
    )
    mt_lof = mt_lof.filter_rows(
        mt_lof.vep.transcript_consequences.canonical == CANONICAL
    )
    mt_lof = mt_lof.filter_rows(
        hl.is_defined(mt_lof.vep.transcript_consequences.lof)
    )

    # Filter VCF
    mt_filter = VCFFilter()
    mt_lof = mt_filter.mean_read_depth(mt_lof, min_depth=7)
    mt_lof = hl.variant_qc(mt_lof)
    mt_lof = mt_filter.variant_missingness(mt_lof, min_ratio=0.1)
    mt_lof = mt_filter.hardy_weinberg(mt_lof, min_p_value=1e-15)
    mt_lof = mt_filter.allele_balance(mt_lof, n_sample=1, min_ratio=0.15)

    mt_lof_grouped = mt_lof.group_rows_by(mt_lof.vep.transcript_consequences.gene_id, mt_lof.vep.transcript_consequences.gene_symbol)
    result = mt_lof_grouped.aggregate_entries(
        hc_lof_hom=hl.agg.max(mt_lof.GT.n_alt_alleles() * (mt_lof.vep.transcript_consequences.lof == 'HC')) == 2,
        hc_lof_n_het=hl.agg.sum(mt_lof.GT.is_het() * (mt_lof.vep.transcript_consequences.lof == 'HC')),
        any_lof_hom=hl.agg.max(mt_lof.GT.n_alt_alleles()) == 2,
        any_lof_n_het=hl.agg.sum(mt_lof.GT.is_het()),
        n_non_ref=hl.agg.sum(mt_lof.GT.is_non_ref()),
    ).result()
    
    result = result.checkpoint((tmp_path / 'result').rstr, overwrite=True)
    df = result.entries().to_pandas()
    df[['gene_symbol', 's', 'hc_lof_hom', 'hc_lof_n_het']].pivot(index='s', columns='gene_symbol', values='hc_lof_n_het').to_csv('file:/opt/notebooks/out.csv')
