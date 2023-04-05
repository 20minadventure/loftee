import re
import sys
import random
from datetime import datetime

import hail as hl
import gzip
import dxpy
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

def main():
    db_ref = 'wes_mt'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_ref} LOCATION 'dnax://'")
    
#     tmp_path = PathDx(database=db_ref)
    tmp_path = PathDx(database_id='database-GPffV5jJ3B1y1JKpGvJ5b7pZ')

    db_hail_tmp = 'hail_tmp'
    SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_hail_tmp} LOCATION 'dnax://'")
    hail_tmp_path = PathDx(database='hail_tmp')

    log_path = f'/tmp/{datetime.now().strftime("%Y%m%d-%H%M")}-{random.randrange(16 ** 6):04x}.log'

    hl_init(tmp_dir=hail_tmp_path.rstr, log=log_path)

    chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']
    blocks = [str(i) for i in range(100)]
    if len(sys.argv) > 1:
        chrs = [sys.argv[1]]
    if len(sys.argv) > 2:
        blocks = sys.argv[2].split(',')
        

    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    out_mts = []
    for p in b_vcf.listdir():
        m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
        if m:
            contig = m.group(1)
            block = m.group(2)
            chr_b_path = tmp_path / f'chr-{contig}-b{block}.mt'
            if contig in chrs and block in blocks:
                print(chr_b_path, flush=True)
                try:
                    tmp_paths_list = tmp_path.listdir()
                except Exception as e:
                    tmp_paths_list = []
                if  chr_b_path in tmp_paths_list:
                    try:
                        _mt = hl.read_matrix_table(chr_b_path.rstr)
                    except Exception as e:
                        print('Saved matrix table corrupted: rerunning with permit_shuffle=True', flush=True)
                        split_annotate(p, chr_b_path, permit_shuffle=True)
                else:
                    try:
                        split_annotate(p, chr_b_path)
                    except Exception as e:
                        print('ERROR: ', p)
                        print(e)
                        print('Rerunning with permit_shuffle=True', flush=True)
                        try:
                            split_annotate(p, chr_b_path, permit_shuffle=True)
                        except Exception as x:
                            print(x, flush=True)
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
        hc_lof_hom_n=hl.agg.max(mt_lof.GT.n_alt_alleles() * (mt_lof.vep.transcript_consequences.lof == 'HC')),
        hc_lof_n_het=hl.agg.sum(mt_lof.GT.is_het() * (mt_lof.vep.transcript_consequences.lof == 'HC')),
        any_lof_hom_n=hl.agg.max(mt_lof.GT.n_alt_alleles()),
        any_lof_n_het=hl.agg.sum(mt_lof.GT.is_het()),
        n_non_ref=hl.agg.sum(mt_lof.GT.is_non_ref()),
    ).result()
    result = result.transmute_entries(
        hc_lof_hom=hl.if_else(result.hc_lof_hom_n == 2, True, False, missing_false=True),
        any_lof_hom=hl.if_else(result.any_lof_hom_n == 2, True, False, missing_false=True),
    )

    result = result.checkpoint((hail_tmp_path / 'result').rstr, overwrite=True)
    result = result.annotate_entries(
        value=hl.if_else(
            result.hc_lof_hom,
            2,
            hl.min(result.hc_lof_n_het, 2)
        )
    )

    result_bm_path = hail_tmp_path / 'result.bm'
    block_size = 512
    print('Save as block matrix', flush=True)
    hl.linalg.BlockMatrix.write_from_entry_expr(result.value, result_bm_path.rstr, block_size=block_size, overwrite=True)
    arr = hl.linalg.BlockMatrix.read(result_bm_path.rstr)

    print('Export to csv', flush=True)
    patients = result.s.collect()
    gene_symbols = result.gene_symbol.collect()

    arr_t = arr.T
    blocks = 512 * 100
    out_path = f'/opt/notebooks/out-{"-".join(chrs)}-{random.randrange(16 ** 6):04x}.csv.gz'
    with gzip.open(out_path, 'wt') as f:
        assert len(patients) == arr_t.shape[0]
        assert len(gene_symbols) == arr_t.shape[1]
        f.write(f"gene_symbol,{','.join(gene_symbols)}\n")
        for b in range(ceil(arr_t.shape[0] / blocks)):
            start = b * blocks
            end = min((b + 1) * blocks, arr_t.shape[0])
            arr_np = arr_t[start:end, :].to_numpy().astype(np.int8)
            for i in range(arr_np.shape[0]):
                line = f"{patients[i]},{','.join(map(str, arr_np[i, :]))}\n"
                f.write(line)
