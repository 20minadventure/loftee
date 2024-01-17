import re
import sys
import random
import gzip
import subprocess
from math import ceil, floor
from itertools import chain
from functools import reduce
from datetime import datetime

import hail as hl
import dxpy
import numpy as np
from analysis.utils.load_spark import hl_init, SC
from analysis.utils.dxpathlib import PathDx
from analysis.utils.variant_filtering import VCFFilter

import pkg_resources
file_path = pkg_resources.resource_filename('analysis', 'utils/vep-config.json')


# Initialize hail
db_ref = 'wes_mt'
SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_ref} LOCATION 'dnax://'")

tmp_path = PathDx(database=db_ref)

db_hail_tmp = 'hail_tmp'
SC.sql(f"CREATE DATABASE IF NOT EXISTS {db_hail_tmp} LOCATION 'dnax://'")
hail_tmp_path = PathDx(database=db_hail_tmp)

log_path = f'/tmp/{datetime.now().strftime("%Y%m%d-%H%M")}-{random.randrange(16 ** 6):04x}.log'

hl_init(tmp_dir=hail_tmp_path.rstr, log=log_path)

# Get parameters
chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']
eids = None
if len(sys.argv) > 1:
    chrs = sys.argv[1].split(',')
if len(sys.argv) > 2:
    eids_path = PathDx('/mnt/project/') / sys.argv[2]
    with open(eids_path) as f:
        eids = [eid.rstrip() for eid in f.readlines()]


def try_to_int(x):
    try:
        return int(x)
    except Exception:
        return x


def nsort(name):
    sub_names = re.split(r'(\d+)', name)
    return [try_to_int(n) for n in sub_names]


def mt_name(contig, block):
    return f'chr-{contig}-b{block}.mt'


def match(a, b):
    b_dict = {x: i for i, x in enumerate(b)}
    return [b_dict.get(x, None) for x in a]


def split_list(n, k):
    """Iterate through slices spliting list of length n into k lists of equal
    size, eg:
    n = 5, k = 3
    return: (0, 2), (2, 4), (4, 5)"""
    avg_length, remainder = divmod(n, k)
    start = 0

    for i in range(k):
        end = start + avg_length + (1 if i < remainder else 0)
        yield start, end
        start = end


def split_annotate(p, out, permit_shuffle=False, vep_config_path=PathDx(file_path)):
    mt = hl.import_vcf(
        p.rstr,
        force_bgz=True,
        array_elements_required=False,
        block_size=128,
    )
    mt = hl.split_multi_hts(mt, permit_shuffle=permit_shuffle)
    mt = hl.vep(mt, vep_config_path.rstr)

    mt.write(out.rstr, overwrite=True)


def annotate_vcf():
    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    b_vcf_files = sorted(b_vcf.listdir(), key=lambda f: nsort(f.name))
    for p in b_vcf_files:
        m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
        if m:
            contig = m.group(1)
            block = m.group(2)
            chr_b_path = tmp_path / mt_name(contig, block)
            if contig in chrs:
                print(chr_b_path, flush=True)
                p_local = PathDx(f'/cluster/{p.name}')
                try:
                    tmp_paths_list = tmp_path.listdir()
                except Exception as e:
                    tmp_paths_list = []
                if  chr_b_path not in tmp_paths_list:
                    print(f'Copying {p_local.rstr}...', flush=True)
                    subprocess.run(['hdfs', 'dfs', '-cp', p.rstr, p_local.rstr], check=True)
                    p = p_local
                    try:
                        split_annotate(p, chr_b_path)
                    except Exception as e:
                        print('ERROR: ', p, flush=True)
                        print(e, flush=True)
                        print('Rerunning with permit_shuffle=True', flush=True)
                        try:
                            split_annotate(p, chr_b_path, permit_shuffle=True)
                        except Exception as x:
                            print(x, flush=True)
                            print('SECOND TRY FAILED, PLEASE CHECK FILE MANUALLY', flush=True)
                            continue
                    subprocess.run(
                        ['hdfs', 'dfs', '-rm', '-r', '-skipTrash', p_local.rstr],
                        check=True
                    )


def rare_variants_table():
    try:
        tmp_paths_list = tmp_path.listdir()
    except Exception as e:
        print(f'No VCF file is annotated', flush=True)
        return

    b_vcf = PathDx('/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release')
    b_vcf_files = sorted(b_vcf.listdir(), key=lambda f: nsort(f.name))
    for chrom in chrs:
        print(f'Chr {chrom}')
        out_mts = []
        for p in b_vcf_files:
            m = re.fullmatch(r'ukb23157_c(\d{1,2}|X|Y)_b(\d{1,3})_v1.vcf.gz', p.name)
            if m:
                contig = m.group(1)
                block = m.group(2)
                chr_b_path = tmp_path / mt_name(contig, block)
                if contig == chrom:
                    if chr_b_path in tmp_paths_list:
                        try:
                            has_success = any(
                                '_SUCCESS' in file.name
                                for file in chr_b_path.listdir()
                            )
                        except Exception as e:
                            has_success = False
                        if has_success:
                            print(f'{chr_b_path} OK', flush=True)
                            out_mts.append(chr_b_path)
                        else:
                            print(f'{chr_b_path} FAIL (no _SUCCESS)', flush=True)
                            out_mts.append(False)
                    else:
                        print(f'{chr_b_path} FAIL (no file)', flush=True)
                        out_mts.append(False)
        if all(out_mts):
            _chr_table(chrom, out_mts, eids)
        else:
            print(f'Some VCF files are not ready', flush=True)


def _chr_table(chrom, mts, eids):
    print('Unifying colnames...', flush=True)
    mts_dict = {b: hl.read_matrix_table(b.rstr) for b in mts}
    mts_patients = {b: mt.s.collect() for b, mt in mts_dict.items()}
    common_pats = reduce(lambda x, y: x & y, (set(p) for p in mts_patients.values()))
    if eids:
        common_pats &= set(eids)
    mts_unified = []
    for b, pats in mts_patients.items():
        pat_indices = match(list(common_pats), pats)
        mts_unified.append(mts_dict[b].choose_cols(pat_indices))

    # out table
    min_batch = 19
    n, k = len(mts), max(1, floor(len(mts) / min_batch))
    all_gene_names = set()
    for i, (start, end) in enumerate(split_list(n, k)):
        print(f'Part {i}: [{start}:{end}]', flush=True)
        mt_lof = hl.MatrixTable.union_rows(*mts_unified[start:end])

        CANONICAL = 1
        mt_lof = mt_lof.explode_rows(
            mt_lof.vep.transcript_consequences
        )
        mt_lof = mt_lof.filter_rows(
            (mt_lof.vep.transcript_consequences.canonical == CANONICAL)
            & (mt_lof.vep.transcript_consequences.biotype == 'protein_coding')
        )
        mt_lof = mt_lof.annotate_rows(
            gene_name=hl.if_else(
                hl.is_defined(mt_lof.vep.transcript_consequences.gene_symbol),
                mt_lof.vep.transcript_consequences.gene_symbol,
                mt_lof.vep.transcript_consequences.gene_id
            )
        )
        print('....aggregating names', flush=True)
        all_gene_names |= mt_lof.aggregate_rows(hl.agg.collect_as_set(mt_lof.gene_name))

        print('....filtering', flush=True)
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
        mt_lof = mt_lof.filter_rows(~mt_lof.was_split)
        mt_lof.write(PathDx(f'/cluster/result-{chrom}-0-p{i}').rstr, overwrite=True)

    print('Unioning all', flush=True)
    mt_lof = hl.MatrixTable.union_rows(
        *[
            hl.read_matrix_table(PathDx(f'/cluster/result-{chrom}-0-p{i}').rstr)
            for i in range(k)
        ]
    )
    mt_lof = mt_lof.checkpoint(PathDx(f'/cluster/result-{chrom}-0').rstr, overwrite=True)

    mt_lof_grouped = mt_lof.group_rows_by(mt_lof.gene_name)
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

    result = result.annotate_entries(
        value=hl.if_else(
            result.hc_lof_hom,
            2,
            hl.min(result.hc_lof_n_het, 2)
        )
    )
    result = result.checkpoint(PathDx(f'/cluster/result-{chrom}-1b').rstr, overwrite=True)

    result_bm_path = PathDx(f'/cluster/result-{chrom}.bm')
    block_size = 512
    print('Save as block matrix', flush=True)
    hl.linalg.BlockMatrix.write_from_entry_expr(result.value, result_bm_path.rstr, block_size=block_size, overwrite=True)
    arr = hl.linalg.BlockMatrix.read(result_bm_path.rstr)

    print('Export to csv', flush=True)
    patients = result.s.collect()
    gene_names = result.gene_name.collect()
    zero_genes = all_gene_names - set(gene_names)

    arr_t = arr.T
    blocks = 512 * 100
    out_path = f'/opt/notebooks/out-{chrom}-{random.randrange(16 ** 6):04x}.csv.gz'
    with gzip.open(out_path, 'wt') as f:
        assert len(patients) == arr_t.shape[0]
        assert len(gene_names) == arr_t.shape[1]
        f.write(f"s,{','.join(chain(gene_names, zero_genes))}\n")
        for b in range(ceil(arr_t.shape[0] / blocks)):
            start = b * blocks
            end = min((b + 1) * blocks, arr_t.shape[0])
            arr_np = arr_t[start:end, :].to_numpy().astype(np.int8)
            for i in range(arr_np.shape[0]):
                line = f"{patients[i + start]},{','.join(chain(map(str, arr_np[i, :]), '0' * len(zero_genes)))}\n"
                f.write(line)
