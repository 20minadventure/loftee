import hail as hl


class VCFFilter:
    def __init__(self):
        pass

    def _variant_qc(filter_fun):
        def demand_variant_qc(self, mt, *args, qc_col_name='variant_qc', **kwargs):
            if not hasattr(mt, qc_col_name):
                raise ValueError(f"MatrixTable has no attribute '{qc_col_name}'")
            return filter_fun(self, mt, *args, **kwargs)
        return demand_variant_qc

    def _split_multi(filter_fun):
        def demand_split(self, mt, *args, **kwargs):
            split_multi_col_name = 'was_split'
            if not hasattr(mt, split_multi_col_name):
                raise ValueError(
                    f"MatrixTable has no attribute '{split_multi_col_name}'"
                )
            return filter_fun(self, mt, *args, **kwargs)
        return demand_split

    def mean_read_depth(self, mt, min_depth=7):
        mt = mt.filter_rows(
            hl.agg.mean(mt.DP) >= min_depth
        )
        return mt

    @_variant_qc
    def variant_missingness(self, mt, min_ratio=0.1):
        return mt.filter_rows(mt.variant_qc.call_rate >= min_ratio)

    @_split_multi
    def allele_balance(self, mt, n_sample=1, sample_ratio=None, min_ratio=0.15):
        mt = mt.filter_rows(
            hl.agg.any(
                mt.GT.is_het()
                & (hl.min(mt.AD) / hl.sum(mt.AD) >= min_ratio)
            )
            | hl.agg.all(~mt.GT.is_het())
        )
        return mt

    def is_indel(self):
        pass

    @_variant_qc
    @_split_multi
    def hardy_weinberg(self, mt, min_p_value=1e-15):
        return mt.filter_rows(mt.variant_qc.p_value_hwe >= min_p_value)
