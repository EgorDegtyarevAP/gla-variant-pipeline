# Samples (fixed before analysis; denominators always reported)
- S_all: every row of variants_master
- S_mis: variant_class == missense (mnv_missense is kept out)
- S_amen: S_mis with amenability_status in {amenable, non_amenable}  (criterion: Amicus table, see amenability_criterion)
- S_path: S_mis with clinvar_class in {P, LP, P/LP}  (conflicting and VUS excluded)
