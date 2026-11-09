# variants_master schema
One row = one unique variant, key = `hgvs_c` on NM_000169 (syntactically normalised, see src/gla/hgvs.py).
Indels are NOT 3'-shifted (no reference offline): ambiguous indel equivalence across sources is possible and is
reported, not resolved. SNV keys are exact.

| column group | columns | note |
|---|---|---|
| identity | variant_key, hgvs_c, region, variant_class, protein, sources, n_sources, in_clinvar/in_gnomad/in_amicus | `protein` priority ClinVar > gnomAD > Amicus (one-letter, `*` = stop) |
| ClinVar | clinvar_id, clinvar_significance, clinvar_class (P, LP, P/LP, VUS, LB, B, B/LB, conflicting, not_classified), clinvar_review, clinvar_stars, clinvar_phenotype, clinvar_last_evaluated, rsid, clinvar_origin, clinvar_hgvs_source | |
| gnomAD (canonical ENST00000218516.4) | gnomad_id, g_pos/g_ref/g_alt (GRCh38, plus strand), gnomad_vep, gnomad_ac/an/af, gnomad_hom, gnomad_hemi, gnomad_pass (joint filter PASS), gnomad_flags, gnomad_faf(_group), gnomad_clinvar_id/class, cadd, revel, sift, polyphen, spliceai, pangolin, phylop | in-silico scores exist only for gnomAD rows |
| Amicus | amenability_status (`amenable` / `non_amenable` / NA), amenability_criterion, amicus_ids, amicus_variant_type, amicus_n_alt | **NA = not in the table; never non_amenable.** `amicus_n_alt` > 1: the table row listed several nucleotide changes for one protein change, the status is the row's |
| QC | protein_clinvar/gnomad/amicus, protein_cds, codon_no, codon_pos, ref_ok, protein_conflict | `ref_ok`: reference base equals NM_000169.3 (CDS, c.-22..-1, c.*1..*6; NA elsewhere). `protein_cds`: computed from the CDS for CDS SNVs and used as `protein`; reported proteins are only compared against it |
| samples | S_all, S_mis, S_amen, S_path | docs/samples.md |

variant_class values: upstream_flank / downstream_flank (outside c.-22..c.*6), ref_mismatch (reference base disagrees; excluded from samples), missense, nonsense, synonymous, start_lost, stop_lost, frameshift, inframe_indel, mnv_missense (delins with net 0 giving one amino-acid change), delins_unresolved, splice_site (intronic offset <=2), splice_spanning_indel, intronic, utr5, utr3, cds_snv_unannotated (exonic SNV with no protein in any source; resolved at stage 3).
Indel classes use the nucleotide length rule, not the protein string.
