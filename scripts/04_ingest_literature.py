"""Stage 1 (optional). No literature variant table has been supplied yet.
Expected input: data/raw/literature_variants.csv with columns hgvs_c, source_id (PMID/DOI), cohort.
Until it exists the master table is built from ClinVar + gnomAD + Amicus only (stated in audit_report.md)."""
print("literature table not provided - skipped")
