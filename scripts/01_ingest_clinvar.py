"""Stage 1. IN: data/raw/clinvar_GLA_final.csv, clinvar_unparsed_names.csv.
OUT: data/interim/clinvar.parquet, clinvar_unparsed.parquet.
NOTE: the input is the user's pre-filtered ClinVar GLA export (1990 rows), not the full variant_summary;
download date/ClinVar release are not recorded in the file -> fill MANIFEST version by hand."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from gla.common import RAW, RAW_FILES, INTERIM, log_source
from gla.ingest import read_clinvar, read_clinvar_unparsed

if __name__ == "__main__":
    INTERIM.mkdir(parents=True, exist_ok=True)
    for k in ("clinvar", "clinvar_unparsed"):
        log_source(RAW / RAW_FILES[k], "ClinVar (user export)", "release/date NOT RECORDED - to be filled")
    d = read_clinvar(); d.to_parquet(INTERIM / "clinvar.parquet")
    u = read_clinvar_unparsed(); u.to_parquet(INTERIM / "clinvar_unparsed.parquet")
    print(f"clinvar: {len(d)} rows, hgvs_c parsed {d.hgvs_c.notna().sum()}; unparsed names: {len(u)}")
