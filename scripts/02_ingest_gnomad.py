"""Stage 1. IN: data/raw/gnomad_v4.1.1_ENSG00000102393_2026_09_11.csv. OUT: data/interim/gnomad.parquet.
Keeps the canonical transcript only (ENST00000218516.4); other transcripts are counted in the audit."""
import sys, pathlib, json; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from gla.common import RAW, RAW_FILES, INTERIM, log_source
from gla.ingest import read_gnomad

if __name__ == "__main__":
    INTERIM.mkdir(parents=True, exist_ok=True)
    log_source(RAW / RAW_FILES["gnomad"], "gnomAD", "v4.1.1 (export 2026-09-11)", "https://gnomad.broadinstitute.org")
    d, meta = read_gnomad(); d.to_parquet(INTERIM / "gnomad.parquet")
    (INTERIM / "gnomad_meta.json").write_text(json.dumps(meta, indent=1))
    print(meta, "hgvs_c parsed", int(d.hgvs_c.notna().sum()))
