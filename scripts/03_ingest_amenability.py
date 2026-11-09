"""Stage 1. IN: data/raw/galafold_amenability.csv (+ amicus_excluded_multilocus.csv).
OUT: data/interim/amenability.parquet (one row per nucleotide change), amenability_unresolved.parquet."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from gla.common import RAW, RAW_FILES, INTERIM, log_source
from gla.ingest import read_amicus

if __name__ == "__main__":
    INTERIM.mkdir(parents=True, exist_ok=True)
    for k in ("amicus", "amicus_excluded"):
        log_source(RAW / RAW_FILES[k], "Amicus Galafold amenability table", "version/date NOT RECORDED - to be filled")
    ok, bad, n_main, n_excl = read_amicus()
    ok.to_parquet(INTERIM / "amenability.parquet"); bad.to_parquet(INTERIM / "amenability_unresolved.parquet")
    print(f"amicus rows {n_main}, excluded-list {n_excl}; resolved nt changes {len(ok)}; unresolved {len(bad)}")
    print(bad.reason.value_counts().to_string())
