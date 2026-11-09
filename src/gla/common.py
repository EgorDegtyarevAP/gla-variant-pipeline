"""Shared helpers: config, paths, provenance manifest."""
import hashlib, json, datetime, pathlib
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[2]
RAW, INTERIM, PROCESSED = (ROOT / "data" / p for p in ("raw", "interim", "processed"))
RESULTS = ROOT / "results"

def load_config():
    return yaml.safe_load((ROOT / "config.yaml").read_text())

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def log_source(path, source_name, version_or_date, url=""):
    """Append a raw-file record (sha256, version, URL) to data/raw/MANIFEST.jsonl."""
    mf = RAW / "MANIFEST.jsonl"
    digest = sha256(path)
    if mf.exists() and any(json.loads(l)["sha256"] == digest for l in mf.read_text().splitlines() if l.strip()):
        return None          # already registered (idempotent)
    rec = {"file": str(pathlib.Path(path).name), "sha256": digest, "source": source_name,
           "version": version_or_date, "url": url,
           "registered": datetime.date.today().isoformat()}
    with open(mf, "a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


RAW_FILES = {
    "clinvar": "clinvar_GLA_final.csv",
    "clinvar_unparsed": "clinvar_unparsed_names.csv",
    "gnomad": "gnomad_v4.1.1_ENSG00000102393_2026_09_11.csv",
    "amicus": "galafold_amenability.csv",
    "amicus_excluded": "amicus_excluded_multilocus.csv",
}
