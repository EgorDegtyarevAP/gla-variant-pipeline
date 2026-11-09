PY=python
all: stage1 stage2 stage3 stage4 stage5 stage6 stage6b stage7
stage1:
	$(PY) scripts/01_ingest_clinvar.py
	$(PY) scripts/02_ingest_gnomad.py
	$(PY) scripts/03_ingest_amenability.py
	$(PY) scripts/04_ingest_literature.py
stage2:
	$(PY) scripts/05_harmonize.py
	$(PY) -m pytest -q tests
stage3:
	$(PY) scripts/06_classify_cpg.py
stage4:
	$(PY) scripts/07_structure.py
	$(PY) scripts/08_ddg.py
stage5:
	$(PY) scripts/09_conservation.py
stage6:
	$(PY) scripts/10_stats.py
stage6b:
	$(PY) scripts/13_multivariable.py
	$(PY) scripts/14_active_site_sensitivity.py
stage7:
	$(PY) scripts/11_figures.py
	$(PY) scripts/12_pymol.py
test:
	$(PY) -m pytest -q tests
