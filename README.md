# neutron-filter-atlas

A database-wide atlas of single-crystal thermal-neutron filters: for 18,449 inorganic crystals of the Materials Project, the filter thickness that passes 80 % of a thermal beam and the fast-neutron and photon rejection at that thickness, with a model uncertainty, at 293 K and 77 K.

It is for instrument scientists and engineers who choose or compare filter crystals, and for anyone who wants to recompute or extend the atlas.

```
micromamba create -f environment.yml && micromamba activate neutron-filter-atlas && pip install -e .
python scripts/check_data.py --class derived    # after unzipping the results archive (wiki/data.md) into the repository root
python -c "import pandas as pd; a = pd.read_parquet('data/derived/atlas.parquet'); print(a[(a.formula == 'Al2O3') & (a.sg_symbol == 'R-3c')][['T', 'L08', 'A_fast_1E', 'A_gamma_1MeV']])"
```

On Windows, or without conda, install from `requirements-lock.txt` instead ([wiki/installation.md](wiki/installation.md)). The data release is available from the authors on request ([wiki/data.md](wiki/data.md#the-data-release)).

Full documentation: how to read the atlas, run the pipeline and extend it: [wiki/README.md](wiki/README.md).

Cite: see [CITATION.cff](CITATION.cff) and [wiki/citing.md](wiki/citing.md) for the upstream sources to cite alongside.

Licence: code MIT ([LICENSE](LICENSE)); data: see [wiki/data.md](wiki/data.md#licences).
