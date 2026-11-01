# Citing

What this page is for: how to cite the software and the data release, and the upstream sources to cite with them.
Who should read it: anyone who publishes work that uses the atlas or the pipeline.

## The software and the data release

The citation metadata are in `CITATION.cff` at the repository root (GitHub shows them under "Cite this repository"):

> Kaukabi, M. A. & Gilbran, A. (2026). *neutron-filter-atlas: a database-wide atlas of single-crystal neutron
> filters*, version 0.5.0. https://github.com/arsyakaukabi/neutron-filter-atlas

Cite the data release as the data products of the same version:

> Kaukabi, M. A. & Gilbran, A. (2026). *neutron-filter-atlas: pipeline, data products and documentation*,
> version 0.5.0 (data release).

The data release will receive a persistent identifier when it is archived openly; until then it is available from
the authors on request ([Data](data.md)). Name the version you used, and for single values the `material_id`,
temperature and column ([FAQ](faq.md#how-do-i-cite-one-entry)).

## Upstream sources

The atlas is built on the following data and codes. Cite those your use depends on; for any use of atlas values,
that is at least the Materials Project, NCrystal and ENDF/B-VIII.0.

| Source | What it provides | Reference |
|---|---|---|
| Materials Project | candidate structures, energies, elastic moduli | Jain, A. *et al.* (2013). *APL Mater.* **1**, 011002. https://doi.org/10.1063/1.4812323 |
| Materials Project phonons, DFPT | phonon DOS and force constants (`dfpt`) | Petretto, G. *et al.* (2018). *Sci. Data* **5**, 180065. https://doi.org/10.1038/sdata.2018.65 |
| Materials Project phonons, finite displacements | phonon DOS and force constants (`pheasy`) | Sahasrabuddhe, H. *et al.* (2026). ChemRxiv preprint, version 1. https://doi.org/10.26434/chemrxiv.15004632/v1 |
| NCrystal | thermal-neutron cross sections | Cai, X.-X. & Kittelmann, T. (2020). *Comput. Phys. Commun.* **246**, 106851. https://doi.org/10.1016/j.cpc.2019.07.015 ; Kittelmann, T. & Cai, X.-X. (2021). *Comput. Phys. Commun.* **267**, 108082. https://doi.org/10.1016/j.cpc.2021.108082 |
| ENDF/B-VIII.0 | fast-neutron and photon cross sections | Brown, D. A. *et al.* (2018). *Nucl. Data Sheets* **148**, 1–142. https://doi.org/10.1016/j.nds.2018.02.001 |
| EXFOR | measured total cross-sections of the validation | Otuka, N. *et al.* (2014). *Nucl. Data Sheets* **120**, 272–276. https://doi.org/10.1016/j.nds.2014.07.065 ; and the original measurements (authors per dataset in `stage2_exfor_benchmark.parquet`) |
| phonopy | finite-displacement phonons, force-constant screens | Togo, A. (2023). *J. Phys. Soc. Jpn.* **92**, 012001. https://doi.org/10.7566/JPSJ.92.012001 ; Togo, A., Chaput, L., Tadano, T. & Tanaka, I. (2023). *J. Phys. Condens. Matter* **35**, 353001. https://doi.org/10.1088/1361-648X/acd831 |
| SevenNet (SevenNet-MF-ompa) | learned-potential phonons of the atlas | Kim, J. *et al.* (2025). *J. Am. Chem. Soc.* **147**, 1042–1054. https://doi.org/10.1021/jacs.4c14455 |
| OpenMC | slab transport, nuclear-data extraction | Romano, P. K. *et al.* (2015). *Ann. Nucl. Energy* **82**, 90–97. https://doi.org/10.1016/j.anucene.2014.07.048 |
| OpenMC depletion | activation of the shortlist | Romano, P. K., Josey, C. J., Johnson, A. E. & Liang, J. (2021). *Ann. Nucl. Energy* **152**, 107989. https://doi.org/10.1016/j.anucene.2020.107989 |
| pymatgen | structures and compositions | Ong, S. P. *et al.* (2013). *Comput. Mater. Sci.* **68**, 314–319. https://doi.org/10.1016/j.commatsci.2012.10.028 |
| MACE-MP-0 (benchmark only) | learned-potential benchmark | Batatia, I. *et al.* (2025). *J. Chem. Phys.* **163**, 184110. https://doi.org/10.1063/5.0297006 |
| ORB v3 (benchmark only) | learned-potential benchmark | Rhodes, B. *et al.* (2025). arXiv:2504.06231. https://doi.org/10.48550/arXiv.2504.06231 |

The Materials Project data are under CC BY 4.0, which requires attribution; EXFOR asks that the original
measurements be cited.
