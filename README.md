# f_esc

Code accompanying the paper:

> **Ionizing photon production and escape fractions during cosmic reionization in the TNG50 simulation**
> Ivan Kostyuk, Dylan Nelson, Benedetta Ciardi, Martin Glatzle, Annalisa Pillepich
> *MNRAS 521, 3077 (2023)* — [DOI: 10.1093/mnras/stad677](https://doi.org/10.1093/mnras/stad677) · [arXiv:2207.11278](https://arxiv.org/abs/2207.11278)

## Science context

Understanding which galaxies supply the bulk of ionizing photons during the epoch of reionization (6 < z < 10) requires knowing not just how many ionizing photons are produced, but what fraction escape from their host halos. Measuring this escape fraction *f*<sub>esc</sub> directly at high redshift is observationally out of reach, and theoretical predictions differ substantially depending on how radiative transfer and sub-grid absorption are modelled.

This repository contains the analysis code for a study that post-processes ~4 000 TNG50 galaxies (10<sup>6</sup> ≤ M<sub>★</sub>/M<sub>☉</sub> ≤ 10<sup>8</sup>) with the 3D multi-frequency radiative transfer code [CRASH](https://crash.cita.utoronto.ca/), computing *f*<sub>esc</sub> self-consistently along thousands of lines of sight per halo.

Key results:
- *f*<sub>esc</sub> rises with stellar mass from ~0.3 at M<sub>★</sub> = 10<sup>6</sup> M<sub>☉</sub> to ~0.6 at M<sub>★</sub> = 10<sup>7.5</sup> M<sub>☉</sub>, with hints of a turnover at higher masses
- Significant scatter at fixed mass, driven by diversity in the ionizing photon rate and the spatial relationship between stellar sources and the gas density field
- Dust reduces *f*<sub>esc</sub> by a few percent at low masses and up to 10% for M<sub>★</sub> ≳ 10<sup>6.5</sup> M<sub>☉</sub>
- *f*<sub>esc</sub> is energy-dependent: photons above 54.4 eV (relevant for He reionization and binary-star models) escape at a substantially lower rate than lower-energy ionizing photons
- Halos with M<sub>★</sub> ≲ 10<sup>7.5</sup> M<sub>☉</sub> dominate the global ionizing emissivity at all redshifts studied

## Repository structure

The notebooks in this repo are **load-bearing**: each one produces specific figures in the paper and contains the full analysis pipeline for that result. They are not exploratory scratch pads — they are the analysis.

### Python scripts

These scripts interact directly with the CRASH post-processing runs or TNG50 snapshot data and must be run before the analysis notebooks.

| Script | Purpose |
|---|---|
| `build_df.py` | Reads CRASH output folders and assembles a pandas dataframe of halo properties |
| `update_df.py` | Adds stellar ages, gas clumping, and surface density columns to the dataframe |
| `add_metallicities.py` | Appends stellar metallicities to the dataframe |
| `prepare_new_runs.py` | Tools for setting up a new CRASH post-processing run |
| `synchronize_folders.py` | Moves halos from reduced-source RT runs to full-source *f*<sub>esc</sub> calculations |
| `phase_diagram.py` | Generates phase diagram data |

### Analysis notebooks (paper figures)

These notebooks read the compiled dataframe and produce the paper figures. They can be run independently of the simulation data, provided the dataframe is available (see [Data](#data) below).

| Notebook | Figures | What it does |
|---|---|---|
| `counts_histogram.ipynb` | Fig. 1 | Halo population summary plots |
| `median_fesc.ipynb` | Figs. 2, 3, 13 | *f*<sub>esc</sub> vs stellar/halo mass; dust comparison; literature comparison |
| `fesc_vs_quant.ipynb` | — | *f*<sub>esc</sub> as a function of various galaxy properties |
| `histograms.ipynb` | Fig. 6 | 2D histograms of halo properties vs *f*<sub>esc</sub>; bimodal escape fraction analysis |
| `loc_esc.ipynb` | Fig. 7 | Effect of different local (cloud-scale) escape fraction assumptions |
| `uv_emissivity.ipynb` | Figs. 8, 9 | Escaped ionizing photon density in TNG50; comparison to literature |
| `r_fesc.ipynb` | Fig. 10 | Properties of individual stellar particles and their contribution to *f*<sub>esc</sub> |
| `spectra.ipynb` | Figs. 11, 12 | Spectral (energy) dependence of *f*<sub>esc</sub> |
| `numerical_convergence_tests.ipynb` | Fig. A1 | Source number reduction convergence tests |
| `large_radii.ipynb` | — | *f*<sub>esc</sub> as a function of aperture radius (not in final paper) |

### Notebooks requiring full simulation data

These two notebooks need the raw CRASH density and ionization maps, not just the dataframe:

| Notebook | Figures | What it does |
|---|---|---|
| `esc_fraction.ipynb` | Fig. 5 | Line-of-sight *f*<sub>esc</sub> maps and density projections for individual halos |
| `halo_image.ipynb` | Fig. 4 | Projected images of halo gas and stellar distributions |

### Simulation management notebooks

Not used for data analysis — for preparing and cleaning CRASH runs:

| Notebook | Purpose |
|---|---|
| `merge_sources.ipynb` | Source reduction pre-processing; includes validation tests |
| `clean_up.ipynb` | Deletes intermediate simulation files to reduce disk usage |

## Data

The analysis notebooks require the compiled halo dataframe, which is not included in this repository due to size. To obtain the dataframe, please open an issue or reach out via the contact on the [arXiv page](https://arxiv.org/abs/2207.11278).

## Requirements

Python ≥ 3.7, NumPy, pandas, matplotlib, h5py, astropy. Access to IllustrisTNG simulation data and CRASH outputs is required to run the simulation-facing scripts.

## Citation

If you use this code, please cite:

```bibtex
@article{Kostyuk2023fesc,
  author  = {Kostyuk, Ivan and Nelson, Dylan and Ciardi, Benedetta and Glatzle, Martin and Pillepich, Annalisa},
  title   = {Ionizing photon production and escape fractions during cosmic reionization in the {TNG50} simulation},
  journal = {MNRAS},
  year    = {2023},
  volume  = {521},
  pages   = {3077},
  doi     = {10.1093/mnras/stad677},
  eprint  = {2207.11278},
  archivePrefix = {arXiv}
}
```
