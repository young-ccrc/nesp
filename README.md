# nesp — README

This repository contains scripts, notebooks and helper code used to
identify extreme weather events, propose minimum WRF domains, and
generate WRF/WPS intermediate files from BARRA2 data.

Contents
- `step1_data_mask_daily_5min_rv_Kim.ipynb`, `step3_highest_10_yk.ipynb`,
  `step4_barra_aws_extreme_test_Kim.ipynb` — notebooks used for selection
  and ranking of extreme events (Step 1, 3, 4).
- `aws_map_shared.ipynb`, `auto_generate_minimum_domain_maps.py`,
  `minimum_domain_proposal.ipynb` — tools and notebooks to propose
  minimum WRF domains based on station/AWS locations.
- `wrf/` — WRF helper scripts and the NCL script to create WPS/WRF
  intermediate files from BARRA2 (`wrf/wrfint_barra_v4.15.ncl`).

Quick start
-----------

1) Selection of extreme events

- Run the Step 1 notebook(s) to prepare daily/5-min masked data.
- Run `step3_highest_10_yk.ipynb` to identify and rank the top-N extreme
  events for each site/region. This produces CSV summaries (see `figure/`
  for examples such as `Adelaide_top10_extreme_events_ranked_with_latlon.csv`).
- Use `step4_barra_aws_extreme_test_Kim.ipynb` to compare station (AWS)
  and BARRA fields for selected events.

Notes:
- Notebooks assume you have the required data available locally (AWS
  station files, BARRA2 NetCDF). Open the notebooks in JupyterLab and run
  cells top-to-bottom for reproducible results.

2) Proposed domain (minimum domain generation)

- The notebooks and script used for domain proposals are:
  - `aws_map_shared.ipynb` — visualisation helper for AWS/station locations
  - `auto_generate_minimum_domain_maps.py` — automated domain-generation
    script that attempts to create minimum bounding domains for selected
    events/stations
  - `minimum_domain_proposal.ipynb` — interactive notebook for reviewing
    and refining proposed domains

How to use:

- Prepare the event selection CSV (output from Step 1/3 above).
- Run the automated script to produce initial domain maps:

```bash
cd <repo-root>
python3 auto_generate_minimum_domain_maps.py
```

- Inspect and refine the proposals in `minimum_domain_proposal.ipynb`.

Assumptions / dependencies:
- Python 3.8+ with typical geospatial libraries (geopandas, shapely,
  rasterio, cartopy, matplotlib). The notebooks contain explicit
  import cells — install the packages shown there into your environment.

3) Create WRF intermediate files using BARRA2

### Python converter 

Script:
- `wrf/barra_to_int.py`

Helper modules (bundled here):
- `wrf/WPSUtils.py`, `wrf/fortran_io.py`

What it does:
- Robust file/variable matching across BARRA2 naming variants
- Reads pressure-level fields (TT, UU, VV, SPECHUMD, GHT) and single-level
  PMSL and PSFC; optionally 2-m/10-m fields
- Writes correct WPS INTERMEDIATE records with LatLon projection metadata

Quick usage (example):

```bash
export BARRA2_ROOT="/g/data/w28/chs548/BARRA2_For_WRF"

python wrf/barra_to_int.py \
  --root "$BARRA2_ROOT" \
  --vars TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC \
  --outdir /path/to/wps/int \
  --prefix BARRA2 \
  --alias-hh \
  --verbose \
  2016-01-28_00 2016-02-01_00 6
```

Arguments of interest:
- `--root`: BARRA2 root directory
- `--outdir`: destination for INTERMEDIATE files
- `--prefix`: prefix used by metgrid (e.g., BARRA2)
- `--vars`: comma-separated selection; defaults to `TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC`
- `start end interval`: start/end (inclusive) as `YYYY-MM-DD_HH` and step in hours
- `--alias-hh`: also write hour-only aliases `PREFIX:YYYY-MM-DD_HH`

PBS batch example (template provided):
- See `wrf/wrfint_barra_py_v1.pbs` for a ready-to-adapt job script that
  loads a Python environment and calls `barra_to_int.py` over a time range.

Python prerequisites:
- Python 3.8+
- `netCDF4`, `numpy` available in your environment
- The bundled `WPSUtils.py` and `fortran_io.py` must be importable (they
  live alongside the script inside `wrf/`). If needed, run from that folder
  or add it to `PYTHONPATH`.

Attribution:
- The converters and INTERMEDIATE writing logic in this repo are adapted
  from NCAR's era5_to_int project: https://github.com/NCAR/era5_to_int/tree/main
  with local adjustments for BARRA2 variable names, directory layout, and
  additional robustness for file/variable discovery.

Repository layout (high level)

- `wrf/` — NCL scripts, helpers and decks used for WRF/WPS conversions.
- `figure/` — plotted outputs and CSV summaries of extreme events (excluded).
- Notebooks and helper scripts (step1/step3/step4, aws_map, minimum domain tools).

Contact / Authors

For questions about these scripts or running the workflows, open an
issue or contact the repository owner.

-- end --
