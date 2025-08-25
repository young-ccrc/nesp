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

This repository contains an NCL script that converts BARRA2 NetCDF
fields into WPS/WRF intermediate files suitable for use with `metgrid`
and WRF real-data runs.

Script of interest:
- `wrf/wrfint_barra_v4.15.ncl`

Usage (example)

Set the environment variables expected by the NCL script and run it with NCL:

```bash
export YEAR=2016
export MONTH=01
export OUTDIR="/path/to/wps/int"
export FGNAME="BARRA2"
export BARRA2_ROOT="/g/data/w28/chs548/BARRA2_For_WRF"

# Run the NCL conversion (requires NCL installed and NCARG root configured)
ncl wrf/wrfint_barra_v4.15.ncl
```

Key notes about the NCL script (`wrf/wrfint_barra_v4.15.ncl`):
- Expects BARRA2 files arranged as `/.../BARRA2_For_WRF/YYYY/MM/<CYCLE>/nc/`
  with subfolders `WRFPRS1`, `WRFPRS2`, `WRFSLV`, `WRFSURF` containing
  variables such as `air_temp_uv`, `geop_ht_uv`, `wnd_ucmp`, `wnd_vcmp`,
  `mslp`, `sfc_pres`, `temp_scrn`, `qsair_scrn`, `uwnd10m`, `vwnd10m`.
- Environment variables control input month/year, output directory and
  the `fgname` used in `namelist.wps`.
- The script normalizes lat/lon grids, unpacks scaled NetCDF variables,
  computes relative humidity if only specific humidity is present, and
  writes one WPS intermediate file per variable/time slice.

Prerequisites for running the NCL script
- NCL (NCAR Command Language) installed and `NCARG_ROOT` set.
- NetCDF libraries accessible to NCL.
- Access to the BARRA2 data tree described above.

Repository layout (high level)

- `wrf/` — NCL scripts, helpers and decks used for WRF/WPS conversions.
- `figure/` — plotted outputs and CSV summaries of extreme events (excluded).
- Notebooks and helper scripts (step1/step3/step4, aws_map, minimum domain tools).

Git & large files

- This repository intentionally ignores large data outputs. See the
  top-level `.gitignore` which contains entries such as `*.nc`, station
  and figure directories.
- If you need to store large binary/dataset artifacts in the remote
  repository, use Git LFS. For one-off uploads, keep large datasets
  outside the git repository and provide download links in the README.

Troubleshooting

- Push very slow or stalled? Check `git status` and `git ls-files` to see
  which large files are tracked. Use `git rm --cached <path>` to stop
  tracking an accidentally added large file, commit, then push.
- If large files are already present in history and you need to remove
  them permanently, use `git filter-repo` or BFG and then force-push.
  This rewrites history and needs coordination with collaborators.

Contact / Authors

For questions about these scripts or running the workflows, open an
issue or contact the repository owner.

-- end --
