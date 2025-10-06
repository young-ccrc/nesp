# Why I chose Python (and not NCL) for BARRA → WRF INTERMEDIATE

## Summary

I moved the BARRA pre-processing workflow to **Python** because it’s simpler to maintain, easier to extend, and matches WPS’s input expectations with fewer moving parts. 
The Python script (`barra_to_int.py`) reads BARRA pressure-level and single-level NetCDF, fixes level ordering, writes WRF **INTERMEDIATE** records directly, and (crucially) populates the **surface proxy level (XLVL=200100)** so metgrid doesn’t see missing values at “lev 0”.

## Why not NCL?

Using NCL for this conversion typically meant extra steps and more complexity, e.g.:

* Additional steps

  * Manually normalizing or infilling below-ground values (e.g., using `frac_time_p_above`)
  * Reversing/ordering pressure levels for WPS (top→bottom)
  * Adding surface proxies for TT/Q/U/V at **200100 Pa** to satisfy `METGRID.TBL`
  * Building the exact projection metadata and grid spacing for WPS

* More complex to keep correct

  * CF-time handling and calendar nuances
  * Matching varied BARRA file/variable names across dated directory trees
  * Writing WPS INTERMEDIATE (Fortran-style records) without a robust, reusable writer

* Sustainability

  * Fewer maintainers and libraries; harder to integrate with modern packaging, linting, testing
  * Slower iteration when we need to add variables, patterns, or diagnostics

## Why Python works better here

* One pass, end-to-end: reads, fixes, and writes INTERMEDIATE in a single run (no handoff files).
* Robust I/O: `netCDF4` handles CF-time; we match multiple filename and variable patterns.
* WPS-aware: uses a proven IntermediateFile writer; sets projection from the files; ensures descending pressure levels; writes the **surface proxy slabs** (T2→TT, Q2→SPECHUMD, U10→UU, V10→VV at 200100).
* Ergonomic CLI:

  * `--layout dated` understands `ROOT/YYYY/MM/YYYYMMDDTHHZ/nc/...`
  * `--outdir` chooses output folder
  * `--alias-hh` also writes hour-only filenames (what metgrid opens)
  * `--emit-2m10m` optionally writes explicit `T2/Q2/U10/V10` in addition to proxies
* Easy to maintain: simple mapping table for new variables; clear logging (`--verbose`).

## How to use the Python script (keep this in mind)

1. Prereqs: Python 3, `numpy`, `netCDF4`, and NCAR’s `WPSUtils.py` / `fortran_io.py` alongside the script.

2. Use infilled pressure-levels: prefer `*.infilled.nc` for TT/UU/VV/GHT/SPECHUMD to avoid below-ground holes.

3. Command to run (typical):

   ```bash
   python barra_to_int.py \
     --root /g/data/w28/chs548/BARRA2_For_WRF \
     --vars TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC \
     --outdir /g/data/w28/yk8692/nesp/wrf \
     --prefix BARRA \
     --alias-hh \
     --verbose \
     2016-01-28_00 2016-01-28_12 6
   ```

   * This writes files like `BARRA:2016-01-28_00:00:00` plus an alias `BARRA:2016-01-28_00` for metgrid.
   * Lev 0 is **filled** because the script writes surface proxies from `WRFSLV` (`temp_scrn`, `qsair_scrn`, `uwnd10m`, `vwnd10m`).
   * Add `--emit-2m10m` if you also want explicit `T2/Q2/U10/V10` fields in `met_em*` (not just proxies).

4. Validation checklist

   * `ncview`/`ncdump` show **no missing values** at lev 0 for TT/UU/VV/SPECHUMD.
   * `PRES` at lev 0 ≈ PSFC; lev 1..N are the constant isobars (100000, 97500, … Pa).
   * `metgrid.log` no longer warns about missing fields or files (`BARRA:YYYY-MM-DD_HH` opens cleanly).

## Bottom line

Python cuts out fragile glue steps, codifies the WPS rules we actually need (level ordering + surface proxies), and gives us a predictable, reproducible path from BARRA NetCDF to metgrid-ready INTERMEDIATE files. 
Keep using `barra_to_int.py` with the flags above; it’s the most reliable way to get fully-populated `met_em*` inputs from BARRA.
