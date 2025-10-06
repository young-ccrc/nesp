#!/usr/bin/env python3
"""
BARRA -> WRF INTERMEDIATE converter (barra_to_int.py)

Robust features:
  - Layouts: ROOT/YYYY/MM/YYYYMMDDTHHZ/nc/<WRFPRS1|WRFPRS2|WRFSLV>/... (default) or flat ROOT/...
  - Pressure-level inputs default to *.infilled.nc; fall back to *.nc if needed
  - Multiple filename patterns per variable (handles small naming differences)
  - Multiple candidate variable names inside each file (handles var-name differences)
  - Writes GHT (height in meters) for geopotential height-like input

Requirements:
  - netCDF4, numpy
  - WPSUtils.py and fortran_io.py (from NCAR helpers) on PYTHONPATH

Example:
  python barra_to_int.py --root /g/data/w28/chs548/BARRA2_For_WRF \
    --vars TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC \
    2016-01-28_00 2016-01-28_12 6
"""
from __future__ import annotations
import argparse, os, sys, glob
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List, Sequence
from netCDF4 import Dataset, num2date
import numpy as np
import WPSUtils
import shutil, contextlib

# --------------------------- Small utilities ---------------------------

def intdate_to_string(utc_int: int) -> str:
    year = utc_int // 1000000
    month = (utc_int // 10000) % 100
    day = (utc_int // 100) % 100
    hour = utc_int % 100
    return f"{year:04d}-{month:02d}-{day:02d}_{hour:02d}:00:00"

def string_to_yyyymmddhh(utc_date: str) -> Tuple[int,int,int,int]:
    # Accepts 'YYYY-MM-DD_HH' or 'YYYYMMDDHH'
    if "_" in utc_date:
        date, hour = utc_date.split("_")
        y,m,d = [int(x) for x in date.split("-")]
        h = int(hour)
        return y,m,d,h
    else:
        utc_date = utc_date.strip()
        y = int(utc_date[0:4]); m=int(utc_date[4:6]); d=int(utc_date[6:8]); h=int(utc_date[8:10])
        return y,m,d,h

def yyyymmddhh_to_int(y,m,d,h) -> int:
    return y*10**6 + m*10**4 + d*10**2 + h

def add_hours(utc_int: int, hours: int) -> int:
    from datetime import datetime, timedelta
    y = utc_int // 10**6
    m = (utc_int // 10**4) % 100
    d = (utc_int // 10**2) % 100
    h = utc_int % 100
    dt = datetime(y,m,d,h) + timedelta(hours=hours)
    return yyyymmddhh_to_int(dt.year, dt.month, dt.day, dt.hour)

def stamp_yyyymmddThhz(y:int,m:int,d:int,h:int) -> str:
    return f"{y:04d}{m:02d}{d:02d}T{h:02d}00Z"

# --------------------------- MapProjection ---------------------------

class MapProjection:
    """ Parameters for WPS intermediate file projection metadata. """
    def __init__(self, projType, startLat, startLon, startI, startJ,
                 deltaLat, deltaLon, dx=0.0, dy=0.0, truelat1=0.0, truelat2=0.0, xlonc=0.0):
        self.projType = projType
        self.startLat = float(startLat)
        self.startLon = float(startLon)
        self.startI = float(startI)
        self.startJ = float(startJ)
        self.deltaLat = float(deltaLat)
        self.deltaLon = float(deltaLon)
        self.dx = float(dx)
        self.dy = float(dy)
        self.truelat1 = float(truelat1)
        self.truelat2 = float(truelat2)
        self.xlonc = float(xlonc)

def build_latlon_proj_from_file(nc: Dataset) -> MapProjection:
    lat = nc.variables['lat'][:] if 'lat' in nc.variables else nc.variables['latitude'][:]
    lon = nc.variables['lon'][:] if 'lon' in nc.variables else nc.variables['longitude'][:]
    dlat = float(lat[1] - lat[0])
    dlon = float(lon[1] - lon[0])
    startLat = float(lat[-1])   # northernmost row first in WPS LatLon
    startLon = float(lon[0])    # westernmost
    return MapProjection(WPSUtils.Projections.LATLON, startLat, startLon, 1.0, 1.0, -abs(dlat), abs(dlon))

# --------------------------- Writing helper ---------------------------

def write_slab(intfile: WPSUtils.IntermediateFile, slab2d: np.ndarray, xlvl: float,
               proj: MapProjection, WPSname: str, hdate: str, units: str,
               map_source: str, desc: str):
    missing_value = -1.0e30
    if np.ma.isMaskedArray(slab2d):
        arr = slab2d.filled(missing_value)
    else:
        arr = np.where(np.isfinite(slab2d), slab2d, missing_value)
    intfile.write_next_met_field(
        5, arr.shape[1], arr.shape[0], proj.projType, 0.0, xlvl,
        proj.startLat, proj.startLon, proj.startI, proj.startJ,
        proj.deltaLat, proj.deltaLon, proj.dx, proj.dy, proj.xlonc,
        proj.truelat1, proj.truelat2, 6371229.0, 0, WPSname, hdate, units,
        map_source, desc, arr
    )

# --------------------------- Variable map ---------------------------

@dataclass
class VarInfo:
    WPS: str
    patterns: Sequence[str]      # list of file patterns under the base dir
    nc_candidates: Sequence[str] # candidate variable names inside the file
    kind: str                    # 'pl' or 'sl'
    units: str
    desc: str

# Default map includes multiple patterns and var-name candidates
DEFAULT_MAP = {
    # -------- Pressure-level (use *.infilled.nc, fall back to .nc) --------
    'TT': VarInfo('TT',
        patterns=(
            'WRFPRS1/air_temp_uv-barra_r2-hres-*.infilled.nc',
            'WRFPRS1/air_temp_uv-barra_r2-hres-*.nc',
        ),
        nc_candidates=('air_temp_uv','air_temp','t','ta'),
        kind='pl', units='K', desc='Temperature'),

    'SPECHUMD': VarInfo('SPECHUMD',
        patterns=(
            'WRFPRS1/spec_hum_uv-barra_r2-hres-*.infilled.nc',
            'WRFPRS1/spec_hum_uv-barra_r2-hres-*.nc',
        ),
        nc_candidates=('spec_hum_uv','spec_hum','q','qv','specific_humidity'),
        kind='pl', units='kg kg-1', desc='Specific humidity'),

    'UU': VarInfo('UU',
        patterns=(
            'WRFPRS2/wnd_ucmp-barra_r2-hres-*.infilled.nc',
            'WRFPRS2/wnd_ucmp-barra_r2-hres-*.nc',
        ),
        nc_candidates=('wnd_ucmp','u','uwnd','u_wind'),
        kind='pl', units='m s-1', desc='Zonal wind'),

    'VV': VarInfo('VV',
        patterns=(
            'WRFPRS2/wnd_vcmp-barra_r2-hres-*.infilled.nc',
            'WRFPRS2/wnd_vcmp-barra_r2-hres-*.nc',
        ),
        nc_candidates=('wnd_vcmp','v','vwnd','v_wind'),
        kind='pl', units='m s-1', desc='Meridional wind'),

    'GHT': VarInfo('GHT',
        patterns=(
            'WRFPRS2/geop_ht_uv-barra_r2-hres-*.infilled.nc',
            'WRFPRS2/geop_ht_uv-barra_r2-hres-*.nc',
        ),
        nc_candidates=('geop_ht_uv','geop_ht','z','gh'),
        kind='pl', units='m', desc='Geopotential height (m)'),

    # -------- Single-level (plain .nc in WRFSLV) --------
    'PMSL': VarInfo('PMSL',
        patterns=('WRFSLV/mslp-barra_r2-hres-*.nc',),
        nc_candidates=('mslp','mean_sea_level_prs','mean_sea_level_pressure','prmsl'),
        kind='sl', units='Pa', desc='Mean sea level pressure'),

    'PSFC': VarInfo('PSFC',
        patterns=('WRFSLV/sfc_pres-barra_r2-hres-*.nc',),
        nc_candidates=('sfc_pres','surface_pressure','psfc'),
        kind='sl', units='Pa', desc='Surface pressure'),

    # Optional extras if you want them later
    'TT2': VarInfo('TT',
        patterns=('WRFSLV/temp_scrn-barra_r2-hres-*.nc',),
        nc_candidates=('temp_scrn','air_temp_2m','t2m','t2'),
        kind='sl', units='K', desc='2-m temperature'),

    'Q2': VarInfo('SPECHUMD',
        patterns=('WRFSLV/qsair_scrn-barra_r2-hres-*.nc',),
        nc_candidates=('qsair_scrn','q2','q_2m'),
        kind='sl', units='kg kg-1', desc='2-m specific humidity'),

    'RH2': VarInfo('RH',
        patterns=('WRFSLV/rh_scrn-barra_r2-hres-*.nc',),
        nc_candidates=('rh_scrn',),
        kind='sl', units='%', desc='2-m relative humidity'),

    'U10': VarInfo('UU',
        patterns=('WRFSLV/uwnd10m-barra_r2-hres-*.nc',),
        nc_candidates=('uwnd10m','u10','u_10m','wnd_ucmp_10m'),
        kind='sl', units='m s-1', desc='10-m U wind'),

    'V10': VarInfo('VV',
        patterns=('WRFSLV/vwnd10m-barra_r2-hres-*.nc',),
        nc_candidates=('vwnd10m','v10','v_10m','wnd_vcmp_10m'),
        kind='sl', units='m s-1', desc='10-m V wind'),
}

# --- Single-level helper patterns for writing surface proxy slabs (XLVL=200100)
SL_HELPERS = {
    # from your WRFSLV listing
    'T2': {
        'patterns': ('WRFSLV/temp_scrn-barra_r2-hres-*.nc',
                     'WRFSLV/2m_air_temp-*.nc'),
        'candidates': ('temp_scrn','air_temp_2m','t2m','t2'),
        'wps': 'TT',
        'units': 'K',
        'desc': '2-m temp as surface TT'
    },
    'Q2': {
        'patterns': ('WRFSLV/qsair_scrn-barra_r2-hres-*.nc',),
        'candidates': ('qsair_scrn','q2','q_2m'),
        'wps': 'SPECHUMD',
        'units': 'kg kg-1',
        'desc': '2-m specific humidity as surface Q'
    },
    'U10': {
        'patterns': ('WRFSLV/uwnd10m-barra_r2-hres-*.nc',),
        'candidates': ('uwnd10m','u10','u_10m','wnd_ucmp_10m'),
        'wps': 'UU',
        'units': 'm s-1',
        'desc': '10-m U as surface wind'
    },
    'V10': {
        'patterns': ('WRFSLV/vwnd10m-barra_r2-hres-*.nc',),
        'candidates': ('vwnd10m','v10','v_10m','wnd_vcmp_10m'),
        'wps': 'VV',
        'units': 'm s-1',
        'desc': '10-m V as surface wind'
    },
    'SKINTEMP': {
        'patterns': ('WRFSLV/sfc_temp-barra_r2-hres-*.nc',),
        'candidates': ('sfc_temp','skin_temperature','skt'),
        'wps': 'SKINTEMP',
        'units': 'K',
        'desc': 'skin temperature'
    }
}

# --------------------------- Path resolution ---------------------------

def base_dir_for_time(root: str, layout: str, y: int, m: int, d: int, h: int) -> str:
    if layout == 'dated':
        stamp = stamp_yyyymmddThhz(y,m,d,h)
        return os.path.join(root, f"{y:04d}", f"{m:02d}", stamp, "nc")
    else:
        return root

def find_file_multi(root: str, patterns: Sequence[str], layout: str, y: int, m: int, d: int, h: int) -> Optional[str]:
    base = base_dir_for_time(root, layout, y, m, d, h)
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(base, pat)))
        if hits:
            return hits[0]
    # Fallbacks for dated layout (non-padded month, etc.)
    if layout == 'dated':
        alt_base = os.path.join(root, f"{y:04d}", str(m), stamp_yyyymmddThhz(y,m,d,h), "nc")
        for pat in patterns:
            hits = sorted(glob.glob(os.path.join(alt_base, pat)))
            if hits:
                return hits[0]
    return None

# --------------------------- netCDF helpers ---------------------------

def find_time_index(nc: Dataset, y: int, m: int, d: int, h: int) -> int:
    if 'time' in nc.variables:
        t = nc.variables['time']
        cal = getattr(t, 'calendar', 'standard')
        times = num2date(t[:], t.units, cal)
        for i, dt in enumerate(times):
            if (dt.year, dt.month, dt.day, dt.hour) == (y, m, d, h):
                return i
        return -1
    if 'utc_date' in nc.variables:
        utc = [''.join(ch.decode() if isinstance(ch, bytes) else ch for ch in row).strip() for row in nc.variables['utc_date'][:]]
        stamp = f"{y:04d}-{m:02d}-{d:02d}_{h:02d}"
        try:
            return utc.index(stamp)
        except ValueError:
            return -1
    return -1

def read_var_2d_by_candidates(nc: Dataset, candidates: Sequence[str], tidx: int) -> Tuple[np.ndarray, str]:
    for name in candidates:
        if name in nc.variables and nc.variables[name].ndim >= 2:
            arr = nc.variables[name][tidx, ...]
            return np.array(arr, dtype=np.float32), name
    raise KeyError("None of candidate names found: " + ",".join(candidates))

def read_var_3d_by_candidates(nc: Dataset, candidates: Sequence[str], tidx: int) -> Tuple[np.ndarray, np.ndarray, str, str]:
    for name in candidates:
        if name in nc.variables:
            var = nc.variables[name][tidx, ...]
            # find a level dim
            for levname in ('level','isobaricInhPa','plev','pressure'):
                if levname in nc.dimensions or levname in nc.variables:
                    lev = nc.variables[levname][:]
                    if var.ndim == 3 and len(lev) == var.shape[0]:
                        return np.array(var, dtype=np.float32), np.array(lev, dtype=np.float32), levname, name
    # fallback try first candidate anyway
    raise KeyError("None of candidate 3D vars found: " + ",".join(candidates))

def ensure_descending_levels(arr3d: np.ndarray, levvals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if levvals[0] < levvals[-1]:
        return arr3d[::-1, :, :], levvals[::-1]
    return arr3d, levvals

# --------------------------- Main conversion ---------------------------

def convert_time(root: str, layout: str, var_select: List[str], utc_date: str, prefix: str,
                 verbose: bool=False, outdir: str=".", alias_hh: bool=False):

    import shutil  # for alias copy fallback

    y, m, d, h = string_to_yyyymmddhh(utc_date)
    hdate_full = intdate_to_string(yyyymmddhh_to_int(y, m, d, h))   # 'YYYY-MM-DD_HH:00:00'

    # write in chosen directory
    cwd = os.getcwd()
    try:
        os.chdir(outdir)
        out = WPSUtils.IntermediateFile(prefix, hdate_full)

        basedir = base_dir_for_time(root, layout, y, m, d, h)
        last_proj = None  # reuse if a proxy file lacks coords

        # ---------- MAIN: user-requested variables ----------
        for key in var_select:
            if key not in DEFAULT_MAP:
                print(f"[WARN] Unknown variable key '{key}' – skipping."); continue
            info = DEFAULT_MAP[key]
            path = find_file_multi(root, info.patterns, layout, y, m, d, h)
            if not path:
                print(f"[WARN] No file for {key} under {basedir} matching {info.patterns} – skipping."); continue

            with Dataset(path) as nc:
                tidx = find_time_index(nc, y, m, d, h)
                if tidx < 0:
                    print(f"[WARN] Time {utc_date} not found in {os.path.basename(path)} – skipping."); continue

                proj = build_latlon_proj_from_file(nc)
                last_proj = proj
                units = info.units
                desc = info.desc
                map_source = "BARRA reanalysis grid"

                if info.kind == 'sl':
                    slab, used = read_var_2d_by_candidates(nc, info.nc_candidates, tidx)
                    xlvl = 200100.0  # conventional surface code for SL variables
                    write_slab(out, slab, xlvl, proj, info.WPS, hdate_full, units, map_source, desc)
                    if verbose: print(f"[OK] {key}: file={os.path.basename(path)} var={used}")

                elif info.kind == 'pl':
                    arr3d, levvals, levname, used = read_var_3d_by_candidates(nc, info.nc_candidates, tidx)
                    arr3d, levvals = ensure_descending_levels(arr3d, levvals)
                    for k in range(arr3d.shape[0]):
                        slab = arr3d[k, :, :]
                        xlvl = float(levvals[k]) * 100.0  # hPa -> Pa
                        write_slab(out, slab, xlvl, proj, info.WPS, hdate_full, units, map_source, desc)
                    if verbose: print(f"[OK] {key}: file={os.path.basename(path)} var={used} lev={levname} nlev={arr3d.shape[0]}")

        # ---------- NEW: surface proxy slabs at XLVL=200100 ----------
        # Skip proxies that the user already wrote via --vars:
        # (TT2 writes TT@200100, U10->UU@200100, V10->VV@200100, Q2->SPECHUMD@200100, SKINTEMP unchanged)
        skip_if_present = {'T2': 'TT2', 'Q2': 'Q2', 'U10': 'U10', 'V10': 'V10', 'SKINTEMP': 'SKINTEMP'}
        proxies = ['T2','Q2','U10','V10','SKINTEMP']

        for p in proxies:
            if skip_if_present[p] in var_select:
                if verbose: print(f"[INFO] Skipping proxy {p} because {skip_if_present[p]} already in --vars")
                continue

            cfg = SL_HELPERS[p]
            path = find_file_multi(root, cfg['patterns'], layout, y, m, d, h)
            if not path:
                if verbose: print(f"[WARN] {p}: no file in {basedir} matching {cfg['patterns']}")
                continue

            with Dataset(path) as nc:
                tidx = find_time_index(nc, y, m, d, h)
                if tidx < 0:
                    if verbose: print(f"[WARN] {p}: time {utc_date} not found in {os.path.basename(path)}")
                    continue

                # derive proj; fall back to last good one
                try:
                    proj = build_latlon_proj_from_file(nc)
                except Exception:
                    proj = last_proj
                if proj is None:
                    if verbose: print(f"[WARN] {p}: no projection available, skipping")
                    continue

                try:
                    slab, used = read_var_2d_by_candidates(nc, cfg['candidates'], tidx)
                except KeyError:
                    if verbose: print(f"[WARN] {p}: none of {cfg['candidates']} in {os.path.basename(path)}")
                    continue

                write_slab(out, slab, 200100.0, proj, cfg['wps'], hdate_full, cfg['units'], "BARRA reanalysis grid", cfg['desc'])
                if verbose: print(f"[OK] {p}-> {cfg['wps']}@200100: file={os.path.basename(path)} var={used}")

        out.close()

        # Optional hour-only alias (so metgrid opens PREFIX:YYYY-MM-DD_HH)
        if alias_hh:
            src = f"{prefix}:{hdate_full}"
            dst = f"{prefix}:{hdate_full.split(':')[0]}"  # strip ':MM:SS'
            if not os.path.exists(dst):
                try:
                    os.symlink(src, dst)
                except Exception:
                    shutil.copy2(src, dst)

    finally:
        os.chdir(cwd)

        
def main():
    ap = argparse.ArgumentParser(description="Convert BARRA NetCDF to WRF INTERMEDIATE files (robust patterns)")
    ap.add_argument("--root", required=True, help="Root directory of BARRA files")
    ap.add_argument("--outdir", default=".", help="Directory to write WPS intermediate files")
    ap.add_argument("--alias-hh", action="store_true",
                    help="Also create an hour-only alias (PREFIX:YYYY-MM-DD_HH) for metgrid")
    ap.add_argument("--layout", choices=['dated','flat'], default='dated',
                    help="Directory layout under ROOT: dated=ROOT/YYYY/MM/YYYYMMDDTHHZ/nc/... (default), flat=ROOT/...")
    ap.add_argument("--vars", default="TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC",
                    help="Comma-separated keys: TT,UU,VV,SPECHUMD,GHT,PMSL,PSFC,TT2,U10,V10 (and others you add)")
    ap.add_argument("--prefix", default="FILE", help="Prefix for WPS INTERMEDIATE output (default: FILE)")
    ap.add_argument("--verbose", action="store_true", help="Print which files/variables were used")
    ap.add_argument("start", help="Start time YYYY-MM-DD_HH or YYYYMMDDHH")
    ap.add_argument("end", help="End time YYYY-MM-DD_HH or YYYYMMDDHH (inclusive)")
    ap.add_argument("interval", type=int, help="Interval hours (e.g., 6)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    var_select = [v.strip() for v in args.vars.split(",") if v.strip()]
    start_y, start_m, start_d, start_h = string_to_yyyymmddhh(args.start)
    end_y, end_m, end_d, end_h = string_to_yyyymmddhh(args.end)
    curr = yyyymmddhh_to_int(start_y, start_m, start_d, start_h)
    endi = yyyymmddhh_to_int(end_y, end_m, end_d, end_h)

    if not os.path.isdir(args.root):
        print(f"[ERROR] --root not found: {args.root}"); sys.exit(1)

    while True:
        datestr = intdate_to_string(curr)
        print(f"[INFO] Writing INTERMEDIATE for {datestr}")
        # convert_time(args.root, args.layout, var_select, datestr[:13], args.prefix, verbose=args.verbose)
        convert_time(args.root, args.layout, var_select, datestr[:13], args.prefix,
             verbose=args.verbose, outdir=args.outdir, alias_hh=args.alias_hh)
        if curr == endi: break
        curr = add_hours(curr, args.interval)

if __name__ == "__main__":
    main()
