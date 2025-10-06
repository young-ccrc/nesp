; ------------------------------------------------------------------------------
; wrfint_barra_v4.15.ncl
;
; PURPOSE:
;   Convert BARRA2 (BARRA_R2 high-res) NetCDF meteorological fields to WPS/WRF
;   intermediate files suitable for WRF real-data runs. Handles upper-air and
;   surface fields, including pressure, temperature, wind, humidity, and more.
;
; INPUT DATA LAYOUT:
;   /g/data/w28/chs548/BARRA2_For_WRF/YYYY/MM/CYCLE/nc/{WRFPRS1,WRFPRS2,WRFSLV,WRFSURF}/*.nc
;   - Each CYCLE is a forecast cycle (e.g., 20160128T0000Z)
;   - Required variables: air_temp_uv, relhum_uv or spec_hum_uv, geop_ht_uv,
;     wnd_ucmp, wnd_vcmp, frac_time_p_above, mslp, sfc_pres, temp_scrn,
;     qsair_scrn, uwnd10m, vwnd10m
;   - Optional: sfc_temp (for SKINTEMP)
;
; OUTPUT:
;   WPS intermediate files in <OUTDIR>/<FGNAME>:YYYY-MM-DD_HH
;   - Source BARRA2 data at 30‑min resolution are sub-sampled to top-of-hour
;     snapshots (minute == 00). Minutes are discarded to match WPS expectations.
;   - Internal metadata 'date' keeps full YYYY-MM-DD_HH:MM:SS.
;
; USAGE:
;   export YEAR=2016 MONTH=01
;   export OUTDIR="/g/data/w28/yk8692/wps/int"
;   export FGNAME="BARRA2"
;   export BARRA2_ROOT="/g/data/w28/chs548/BARRA2_For_WRF"
;   ncl wrfint_barra_v4.15.ncl
;
;   - YEAR and MONTH select the period to process.
;   - OUTDIR is where WPS files are written.
;   - FGNAME is the prefix for WPS files and should match namelist.wps.
;   - BARRA2_ROOT can be overridden for alternate data locations.
;
; NOTES:
;   - Handles both 1D and 2D lat/lon grids, normalising to 1D for WPS.
;   - Automatically scales/offsets packed NetCDF variables.
;   - Will error if any required variable is missing for a cycle.
;   - Output files are suitable for use with WPS/metgrid; see WRF documentation.
;   - For debugging or extension, see helper functions at top of script.
;
; Author: Youngil Kim
; Contact: youngil.kim@unsw.edu.au
; Last update: 2025-09-15
; ------------------------------------------------------------------------------


; Load common contributed NCL routines
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; ---------- helpers ----------
; ---------------------------------------------------------------------------
; onefile(dir, pattern)
; Purpose: Return the first filename in `dir` matching `pattern` or fatal out.
; Inputs: dir (string) - directory to search
;         pattern (string) - shell glob pattern to match filenames
; Returns: string - first matching filename (absolute or relative)
; Side effects: prints error and exits if no match found
; ---------------------------------------------------------------------------
function onefile(dir:string, pattern:string)
begin
  cmd = "ls " + dir + "/" + pattern + " 2>/dev/null | head -n 1"
  fn = systemfunc(cmd)
  if (ismissing(fn)) then
    fn = ""
  end if
  if (strlen(fn) .eq. 0) then
    print("FATAL: Missing '" + pattern + "' in " + dir)
    exit
  end if
  return fn
end

; ---------------------------------------------------------------------------
; prefer_infilled(dir, stem)
; Purpose: Return an infilled file if available, else a regular .nc file.
; Inputs: dir (string) - directory to search
;         stem (string) - filename stem like "air_temp_uv-*"
; Returns: string - first matching filename
; ---------------------------------------------------------------------------
function prefer_infilled(dir:string, stem:string)
begin
  cmd1 = "ls " + dir + "/" + stem + ".infilled.nc 2>/dev/null | sort | head -n 1"
  fn1 = systemfunc(cmd1)
  if (.not. ismissing(fn1)) then
    if (strlen(fn1) .gt. 0) then
      return fn1
    end if
  end if
  cmd2 = "ls " + dir + "/" + stem + ".nc 2>/dev/null | sort | head -n 1"
  fn2 = systemfunc(cmd2)
  if (ismissing(fn2)) then
    fn2 = ""
  end if
  if (strlen(fn2) .eq. 0) then
    print("FATAL: Missing '" + stem + "(.infilled).nc' in " + dir)
    exit
  end if
  return fn2
end

; ---------------------------------------------------------------------------
; unpack_to_float(x)
; Purpose: Convert a possibly packed NetCDF variable to floating point
; and apply scale_factor/add_offset attributes if present.
; Inputs: x (numeric) - NetCDF variable (may have scale_factor/add_offset)
; Returns: float array with applied scaling and a float _FillValue
; ---------------------------------------------------------------------------
function unpack_to_float(x:numeric)
local y, sf, ao, has_sf, has_ao, has_fv, fv_in, fv_out
begin
  has_sf = isatt(x, "scale_factor")
  has_ao = isatt(x, "add_offset")
  has_fv = isatt(x, "_FillValue")
  fv_out = default_fillvalue("float")

  if (has_sf) then
    sf = tofloat(x@scale_factor)
  else
    sf = 1.0
  end if
  if (has_ao) then
    ao = tofloat(x@add_offset)
  else
    ao = 0.0
  end if

  if (has_fv) then
    fv_in = x@_FillValue
    ; Mask out original fill values BEFORE scaling/offset
    y = where(x.ne.fv_in, tofloat(x)*sf + ao, fv_out)
  else
    y = tofloat(x)*sf + ao
  end if
  y@_FillValue = fv_out
  return y
end

; ---------------------------------------------------------------------------
; to_date_char(tscalar)
; Purpose: Convert a time scalar with proper time units/calendar into a
; human-readable char string "YYYY-MM-DD_HH:MM:SS" using cd_calendar.
; Inputs: tscalar (numeric) - NetCDF time variable (single element)
; Returns: string representing the date/time
; ---------------------------------------------------------------------------
function to_date_char(tscalar:numeric)
local ymdhms, y,mo,d,hh,mm,ss, s
begin
  ; Requires tscalar to have time@units (and optionally @calendar)
  ymdhms = cd_calendar(tscalar, 0)
  ; cd_calendar returns a 2-D array of shape [n x 6], even for a single time (1 x 6).
  ; Flatten to a 1-D 6-element vector to avoid subscript rank mismatches.
  ymdhms1d = ndtooned(ymdhms)
  y  = toint(ymdhms1d(0))
  mo = toint(ymdhms1d(1))
  d  = toint(ymdhms1d(2))
  hh = toint(ymdhms1d(3))
  mm = toint(ymdhms1d(4))
  ss = toint(ymdhms1d(5))

  s  = sprinti("%0.4i",y) + "-" + sprinti("%0.2i",mo) + "-" + sprinti("%0.2i",d) \
     + "_" + sprinti("%0.2i",hh) + ":" + sprinti("%0.2i",mm) + ":" + sprinti("%0.2i",ss)
  return s
end

; ---------------------------------------------------------------------------
; ensure_dir(path)
; Purpose: Create directory `path` (mkdir -p) unless path is empty or '.'
; Inputs: path (string) - directory path to ensure exists
; Side effects: creates directories on filesystem
; ---------------------------------------------------------------------------
procedure ensure_dir(path:string)
begin
  if (path .ne. "" .and. path .ne. ".") then
    system("mkdir -p " + path)
  end if
end

; ---------------------------------------------------------------------------
; die(msg)
; Purpose: Print a fatal message and exit the script.
; Inputs: msg (string) - error message to display
; Side effects: prints to stdout and terminates execution
; ---------------------------------------------------------------------------
procedure die(msg:string)
begin
  print("FATAL: " + msg)
  exit
end

; ---------------------------------------------------------------------------
; build_opt_equil(dchar, lat, lon, level)
; Purpose: Construct an options record used by wrf_wps_write_int for
;          equidistant lat-lon grids (metadata for WPS intermediate files).
; Inputs: dchar (string) - date/time char for filename metadata
;         lat, lon - 1D or 2D latitude/longitude arrays
;         level (float) - pressure/level code (Pa) for this field
; Returns: opt (record) - options record with projection, startloc, etc.
; ---------------------------------------------------------------------------
function build_opt_equil(dchar:string, lat[*], lon[*], level:float)
local opt, dlon, dlat, la, lo, di_la, di_lo
begin
  ; Ensure 1D coordinate vectors
  di_la = dimsizes(lat)
  if (dimsizes(di_la).eq.2) then
    la = tofloat(lat(:,0))
  else
    la = tofloat(lat)
  end if
  di_lo = dimsizes(lon)
  if (dimsizes(di_lo).eq.2) then
    lo = tofloat(lon(0,:))
  else
    lo = tofloat(lon)
  end if

  ; If decreasing, flip so that startloc SWCORNER is valid
  if (dimsizes(la).gt.1 .and. la(1).lt.la(0)) then
    la = la(::-1)
  end if
  if (dimsizes(lo).gt.1 .and. lo(1).lt.lo(0)) then
    lo = lo(::-1)
  end if

  dlon = lo(1) - lo(0)
  dlat = la(1) - la(0)
  if (dlon .lt. 0 .or. dlat .lt. 0) then
    dlon = abs(dlon)
    dlat = abs(dlat)
  end if

  opt = True
  opt@projection        = 0                      ; Equidistant_Lat_Lon
  opt@date              = dchar
  opt@map_source        = "BARRA2_R2_HRES"
  opt@startloc          = "SWCORNER"             ; 8 chars
  opt@startlon          = lo(0)
  opt@startlat          = la(0)
  opt@deltalon          = dlon
  opt@deltalat          = dlat
  opt@is_wind_earth_relative = True
  opt@level             = tofloat(level)
  return opt
end

begin
  ; --- Config ---
  root   = getenv("BARRA2_ROOT")
  if (ismissing(root) .or. strlen(root) .eq. 0) then
    root = "/g/data/w28/chs548/BARRA2_For_WRF"
  end if
  outdir = getenv("OUTDIR")
  if (ismissing(outdir) .or. strlen(outdir) .eq. 0) then
    outdir = "./wrfint_out"
  end if
  fgname = getenv("FGNAME")
  if (ismissing(fgname) .or. strlen(fgname) .eq. 0) then
    fgname = "BARRA2"
  end if
  ensure_dir(outdir)

  year = getenv("YEAR")
  mon  = getenv("MONTH")
  if (ismissing(year) .or. ismissing(mon)) then
    print("FATAL: Provide YEAR and MONTH, e.g., YEAR=2016 MONTH=01")
    exit
  end if

  cyc_dir = root + "/" + year + "/" + mon
  cmd = "ls -1 " + cyc_dir + " | grep -E 'T[0-9]{4}Z' | sort"
  cycles = systemfunc(cmd)
  if (any(ismissing(cycles))) then
    print("FATAL: No cycles found in " + cyc_dir)
    exit
  end if

  ; (QUICK mode removed for clarity)

  print("BARRA2 root : " + root)
  print("Output dir  : " + outdir)
  print("FGNAME      : " + fgname)
  print("Year/Month  : " + year + "/" + mon)
  print("Found " + tostring(dimsizes(cycles)) + " cycles.")

  do ic = 0, dimsizes(cycles)-1
    cyc = cycles(ic)
    base = root + "/" + year + "/" + mon + "/" + cyc + "/nc/"
    dprs1 = base + "WRFPRS1"
    dprs2 = base + "WRFPRS2"
    dslv  = base + "WRFSLV"
    dsurf = base + "WRFSURF"

    print(">>> Cycle: " + cyc)

    f_t  = addfile(prefer_infilled(dprs1, "air_temp_uv-*") , "r")
    f_rh = addfile(prefer_infilled(dprs1, "relhum_uv-*")  , "r")
    f_q  = addfile(prefer_infilled(dprs1, "spec_hum_uv-*") , "r")

    f_gh = addfile(prefer_infilled(dprs2, "geop_ht_uv-*") , "r")
    f_u  = addfile(prefer_infilled(dprs2, "wnd_ucmp-*")   , "r")
    f_v  = addfile(prefer_infilled(dprs2, "wnd_vcmp-*")   , "r")
    f_h  = addfile(prefer_infilled(dprs2, "frac_time_p_above-*") , "r")

    f_mslp = addfile(onefile(dslv, "mslp-*.nc"),      "r")
    f_psfc = addfile(onefile(dslv, "sfc_pres-*.nc"),  "r")
    f_t2   = addfile(onefile(dslv, "temp_scrn-*.nc"), "r")
    f_q2   = addfile(onefile(dslv, "qsair_scrn-*.nc"),"r")
    f_u10  = addfile(onefile(dslv, "uwnd10m-*.nc"),   "r")
    f_v10  = addfile(onefile(dslv, "vwnd10m-*.nc"),   "r")

    ; Required variable presence checks
    if (.not. isfilevar(f_t,  "air_temp_uv")) then
      die("Missing var air_temp_uv in " + f_t@fname)
    end if

    ; Humidity can come from either relhum_uv or spec_hum_uv
    if (.not. isfilevar(f_rh, "relhum_uv") .and. .not. isfilevar(f_q, "spec_hum_uv")) then
      die("Missing humidity: neither relhum_uv in " + f_rh@fname + " nor spec_hum_uv in " + f_q@fname)
    end if

    if (.not. isfilevar(f_gh, "geop_ht_uv")) then
      die("Missing var geop_ht_uv in " + f_gh@fname)
    end if
    if (.not. isfilevar(f_u,  "wnd_ucmp")) then
      die("Missing var wnd_ucmp in " + f_u@fname)
    end if
    if (.not. isfilevar(f_v,  "wnd_vcmp")) then
      die("Missing var wnd_vcmp in " + f_v@fname)
    end if
    if (.not. isfilevar(f_h,  "frac_time_p_above")) then
      die("Missing var frac_time_p_above in " + f_h@fname)
    end if

    if (.not. isfilevar(f_mslp, "mslp")) then
      die("Missing var mslp in " + f_mslp@fname)
    end if
    if (.not. isfilevar(f_psfc, "sfc_pres")) then
      die("Missing var sfc_pres in " + f_psfc@fname)
    end if
    if (.not. isfilevar(f_t2, "temp_scrn")) then
      die("Missing var temp_scrn in " + f_t2@fname)
    end if
    if (.not. isfilevar(f_q2, "qsair_scrn")) then
      die("Missing var qsair_scrn in " + f_q2@fname)
    end if
    if (.not. isfilevar(f_u10, "uwnd10m")) then
      die("Missing var uwnd10m in " + f_u10@fname)
    end if
    if (.not. isfilevar(f_v10, "vwnd10m")) then
      die("Missing var vwnd10m in " + f_v10@fname)
    end if

    if (.not. isfilevar(f_mslp, "mslp"))            then die("Missing var mslp in " + f_mslp@fname)      end if
    if (.not. isfilevar(f_psfc, "sfc_pres"))        then die("Missing var sfc_pres in " + f_psfc@fname)  end if
    if (.not. isfilevar(f_t2, "temp_scrn"))         then die("Missing var temp_scrn in " + f_t2@fname)   end if
    if (.not. isfilevar(f_q2, "qsair_scrn"))        then die("Missing var qsair_scrn in " + f_q2@fname)  end if
    if (.not. isfilevar(f_u10, "uwnd10m"))          then die("Missing var uwnd10m in " + f_u10@fname)    end if
    if (.not. isfilevar(f_v10, "vwnd10m"))          then die("Missing var vwnd10m in " + f_v10@fname)    end if

    ; Optional SKINTEMP (some datasets may not have this)
    tsk_file = systemfunc("ls " + dslv + "/sfc_temp-*.nc 2>/dev/null | head -n 1")
    if (ismissing(tsk_file)) then
      tsk_file = ""
    end if
    have_tsk_file = (strlen(tsk_file) .gt. 0)
    if (have_tsk_file) then
      f_tsk = addfile(tsk_file, "r")
      sfc_temp_exists = isfilevar(f_tsk, "sfc_temp")
    else
      sfc_temp_exists = False
    end if

  ; Coords & time
    time  = f_t->time
    lev   = f_t->pressure           ; hPa or Pa depending on dataset
    lat   = f_t->latitude
    lon   = f_t->longitude
    ; Normalize coordinates to 1D vectors
    di_lat = dimsizes(lat)
    if (dimsizes(di_lat).eq.2) then
      lat = tofloat(lat(:,0))
    else
      lat = tofloat(lat)
    end if
    di_lon = dimsizes(lon)
    if (dimsizes(di_lon).eq.2) then
      lon = tofloat(lon(0,:))
    else
      lon = tofloat(lon)
    end if
    ; Ensure ascending for SWCORNER metadata and track flips
    lat_flipped = False
    lon_flipped = False
    if (dimsizes(lat).gt.1 .and. lat(1).lt.lat(0)) then
      lat = lat(::-1)
      lat_flipped = True
    end if
    if (dimsizes(lon).gt.1 .and. lon(1).lt.lon(0)) then
      lon = lon(::-1)
      lon_flipped = True
    end if

    nt = dimsizes(time)
    nz = dimsizes(lev)
  print("INFO cycle=" + cyc + " dims: time=" + tostring(nt) + " levels=" + tostring(nz))

  ; Read & unpack
    TT   = unpack_to_float(f_t->air_temp_uv)       ; K
    GHT  = unpack_to_float(f_gh->geop_ht_uv)       ; m or geopotential
    ; Normalize unpacked missings before any arithmetic
    miss = -1.e30
    TT   = unpack_to_float(f_t->air_temp_uv)
    TT   = where(ismissing(TT),  miss, TT)         ; K
    GHT  = where(ismissing(GHT), miss, GHT)
    UU   = unpack_to_float(f_u->wnd_ucmp)          ; m/s
    UU   = where(ismissing(UU),  miss, UU)
    VV   = unpack_to_float(f_v->wnd_vcmp)          ; m/s
    VV   = where(ismissing(VV),  miss, VV)
    HFR  = unpack_to_float(f_h->frac_time_p_above) ; fraction (0-1)
    ; Guard HFR: replace missing with 0, clamp to [0,1]
    HFR = where(ismissing(HFR), 0.0, HFR)
    HFR = where(HFR .lt. 0.0, 0.0, HFR)
    HFR = where(HFR .gt. 1.0, 1.0, HFR)

    have_rh = isfilevar(f_rh, "relhum_uv")
    if (have_rh) then
      RH = unpack_to_float(f_rh->relhum_uv)        ; %
      RH = where(ismissing(RH), miss, RH)
    else
      Q  = unpack_to_float(f_q->spec_hum_uv)       ; kg/kg
      Q  = where(ismissing(Q),  miss, Q)
    end if

    PMSL = unpack_to_float(f_mslp->mslp)           ; Pa or hPa (packed)
    if (max(PMSL) .lt. 20000) then
      PMSL = PMSL * 100.0
    end if
    PSFC = unpack_to_float(f_psfc->sfc_pres)       ; Pa (packed)
    T2   = unpack_to_float(f_t2->temp_scrn)        ; K (packed)
    Q2   = unpack_to_float(f_q2->qsair_scrn)       ; kg/kg (packed)
    U10  = unpack_to_float(f_u10->uwnd10m)         ; m/s (packed)
    V10  = unpack_to_float(f_v10->vwnd10m)         ; m/s (packed)
    if (sfc_temp_exists) then
      TSK = unpack_to_float(f_tsk->sfc_temp)       ; K (packed)
    end if

    ; ---------------------------------------------------------------
    ; Interpolate surface (668-lat) fields onto upper-air latitude grid (669)
    ; The upper-air latitudes are half a grid spacing poleward/equatorward.
    ; We perform 1-D linear interpolation along latitude for each (time, lon).
    ; ---------------------------------------------------------------
    lat_prs = f_t->latitude          ; original pressure-level latitude (669)
    if (dimsizes(dimsizes(lat_prs)).eq.2) then
      lat_prs = tofloat(lat_prs(:,0))
    else
      lat_prs = tofloat(lat_prs)
    end if
    ; Ensure ascending orientation for interpolation consistency
    if (lat_prs(1) .lt. lat_prs(0)) then
      lat_prs = lat_prs(::-1)
      ; If pressure lats flipped, arrays already flipped above later, so skip here
    end if

    ; Need original surface latitude (may differ slightly). Re-read from sfc file to avoid
    ; any flips applied to 'lat'. We'll construct local copies.
    lat_sfc_orig = f_psfc->latitude
    if (dimsizes(dimsizes(lat_sfc_orig)).eq.2) then
      lat_sfc_orig = tofloat(lat_sfc_orig(:,0))
    else
      lat_sfc_orig = tofloat(lat_sfc_orig)
    end if
    if (lat_sfc_orig(1) .lt. lat_sfc_orig(0)) then
      lat_sfc_orig = lat_sfc_orig(::-1)
      PMSL = PMSL(::-1,:)
      PSFC = PSFC(::-1,:)
      T2   = T2(::-1,:)
      Q2   = Q2(::-1,:)
      U10  = U10(::-1,:)
      V10  = V10(::-1,:)
      if (sfc_temp_exists) then
        TSK = TSK(::-1,:)
      end if
    end if

    ; ---------------------------------------------------------------
    ; Ensure surface fields have dimension order (time, lat, lon).
    ; Some BARRA-style files may store as (time, lon, lat).
    ; Detect by comparing middle/last dimension sizes with known lat/lon sizes.
    ; ---------------------------------------------------------------
    lon_sfc_orig = f_psfc->longitude
    if (dimsizes(dimsizes(lon_sfc_orig)).eq.2) then
      lon_sfc_orig = tofloat(lon_sfc_orig(0,:))
    else
      lon_sfc_orig = tofloat(lon_sfc_orig)
    end if
    nlat_sfc_tmp = dimsizes(lat_sfc_orig)
    nlon_sfc_tmp = dimsizes(lon_sfc_orig)
    dims_psfc_chk = dimsizes(PSFC)
    if (dimsizes(dims_psfc_chk) .eq. 3) then
      ; Case A: already (time, lat, lon) -> dims = (nt, nlat, nlon) : do nothing
      if (.not. (dims_psfc_chk(1) .eq. nlat_sfc_tmp .and. dims_psfc_chk(2) .eq. nlon_sfc_tmp)) then
        ; Case B: maybe (time, lon, lat) -> dims = (nt, nlon, nlat)
        if (dims_psfc_chk(1) .eq. nlon_sfc_tmp .and. dims_psfc_chk(2) .eq. nlat_sfc_tmp) then
          print("INFO: Reordering surface fields from (time, lon, lat) to (time, lat, lon)")
          ; 'transpose' in NCL permutes the two rightmost dimensions only (or 1D/2D cases).
          ; For arbitrary reordering of a 3D array we must use dim_permute(var, order).
          PSFC = dim_permute(PSFC, (/0,2,1/))
          PMSL = dim_permute(PMSL, (/0,2,1/))
          T2   = dim_permute(T2  , (/0,2,1/))
          Q2   = dim_permute(Q2  , (/0,2,1/))
          U10  = dim_permute(U10 , (/0,2,1/))
          V10  = dim_permute(V10 , (/0,2,1/))
          if (sfc_temp_exists) then
            TSK = dim_permute(TSK, (/0,2,1/))
          end if
        else
          print("WARN: Unexpected surface field dimensions: " + tostring(dims_psfc_chk) + ". Proceeding but interpolation may fail.")
        end if
      end if
    else
      print("WARN: Surface pressure has unexpected rank " + tostring(dimsizes(dims_psfc_chk)) + ".")
    end if

    nlat_prs = dimsizes(lat_prs)
    nlat_sfc = dimsizes(lat_sfc_orig)
    if (nlat_prs .ne. nlat_sfc + 1) then
      print("WARN: Unexpected latitude size relation: prs=" + tostring(nlat_prs) + " sfc=" + tostring(nlat_sfc))
    end if

    ; Allocate interpolated arrays (time, nlat_prs, lon)
    dims_psfc = dimsizes(PSFC)
    nt_sfc = dims_psfc(0)
    nlon = dims_psfc(2)
    PSFC_i = new( (/nt_sfc, nlat_prs, nlon/), typeof(PSFC), miss)
    PMSL_i = new( (/nt_sfc, nlat_prs, nlon/), typeof(PMSL), miss)
    T2_i   = new( (/nt_sfc, nlat_prs, nlon/), typeof(T2),   miss)
    Q2_i   = new( (/nt_sfc, nlat_prs, nlon/), typeof(Q2),   miss)
    U10_i  = new( (/nt_sfc, nlat_prs, nlon/), typeof(U10),  miss)
    V10_i  = new( (/nt_sfc, nlat_prs, nlon/), typeof(V10),  miss)
    if (sfc_temp_exists) then
      TSK_i  = new( (/nt_sfc, nlat_prs, nlon/), typeof(TSK), miss)
    end if

    ; --- Simplified midpoint interpolation (Option A) ---
    if (nlat_prs .eq. nlat_sfc + 1) then
      ; Endpoints copy
      PSFC_i(:,0,:) = PSFC(:,0,:)
      PSFC_i(:,nlat_prs-1,:) = PSFC(:,nlat_sfc-1,:)
      PMSL_i(:,0,:) = PMSL(:,0,:)
      PMSL_i(:,nlat_prs-1,:) = PMSL(:,nlat_sfc-1,:)
      T2_i(:,0,:) = T2(:,0,:); copy
      T2_i(:,nlat_prs-1,:) = T2(:,nlat_sfc-1,:)
      Q2_i(:,0,:) = Q2(:,0,:)
      Q2_i(:,nlat_prs-1,:) = Q2(:,nlat_sfc-1,:)
      U10_i(:,0,:) = U10(:,0,:)
      U10_i(:,nlat_prs-1,:) = U10(:,nlat_sfc-1,:)
      V10_i(:,0,:) = V10(:,0,:)
      V10_i(:,nlat_prs-1,:) = V10(:,nlat_sfc-1,:)
      if (sfc_temp_exists) then
        TSK_i(:,0,:) = TSK(:,0,:)
        TSK_i(:,nlat_prs-1,:) = TSK(:,nlat_sfc-1,:)
      end if
      do j=1, nlat_prs-2
        PSFC_i(:,j,:) = 0.5*(PSFC(:,j-1,:)+PSFC(:,j,:))
        PMSL_i(:,j,:) = 0.5*(PMSL(:,j-1,:)+PMSL(:,j,:))
        T2_i(:,j,:)   = 0.5*(T2(:,j-1,:)+T2(:,j,:))
        Q2_i(:,j,:)   = 0.5*(Q2(:,j-1,:)+Q2(:,j,:))
        U10_i(:,j,:)  = 0.5*(U10(:,j-1,:)+U10(:,j,:))
        V10_i(:,j,:)  = 0.5*(V10(:,j-1,:)+V10(:,j,:))
        if (sfc_temp_exists) then
          TSK_i(:,j,:) = 0.5*(TSK(:,j-1,:)+TSK(:,j,:))
        end if
      end do
    else
      print("WARN: Cannot apply midpoint interpolation cleanly (prs=" + tostring(nlat_prs) + " sfc=" + tostring(nlat_sfc) + "). Copying nearest rows.")
      ; NCL's min() is a reduction that accepts only ONE array argument, not two scalars.
      ; Choose the smaller dimension explicitly.
      if (nlat_prs .lt. nlat_sfc) then
        nn = nlat_prs
      else
        nn = nlat_sfc
      end if
      PSFC_i(:,0:nn-1,:) = PSFC(:,0:nn-1,:)
      PMSL_i(:,0:nn-1,:) = PMSL(:,0:nn-1,:)
      T2_i(:,0:nn-1,:)   = T2(:,0:nn-1,:)
      Q2_i(:,0:nn-1,:)   = Q2(:,0:nn-1,:)
      U10_i(:,0:nn-1,:)  = U10(:,0:nn-1,:)
      V10_i(:,0:nn-1,:)  = V10(:,0:nn-1,:)
      if (sfc_temp_exists) then
        TSK_i(:,0:nn-1,:) = TSK(:,0:nn-1,:)
      end if
      if (nlat_prs .gt. nlat_sfc) then
        ; extend last
        PSFC_i(:,nn:,:) = PSFC(:,nlat_sfc-1,:)
        PMSL_i(:,nn:,:) = PMSL(:,nlat_sfc-1,:)
        T2_i(:,nn:,:)   = T2(:,nlat_sfc-1,:)
        Q2_i(:,nn:,:)   = Q2(:,nlat_sfc-1,:)
        U10_i(:,nn:,:)  = U10(:,nlat_sfc-1,:)
        V10_i(:,nn:,:)  = V10(:,nlat_sfc-1,:)
        if (sfc_temp_exists) then
          TSK_i(:,nn:,:)  = TSK(:,nlat_sfc-1,:)
        end if
      end if
    end if

    ; Replace originals with interpolated versions and set new latitude
    ; Delete originals first because dimension (lat) length changes from nlat_sfc -> nlat_prs.
    delete(PSFC)
    delete(PMSL)
    delete(T2)
    delete(Q2)
    delete(U10)
    delete(V10)
    if (sfc_temp_exists) then
      delete(TSK)
    end if
    PSFC = PSFC_i
    PMSL = PMSL_i
    T2   = T2_i
    Q2   = Q2_i
    U10  = U10_i
    V10  = V10_i
    if (sfc_temp_exists) then
      TSK  = TSK_i
    end if
    ; Use pressure-level latitude as common grid hereafter
    lat = lat_prs

    ; Ensure pressure-level orientation is consistent across streams
    lev_t = lev
    lev_h = f_h->pressure
    asc_t = (lev_t(1) .gt. lev_t(0))
    asc_h = (lev_h(1) .gt. lev_h(0))
    if (asc_h .ne. asc_t) then
      HFR = HFR(:,::-1,:,:)
    end if
    ; WRFPRS2 winds and heights
    lev_gh = f_gh->pressure
    asc_gh = (lev_gh(1) .gt. lev_gh(0))
    if (asc_gh .ne. asc_t) then
      GHT = GHT(:,::-1,:,:)
    end if
    lev_u = f_u->pressure
    asc_u = (lev_u(1) .gt. lev_u(0))
    if (asc_u .ne. asc_t) then
      UU = UU(:,::-1,:,:)
    end if
    lev_v = f_v->pressure
    asc_v = (lev_v(1) .gt. lev_v(0))
    if (asc_v .ne. asc_t) then
      VV = VV(:,::-1,:,:)
    end if
    ; Humidity on pressure levels (if present as RH/Q)
    if (have_rh) then
      lev_rh = f_rh->pressure
      asc_rh = (lev_rh(1) .gt. lev_rh(0))
      if (asc_rh .ne. asc_t) then
        RH = RH(:,::-1,:,:)
      end if
    else
      lev_q = f_q->pressure
      asc_q = (lev_q(1) .gt. lev_q(0))
      if (asc_q .ne. asc_t) then
        Q = Q(:,::-1,:,:)
      end if
    end if

    ; Optional debug: print first/last 5 pressure levels from each stream
    dbg = getenv("DEBUG")
    dbg_on = (.not. ismissing(dbg) .and. (str_lower(dbg).eq."1" .or. str_lower(dbg).eq."true" .or. str_lower(dbg).eq."yes"))

    if (dbg_on) then
      np = dimsizes(lev_t)
      head_n = 5
      if (np .lt. head_n) then
        head_n = np
      end if
      lt_head = tofloat(lev_t(0:head_n-1))
      lt_tail = tofloat(lev_t(np-head_n:np-1))
      lh_head = tofloat(lev_h(0:head_n-1))
      np2 = dimsizes(lev_h)
      lh_tail = tofloat(lev_h(np2-head_n:np2-1))
      lt_head_s = str_join( sprintf("%7.0f", lt_head), " ")
      lt_tail_s = str_join( sprintf("%7.0f", lt_tail), " ")
      lh_head_s = str_join( sprintf("%7.0f", lh_head), " ")
      lh_tail_s = str_join( sprintf("%7.0f", lh_tail), " ")
      print("DEBUG lev_t[0:" + tostring(head_n-1) + "]=" + lt_head_s)
      print("DEBUG lev_t[-" + tostring(head_n) + ":-1]=" + lt_tail_s)
      print("DEBUG lev_h[0:" + tostring(head_n-1) + "]=" + lh_head_s)
      print("DEBUG lev_h[-" + tostring(head_n) + ":-1]=" + lh_tail_s)
    end if

    ; After normalizing missings, apply GHT unit conversion ONLY if clearly geopotential units
    if (isatt(f_gh->geop_ht_uv, "units")) then
      ustr = tostring(f_gh->geop_ht_uv@units)
      ustr_l = str_lower(ustr)
      if (ustr_l .eq. "m2 s-2") then
        GHT = where(.not. ismissing(GHT), GHT/9.80665, GHT)
        GHT@units = "m"
      end if
    end if

  ; Convert pressure levels to Pa if needed, then build 3-D pressure (Pa)
    if (max(lev) .le. 1100.) then
      lev_pa = tofloat(lev)*100.0
    else
      lev_pa = tofloat(lev)
    end if
    P3D = conform(TT, lev_pa, (/1/))
    ; NOTE: Do NOT apply masking yet. We'll mask below-ground points at write-out using PSFC.

  ; Compute RH if needed (Tetens) and clamp to [0,100]
    if (.not. have_rh) then
      tc = TT - 273.15
      es = 610.94 * exp(17.625*tc/(tc + 243.04))       ; Pa
      eps = 0.622
      qs  = eps*es / (P3D - (1.0-eps)*es)
      RH  = 100.0 * (Q/qs)
      RH@units = "%"
    end if
    RH = where(RH.gt.100.0, 100.0, RH)
    RH = where(RH.lt.0.0,     0.0, RH)

  ; Flip arrays if coordinates were reversed to maintain SWCORNER layout
    if (lat_flipped) then
      TT  = TT(:,:,::-1,:)
      GHT = GHT(:,:,::-1,:)
      UU  = UU(:,:,::-1,:)
      VV  = VV(:,:,::-1,:)
      RH  = RH(:,:,::-1,:)
      HFR = HFR(:,:,::-1,:)
      PMSL = PMSL(::-1,:)
      PSFC = PSFC(::-1,:)
      T2   = T2(::-1,:)
      Q2   = Q2(::-1,:)
      U10  = U10(::-1,:)
      V10  = V10(::-1,:)
      if (sfc_temp_exists) then
        TSK = TSK(::-1,:)
      end if
    end if
    if (lon_flipped) then
      TT  = TT(:,:,:,::-1)
      GHT = GHT(:,:,:,::-1)
      UU  = UU(:,:,:,::-1)
      VV  = VV(:,:,:,::-1)
      RH  = RH(:,:,:,::-1)
      HFR = HFR(:,:,:,::-1)
      PMSL = PMSL(:,::-1)
      PSFC = PSFC(:,::-1)
      T2   = T2(:,::-1)
      Q2   = Q2(:,::-1)
      U10  = U10(:,::-1)
      V10  = V10(:,::-1)
      if (sfc_temp_exists) then
        TSK = TSK(:,::-1)
      end if
    end if

    ; Ensure _FillValue is set so min/max ignore missings in debug
    TT@_FillValue  = miss
    GHT@_FillValue = miss
    UU@_FillValue  = miss
    VV@_FillValue  = miss
    RH@_FillValue  = miss
    P3D@_FillValue = miss
    PMSL@_FillValue = miss
    PSFC@_FillValue = miss
    T2@_FillValue   = miss
    Q2@_FillValue   = miss
    U10@_FillValue  = miss
    V10@_FillValue  = miss
    if (sfc_temp_exists) then
      TSK@_FillValue  = miss
    end if

    ; (Time subsetting removed for clarity)

    ; Determine indices corresponding to minute == 0 (top-of-hour)
    all_cal = cd_calendar(time, 0)         ; shape [nt,6]
    ; Minutes column is index 4 (0-based) already integer; just cast
    mins = toint( all_cal(:,4) )
    hour_ids = ind(mins .eq. 0)
    if (.not. any(ismissing(hour_ids))) then
      print("  Using " + tostring(dimsizes(hour_ids)) + " hourly timesteps out of " + tostring(nt))
    else
      if (nt .gt. 0) then
        n_show = nt - 1
        if (n_show .gt. 11) then
          n_show = 11
        end if
        print("DEBUG minute values (first " + tostring(n_show+1) + ")=" \
              + str_join(tostring(mins(0:n_show)),","))
      else
        print("DEBUG no time steps available (nt=0) to show minutes")
      end if
      die("No hourly (minute==0) timesteps found; aborting")
    end if

    ; --- Write per hourly time step with wrf_wps_write_int ---
    do ii = 0, dimsizes(hour_ids)-1
      it = hour_ids(ii)
      dchar_full = to_date_char(time(it))   ; full timestamp
      dchar_file = str_get_cols(dchar_full, 0, 12) ; YYYY-MM-DD_HH
      rootname = outdir + "/" + fgname

      if (dbg_on) then
        print("DEBUG time=" + dchar_full + ": levels=" + tostring(nz))
        do k = 0, nz-1
          g = GHT(it,k,:,:)
          t = TT(it,k,:,:)
          p = P3D(it,k,:,:)
          g@_FillValue = miss
          t@_FillValue = miss
          p@_FillValue = miss
          msg = "lev=" + sprintf("%7.0f", lev_pa(k)) + " Pa | " \
              + "GHT[min=" + sprintf("%10.2f", min(g)) + " max=" + sprintf("%10.2f", max(g)) + "] " \
              + "TT[min=" + sprintf("%8.2f",  min(t)) + " max=" + sprintf("%8.2f",  max(t)) + "] " \
              + "PRES[min=" + sprintf("%10.2f", min(p)) + " max=" + sprintf("%10.2f", max(p)) + "]"
          print(msg)
        end do
      end if

      ; Upper-air: apply PSFC-based below-ground masking (lev_pa > PSFC => below ground)
      ; Iterate with positive stride and map to descending level index so that kdesc=0 is near-surface level
      do kdesc = 0, nz-1
        k = (nz-1) - kdesc
        opt = build_opt_equil(dchar_full, lat, lon, tofloat(lev_pa(k)))
        if (dbg_on) then
          print("DEBUG writing level k=" + tostring(k) + " levPa=" + sprintf("%7.0f", lev_pa(k)))
        end if
        ; Valid above ground where level pressure <= surface pressure
        valid = (lev_pa(k) .le. PSFC(it,:,:))
        pres_out = where(valid, lev_pa(k), miss)
        tt_full  = TT(it,k,:,:)
        uu_full  = UU(it,k,:,:)
        vv_full  = VV(it,k,:,:)
        ght_full = GHT(it,k,:,:)
        rh_full  = RH(it,k,:,:)
        tt_out   = where(valid, tt_full,  miss)
        uu_out   = where(valid, uu_full,  miss)
        vv_out   = where(valid, vv_full,  miss)
        ght_out  = where(valid, ght_full, miss)
        rh_out   = where(valid, rh_full,  miss)
        wrf_wps_write_int(rootname, "PRES", "Pa",  "PRESSURE",    pres_out, opt)
        wrf_wps_write_int(rootname, "TT",   "K",   "TEMPERATURE", tt_out,   opt)
        wrf_wps_write_int(rootname, "UU",   "m s-1", "U COMPONENT", uu_out, opt)
        wrf_wps_write_int(rootname, "VV",   "m s-1", "V COMPONENT", vv_out, opt)
        wrf_wps_write_int(rootname, "RH",   "%",  "REL HUMIDITY",  rh_out,  opt)
        wrf_wps_write_int(rootname, "GHT",  "m",  "GEOP HEIGHT",   ght_out, opt)
        if (dbg_on) then
          di = dimsizes(valid)
          ntot   = di(0)*di(1)
          nvalid = num(valid)
          cov    = 100.0*tofloat(nvalid)/tofloat(ntot)
          masked = ntot - nvalid
          mean_t = avg( where(valid, tt_full, miss) )
          mean_g = avg( where(valid, ght_full, miss) )
          dbg_msg = "DEBUG level=" + sprintf("%7.0f", lev_pa(k)) + "Pa cov=" + sprintf("%6.2f", cov) \
                  + "% masked=" + sprintf("%7i", masked) \
                  + " meanTT=" + sprintf("%8.2f", mean_t) \
                  + " meanGHT=" + sprintf("%9.2f", mean_g)
          print(dbg_msg)
          if ( (lev_pa(k) .gt. 10000.) .and. any( where(valid, ght_full, miss) .gt. 25000.) ) then
            print("WARN: Suspicious large GHT (>25 km) at pressure > 100 hPa")
          end if
        end if
      end do

      ; Surface/single-levels
      opt_sfc = build_opt_equil(dchar_full, lat, lon, 200100.)     ; surface code
      wrf_wps_write_int(rootname, "PSFC", "Pa", "PSFC",     PSFC(it,:,:), opt_sfc)
      wrf_wps_write_int(rootname, "T2",   "K",  "T2",       T2(it,:,:),   opt_sfc)
      wrf_wps_write_int(rootname, "Q2",   "kg kg-1", "Q2",  Q2(it,:,:),   opt_sfc)
      wrf_wps_write_int(rootname, "U10",  "m s-1", "U10",   U10(it,:,:),  opt_sfc)
      wrf_wps_write_int(rootname, "V10",  "m s-1", "V10",   V10(it,:,:),  opt_sfc)
      if (sfc_temp_exists) then
        wrf_wps_write_int(rootname, "SKINTEMP", "K", "TSK", TSK(it,:,:),  opt_sfc)
      end if

      ; Mean sea level pressure
      opt_msl = build_opt_equil(dchar_full, lat, lon, 201300.)     ; MSL code
      wrf_wps_write_int(rootname, "PMSL", "Pa", "MSLP", PMSL(it,:,:), opt_msl)
      if (dbg_on) then
        print("DEBUG finished hour " + dchar_full)
      end if
    end do

    print("<<< Done cycle: " + cyc)
  end do

  print("All cycles complete. Intermediate files under: " + outdir)
  print("Remember to set &metgrid fg_name = '" + fgname + "' in namelist.wps")
end
