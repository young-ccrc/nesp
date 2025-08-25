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
;   WPS intermediate files in <OUTDIR>/<FGNAME>:YYYY-MM-DD_HH:MM:SS
;   - Each file contains a single time and variable, with WRF/WPS metadata.
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
; AUTHOR: Youngil Kim
; Last update: 2025-08-21
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
  if (strlen(fn) .eq. 0) then
    print("FATAL: Missing '" + pattern + "' in " + dir)
    exit
  end if
  return fn
end

; ---------------------------------------------------------------------------
; unpack_to_float(x)
; Purpose: Convert a possibly packed NetCDF variable to floating point
; and apply scale_factor/add_offset attributes if present.
; Inputs: x (numeric) - NetCDF variable (may have scale_factor/add_offset)
; Returns: float array with applied scaling and a float _FillValue
; ---------------------------------------------------------------------------
function unpack_to_float(x:numeric)
local y, sf, ao
begin
  if (isatt(x, "scale_factor")) then
    sf = tofloat(x@scale_factor)
    ao = 0.0
    if (isatt(x, "add_offset")) then
      ao = tofloat(x@add_offset)
    end if
    y = tofloat(x)
    y = y*sf + ao
    y@_FillValue = default_fillvalue("float")
    return y
  else
    return tofloat(x)
  end if
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
  opt@is_wind_earth_rel = False
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

  ; Discover cycles
  cyc_dir = root + "/" + year + "/" + mon
  cmd = "ls -1 " + cyc_dir + " | grep -E 'T[0-9]{4}Z' | sort"
  cycles = systemfunc(cmd)
  if (any(ismissing(cycles))) then
    print("FATAL: No cycles found in " + cyc_dir)
    exit
  end if

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

    f_t  = addfile(onefile(dprs1, "air_temp_uv-*.nc"), "r")
    f_rh = addfile(onefile(dprs1, "relhum_uv-*.nc"),  "r")
    f_q  = addfile(onefile(dprs1, "spec_hum_uv-*.nc"),"r")

    f_gh = addfile(onefile(dprs2, "geop_ht_uv-*.nc"), "r")
    f_u  = addfile(onefile(dprs2, "wnd_ucmp-*.nc"),   "r")
    f_v  = addfile(onefile(dprs2, "wnd_vcmp-*.nc"),   "r")
    f_h  = addfile(onefile(dprs2, "frac_time_p_above-*.nc"), "r")

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
    have_tsk_file = (strlen(tsk_file) .gt. 0)
    if (have_tsk_file) then
      f_tsk = addfile(tsk_file, "r")
      sfc_temp_exists = isfilevar(f_tsk, "sfc_temp")
    else
      sfc_temp_exists = False
    end if

    ; Coords & time
    time  = f_t->time
    lev   = f_t->pressure           ; hPa
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
    ; Ensure ascending for SWCORNER metadata
    if (dimsizes(lat).gt.1 .and. lat(1).lt.lat(0)) then
      lat = lat(::-1)
      ; Flip the data later if needed before writes (upper-air handled per-slice)
    end if
    if (dimsizes(lon).gt.1 .and. lon(1).lt.lon(0)) then
      lon = lon(::-1)
    end if

    nt = dimsizes(time)
    nz = dimsizes(lev)

    ; Read & unpack
    TT   = unpack_to_float(f_t->air_temp_uv)       ; K
    HGTS = unpack_to_float(f_gh->geop_ht_uv)       ; m
    UU   = unpack_to_float(f_u->wnd_ucmp)          ; m/s
    VV   = unpack_to_float(f_v->wnd_vcmp)          ; m/s
    HFR  = tofloat(f_h->frac_time_p_above)

    have_rh = isfilevar(f_rh, "relhum_uv")
    if (have_rh) then
      RH = unpack_to_float(f_rh->relhum_uv)        ; %
    else
      Q  = unpack_to_float(f_q->spec_hum_uv)       ; kg/kg
    end if

    PMSL = tofloat(f_mslp->mslp)                   ; Pa or hPa
    if (max(PMSL) .lt. 20000) then
      PMSL = PMSL * 100.0
    end if
    PSFC = tofloat(f_psfc->sfc_pres)               ; Pa
    T2   = tofloat(f_t2->temp_scrn)                ; K
    Q2   = tofloat(f_q2->qsair_scrn)               ; kg/kg
    U10  = tofloat(f_u10->uwnd10m)                 ; m/s
    V10  = tofloat(f_v10->vwnd10m)                 ; m/s
    if (sfc_temp_exists) then
      TSK = tofloat(f_tsk->sfc_temp)               ; K
    end if

    ; Build 3-D pressure (Pa) and mask sub-surface levels
    P3D = conform(TT, lev*100.0, (/1/))
    TT   = where(HFR .gt. 0.0, TT,   TT@_FillValue)
    HGTS = where(HFR .gt. 0.0, HGTS, HGTS@_FillValue)
    UU   = where(HFR .gt. 0.0, UU,   UU@_FillValue)
    VV   = where(HFR .gt. 0.0, VV,   VV@_FillValue)
    if (have_rh) then
      RH = where(HFR .gt. 0.0, RH, RH@_FillValue)
    else
      Q  = where(HFR .gt. 0.0, Q,  Q@_FillValue)
    end if

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

    ; --- Write per time step with wrf_wps_write_int ---
    do it = 0, nt-1
      dchar = to_date_char(time(it))
      rootname = outdir + "/" + fgname + ":" + dchar

      ; Upper-air: level in Pa
      do k = 0, nz-1
        opt = build_opt_equil(dchar, lat, lon, tofloat(P3D(it,k,0,0)))
        wrf_wps_write_int(rootname, "PRES", "Pa", "PRESSURE",        P3D(it,k,:,:), opt)
        wrf_wps_write_int(rootname, "TT",   "K",  "TEMPERATURE",     TT(it,k,:,:),   opt)
        wrf_wps_write_int(rootname, "UU",   "m s-1", "U COMPONENT",  UU(it,k,:,:),   opt)
        wrf_wps_write_int(rootname, "VV",   "m s-1", "V COMPONENT",  VV(it,k,:,:),   opt)
        wrf_wps_write_int(rootname, "RH",   "%",  "REL HUMIDITY",    RH(it,k,:,:),   opt)
        wrf_wps_write_int(rootname, "HGTS", "m",  "GEOP HEIGHT",     HGTS(it,k,:,:), opt)
      end do

      ; Surface/single-levels
      opt_sfc = build_opt_equil(dchar, lat, lon, 200100.)     ; surface code
      wrf_wps_write_int(rootname, "PSFC", "Pa", "PSFC",     PSFC(it,:,:), opt_sfc)
      wrf_wps_write_int(rootname, "T2",   "K",  "T2",       T2(it,:,:),   opt_sfc)
      wrf_wps_write_int(rootname, "Q2",   "kg kg-1", "Q2",  Q2(it,:,:),   opt_sfc)
      wrf_wps_write_int(rootname, "U10",  "m s-1", "U10",   U10(it,:,:),  opt_sfc)
      wrf_wps_write_int(rootname, "V10",  "m s-1", "V10",   V10(it,:,:),  opt_sfc)
      if (sfc_temp_exists) then
        wrf_wps_write_int(rootname, "SKINTEMP", "K", "TSK", TSK(it,:,:),  opt_sfc)
      end if

      ; Mean sea level pressure
      opt_msl = build_opt_equil(dchar, lat, lon, 201300.)     ; MSL code
      wrf_wps_write_int(rootname, "PMSL", "Pa", "MSLP", PMSL(it,:,:), opt_msl)
    end do

    print("<<< Done cycle: " + cyc)
  end do

  print("All cycles complete. Intermediate files under: " + outdir)
  print("Remember to set &metgrid fg_name = '" + fgname + "' in namelist.wps")
end
