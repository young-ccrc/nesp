;==============================================================
; prep_barra2_common.ncl  -  NCL 6.6.2 helpers for BARRA2->WRFINT
;==============================================================

undef("b2_scale_to_float")
function b2_scale_to_float(x:numeric)
begin
  y = tofloat(x)
  if (isatt(x, "scale_factor")) then
    y = y * x@scale_factor
  end if
  if (isatt(x, "add_offset")) then
    y = y + x@add_offset
  end if
  if (isatt(x, "_FillValue")) then
    y@_FillValue = tofloat(x@_FillValue)
  else
    y@_FillValue = default_fillvalue("float")
  end if
  return(y)
end

undef("b2_time_select_and_strings")
function b2_time_select_and_strings(time[*]:numeric, units:string, every_minutes:integer)
begin
  ; cd_calendar(time, units, option) - option=0 -> [year,mon,day,hour,minute,second]
  tcal = cd_calendar(time, units, 0)

  nt   = dimsizes(time)
  step = max((/1, toint(every_minutes/30)/))  ; 30-min native BARRA2 -> 6-hr = step 12
  idx  = ispan(0, nt-1, step)

  yy = toint(tcal(:,0)); mm = toint(tcal(:,1)); dd = toint(tcal(:,2))
  hh = toint(tcal(:,3)); mi = toint(tcal(:,4))

  nsel = dimsizes(idx)
  timestr = new(nsel, string)
  do i=0, nsel-1
    ii = idx(i)
    timestr(i) = sprinti("%0.4i",yy(ii)) + "-" + \
                 sprinti("%0.2i",mm(ii)) + "-" + \
                 sprinti("%0.2i",dd(ii)) + "_" + \
                 sprinti("%0.2i",hh(ii)) + ":" + \
                 sprinti("%0.2i",mi(ii)) + ":00"
  end do

  out = True
  out@idx     = idx
  out@timestr = timestr
  return(out)
end

undef("b2_latlon_2d")
function b2_latlon_2d(lat1d[*]:numeric, lon1d[*]:numeric)
begin
  LAT2D = lat1d(:,new) + 0.0*lon1d(new,:)
  LON2D = lon1d(new,:) + 0.0*lat1d(:,new)
  out = True
  out@LAT2D = LAT2D
  out@LON2D = LON2D
  return(out)
end

undef("b2_mask_below_ground")
function b2_mask_below_ground(var:numeric, frac:numeric)
begin
  vf = tofloat(var)
  if (isatt(var,"_FillValue")) then
    miss = tofloat(var@_FillValue)
  else
    miss = default_fillvalue("float")
  end if

  ff = tofloat(frac)
  if (isatt(frac, "scale_factor")) then
    ff = ff * frac@scale_factor
  end if
  ; add_offset is 0 for this var in BARRA2

  where (ff.eq.0.0)
    vf = miss
  end where
  vf@_FillValue = miss
  return(vf)
end

undef("b2_regrid_slv_to_uv")
function b2_regrid_slv_to_uv(slv_3d[*][*][*]:numeric, lat_slv[*]:numeric, lon_slv[*]:numeric, \
                             lat_uv[*]:numeric, lon_uv[*]:numeric)
begin
  nt  = dimsizes(slv_3d(:,0,0))
  nyu = dimsizes(lat_uv)
  nxu = dimsizes(lon_uv)
  out = new((/nt, nyu, nxu/), float, getFillValue(tofloat(slv_3d)))
  do t=0, nt-1
    ; linint2_Wrap(fi(lat,lon), fi_lon, fi_lat, fo_lon, fo_lat, cyclic)
    out(t,:,:) = linint2_Wrap( tofloat(slv_3d(t,:,:)), lon_slv, lat_slv, lon_uv, lat_uv, False )
  end do
  copy_VarMeta(slv_3d, out)
  return(out)
end
