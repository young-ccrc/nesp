; probe_time_and_dims.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
  root = "/g/data/w28/chs548/BARRA2_For_WRF"
  cyc  = "20160128T0000Z"
  base = root + "/2016/01/" + cyc + "/nc/"
  fT = addfile(base + "WRFPRS1/air_temp_uv-201601280000-201601281200.nc","r")

  ; ---- time sanity
  time = fT->time
  printVarSummary(time)
  print("time@units = " + time@units)
  tcal = cd_calendar(time, 0)
  print("nt = " + dimsizes(time))
  do i=0,dimsizes(time)-1
    print(i + ": " + sprinti("%0.4i", toint(tcal(i,0))) + "-" + sprinti("%0.2i",toint(tcal(i,1))) \
             + "-" + sprinti("%0.2i",toint(tcal(i,2))) + " " + sprinti("%0.2i",toint(tcal(i,3))) \
             + ":" + sprinti("%0.2i",toint(tcal(i,4))))
  end do

  ; ---- grid sanity
  lat = fT->latitude
  lon = fT->longitude
  printVarSummary(lat)
  printVarSummary(lon)

  ; ---- a single field slice shape
  TT = fT->air_temp_uv
  printVarSummary(TT)
  printMinMax("TT min/max: ", TT(0,0,:,:), True)
end
