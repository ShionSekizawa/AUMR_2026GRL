;;; File: clim_seasonality.ncl
;;; Description: Calculate monthly rainfall climatology over northern Australia for observations (AGCD) and the d4PDF historical experiment
;;;
;;; Input: AGCD data
;;;        d4PDF data (historical experiment)
;;; Output: clim_seasonality.nc
;;;
;;; Running:
;;; ncl clim_seasonality.ncl
;;;

begin

  ;;; analysis period from y1 to y2
  y1 = 1951
  y2 = 2011
  yyyymm = yyyymm_time(y1,y2,"integer")
  
  ;;; Observation
  in = addfile("${AGCD_dir_path}/v2/land_fixed_0.05.nc","r") ;; AGCD land-sea mask file
  wgt = short2flt(in->land)
  wgt = wgt*conform_dims(dimsizes(wgt),cos(tofloat(in->lat)*atan(1.)/45.),0)
  wgt({:-20},:) = 0.0   ;; wight array for the spatial average
  rain_obs = new(12,float)  ;; monthly rainfall in AGCD
  rain_obs = 0.0
  err_obs = rain_obs        ;; climatological standard deviation
  do i = 0, dimsizes(yyyymm)-1
    ;; read AGCD rainfall data file
    in = addfile("${AGCD_dir_path}/v2/monthly/precip_total_r005_monthly_"+yyyymm(i)+".nc","r")
    a = wgt_areaave2(in->precip,wgt,0)/days_in_month(yyyymm(i)/100,yyyymm(i)%100)  ;; spatial average and monthly mean
    rain_obs(i%12) = rain_obs(i%12) + a
    err_obs(i%12) = err_obs(i%12) + a*a
  end do
  rain_obs = rain_obs/(y2-y1+1)
  err_obs = err_obs/(y2-y1+1)
  err_obs = sqrt(err_obs-rain_obs^2)
  delete(wgt)

  ;;; d4PDF historical
  dir_d4pdf = "${d4PDF_dir_path}/Global/"  ;; directory path to the d4PDF data 
  setfileoption("bin","ReadByteOrder","BigEndian")
  ;; read d4PDF land-sea mask file to prepare a weight array for the spatial average
  wgt = fbindirread("${d4PDF_dir_path}/TopogRatiol_gsmuv_TL319.gd",1,(/320,640/),"float")
  wgt@_FillValue =  -9.99E33
  
  lon = fspan(0,(640-1)*0.5625,640)
  lat = (/-89.570,  -89.013,  -88.453,  -87.892,  -87.331, \
          -86.769,  -86.208,  -85.647,  -85.085,  -84.523, \
          -83.962,  -83.400,  -82.839,  -82.277,  -81.716, \
          -81.154,  -80.592,  -80.031,  -79.469,  -78.908, \
          -78.346,  -77.784,  -77.223,  -76.661,  -76.100, \
          -75.538,  -74.976,  -74.415,  -73.853,  -73.291, \
          -72.730,  -72.168,  -71.607,  -71.045,  -70.483, \
          -69.922,  -69.360,  -68.799,  -68.237,  -67.675, \
          -67.114,  -66.552,  -65.990,  -65.429,  -64.867, \
          -64.306,  -63.744,  -63.182,  -62.621,  -62.059, \
          -61.498,  -60.936,  -60.374,  -59.813,  -59.251, \
          -58.689,  -58.128,  -57.566,  -57.005,  -56.443, \
          -55.881,  -55.320,  -54.758,  -54.196,  -53.635, \
          -53.073,  -52.512,  -51.950,  -51.388,  -50.827, \
          -50.265,  -49.704,  -49.142,  -48.580,  -48.019, \
          -47.457,  -46.895,  -46.334,  -45.772,  -45.211, \
          -44.649,  -44.087,  -43.526,  -42.964,  -42.402, \
          -41.841,  -41.279,  -40.718,  -40.156,  -39.594, \
          -39.033,  -38.471,  -37.909,  -37.348,  -36.786, \
          -36.225,  -35.663,  -35.101,  -34.540,  -33.978, \
          -33.416,  -32.855,  -32.293,  -31.732,  -31.170, \
          -30.608,  -30.047,  -29.485,  -28.924,  -28.362, \
          -27.800,  -27.239,  -26.677,  -26.115,  -25.554, \
          -24.992,  -24.431,  -23.869,  -23.307,  -22.746, \
          -22.184,  -21.622,  -21.061,  -20.499,  -19.938, \
          -19.376,  -18.814,  -18.253,  -17.691,  -17.129, \
          -16.568,  -16.006,  -15.445,  -14.883,  -14.321, \
          -13.760,  -13.198,  -12.636,  -12.075,  -11.513, \
          -10.952,  -10.390,   -9.828,  -9.267,  -8.705, \
          -8.144,  -7.582,  -7.020,  -6.459,  -5.897, \
          -5.335,  -4.774,  -4.212,  -3.651,  -3.089, \
          -2.527,  -1.966,  -1.404,  -0.842,  -0.281, \
          0.281,    0.842,    1.404,    1.966,    2.527, \
          3.089,    3.651,    4.212,    4.774,    5.335, \
          5.897,    6.459,    7.020,    7.582,    8.144, \
          8.705,    9.267,    9.828,   10.390,   10.952, \
          11.513,   12.075,   12.636,   13.198,   13.760, \
          14.321,   14.883,   15.445,   16.006,   16.568, \
          17.129,   17.691,   18.253,   18.814,   19.376, \
          19.938,   20.499,   21.061,   21.622,   22.184, \
          22.746,   23.307,   23.869,   24.431,   24.992, \
          25.554,   26.115,   26.677,   27.239,   27.800, \
          28.362,   28.924,   29.485,   30.047,   30.608, \
          31.170,   31.732,   32.293,   32.855,   33.416, \
          33.978,   34.540,   35.101,   35.663,   36.225, \
          36.786,   37.348,   37.909,   38.471,   39.033, \
          39.594,   40.156,   40.718,   41.279,   41.841, \
          42.402,   42.964,   43.526,   44.087,   44.649, \
          45.211,   45.772,   46.334,   46.895,   47.457, \
          48.019,   48.580,   49.142,   49.704,   50.265, \
          50.827,   51.388,   51.950,   52.512,   53.073, \
          53.635,   54.196,   54.758,   55.320,   55.881, \
          56.443,   57.005,   57.566,   58.128,   58.689, \
          59.251,   59.813,   60.374,   60.936,   61.498, \
          62.059,   62.621,   63.182,   63.744,   64.306, \
          64.867,   65.429,   65.990,   66.552,   67.114, \
          67.675,   68.237,   68.799,   69.360,   69.922, \
          70.483,   71.045,   71.607,   72.168,   72.730, \
          73.291,   73.853,   74.415,   74.976,   75.538, \
          76.100,   76.661,   77.223,   77.784,   78.346, \
          78.908,   79.469,   80.031,   80.592,   81.154, \
          81.716,   82.277,   82.839,   83.400,   83.962, \
          84.523,   85.085,   85.647,   86.208,   86.769, \
          87.331,   87.892,   88.453,   89.013,   89.570/)

  LAT = ind(lat.ge.-20.and.lat.le.-11)
  LON = ind(lon.ge.115.and.lon.le.150)
  wgt = wgt*conform_dims(dimsizes(wgt),latGauWgt(320,"","",""),0)
  
  rain_d4pdf = new(12,float)  ;; monthly rainfal in d4PDF
  rain_d4pdf = 0.0
  err_d4pdf = rain_d4pdf      ;; climatological standard deviation
  do i = 1, 100
    mxxx = "m"+sprinti("%0.3i",i)
    do j = 0, dimsizes(yyyymm)-1
      system("echo read m="+mxxx+", yyyymm="+yyyymm(j))
      in = addfile(dir_d4pdf+"HPB_"+mxxx+"/"+yyyymm(j)+"/sfc_avr_mon_HPB_"+mxxx+"_"+yyyymm(j)+".grib","r")
      a = wgt_areaave2(in->PRECIPI_GDS4_ISBL(LAT,LON),wgt(LAT,LON),0)*60*60*24.   ;; spatial average and monthly mean
      rain_d4pdf(j%12) = rain_d4pdf(j%12) + a
      err_d4pdf(j%12) = err_d4pdf(j%12) + a*a
    end do
  end do
  rain_d4pdf = rain_d4pdf/(100*(y2-y1+1))
  err_d4pdf = err_d4pdf/(100*(y2-y1+1))
  err_d4pdf = sqrt(err_d4pdf-rain_d4pdf^2)

  rain_obs!0 = "month"
  rain_obs&month = ispan(1,12,1)
  copy_VarCoords(rain_obs,err_obs)
  copy_VarCoords(rain_obs,rain_d4pdf)
  copy_VarCoords(rain_obs,err_d4pdf)

  ;; output
  system("rm -f clim_seasonality.nc")
  o = addfile("clim_seasonality.nc","c")
  o->rain_obs = rain_obs
  o->err_obs = err_obs
  o->rain_d4pdf = rain_d4pdf
  o->err_d4pdf = err_d4pdf
  
end