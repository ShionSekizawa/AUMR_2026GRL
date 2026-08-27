;;; File: make_lhf_area.ncl
;;; Description: Calculate JF-mean latent heat flux and related variables averaged over the particular domain
;;;
;;; Input: d4PDF data
;;; Output: lhf_area_${experiment_type_name}.nc
;;;
;;; Running:
;;; ncl make_lhf_area.ncl I=${experiment_type_index}
;;;
;;; ${experiment_type_index} is given as an integer "I"
;;;   I= 0  HPB
;;;   I= 1  HPB_NAT
;;;   I= 2  HFB_2K_CC
;;;   I= 3  HFB_2K_GF
;;;   I= 4  HFB_2K_HA
;;;   I= 5  HFB_2K_MI
;;;   I= 6  HFB_2K_MP
;;;   I= 7  HFB_2K_MR
;;;   I= 8  HFB_4K_CC
;;;   I= 9  HFB_4K_GF
;;;   I=10  HFB_4K_HA
;;;   I=11  HFB_4K_MI
;;;   I=12  HFB_4K_MP
;;;   I=13  HFB_4K_MR
;;;

;;; function for calculating saturation specific humidity with respect to SST
function sst2qs(sst:numeric,ps:numeric)
local t0, a1, a3, a4
begin
  t0 = 273.16 ; K
  a1 = 611.21 ; Pa
  a3 = 17.502
  a4 = 32.19
  esat = a1*exp(a3*((sst-t0)/(sst-a4)))
  eps = 0.621981
  qs = eps*esat/(ps-(1.-eps)*esat)
  return(qs)
end

;;; function for reading d4PDF data, performing spatial average, and calculating the JF-mean
function jf_mean(ifile[1]:string,yyyy[1]:integer,var[1]:string,wgt[*][*]:numeric,LAT[*]:integer,LON[*]:integer)
local in1, in2, a1, a2, d1, d2
begin
  in1 = addfile(str_sub_str(ifile,"##","01"),"r")
  in2 = addfile(str_sub_str(ifile,"##","02"),"r")
  if (var.eq."qs") then
    a1 = wgt_areaave2(sst2qs(in1->TGEF_GDS4_ISBL(LAT,LON),in1->PS_GDS4_ISBL(LAT,LON)),wgt(LAT,LON),0)
    a2 = wgt_areaave2(sst2qs(in2->TGEF_GDS4_ISBL(LAT,LON),in2->PS_GDS4_ISBL(LAT,LON)),wgt(LAT,LON),0)
  else
    a1 = wgt_areaave2(in1->$var$(LAT,LON),wgt(LAT,LON),0)
    a2 = wgt_areaave2(in2->$var$(LAT,LON),wgt(LAT,LON),0)
  end if
  d1 = days_in_month(yyyy,1)
  d2 = days_in_month(yyyy,2)
  return((a1*d1+a2*d2)/(d1+d2))
end


begin
  ;;; directory path referring to the d4PDF
  dir_d4pdf = "${d4PDF_dir_path}/Global/"

  ;;; select experiment type
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(I,0)
  m1 = toint(data_list(I,1))
  m2 = toint(data_list(I,2))
  y1 = toint(data_list(I,3))
  y2 = toint(data_list(I,4))

  ;;; set lat and lon ranges for the spatial average
  s1 = -20
  s2 =  -6
  s3 = 105
  s4 = 150

  ;; prepare the weight array (sea) for the spatial average
  setfileoption("bin","ReadByteOrder","BigEndian")
  ;; read d4PDF land-sea mask file
  land = fbindirread("${d4PDF_dir_path}/TopogRatiol_gsmuv_TL319.gd",1,(/320,640/),"float")
  lon  = fspan(0,(640-1)*0.5625,640)
  lat  = (/-89.570,  -89.013,  -88.453,  -87.892,  -87.331, \
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

  land!0 = "lat"
  land!1 = "lon"
  land&lat = lat
  land&lon = lon
  sea = 1. - land
  sea = sea*conform_dims(dimsizes(sea),cos(lat*atan(1.)/45.),0)
  LAT = ind(lat.ge.s1.and.lat.le.s2)
  LON = ind(lon.ge.s3.and.lon.le.s4)

  ;; read data, and perform JF-mean and spatial average
  M = m2-m1+1
  Y = y2-y1+1
  lhf = new((/M,Y/),float)    ;; surface latent heat flux
  wspd = new((/M,Y/),float)   ;; surface scalar wind speed
  u10  = new((/M,Y/),float)   ;; surface zonal wind
  v10  = new((/M,Y/),float)   ;; surface meridional wind
  qs   = new((/M,Y/),float)   ;; saturation specific humidity at SST
  q2   = new((/M,Y/),float)   ;; surface specific humidity
  do i = 0, M-1
    mxxx = "m"+sprinti("%0.3i",i+m1)
    do j = 0, Y-1
      yyyy = sprinti("%0.4i",j+y1)
      system("echo read m="+mxxx+", y="+yyyy)
      ifile = dir_d4pdf+expr+"_"+mxxx+"/"+yyyy+"##/sfc_avr_mon_"+expr+"_"+mxxx+"_"+yyyy+"##.grib"
      lhf(i,j)  = jf_mean(ifile,j+y1,"FLLH_GDS4_ISBL",sea,LAT,LON)
      wspd(i,j) = jf_mean(ifile,j+y1,"WIND_GDS4_ISBL",sea,LAT,LON)
      u10(i,j)  = jf_mean(ifile,j+y1,"UA_GDS4_ISBL",sea,LAT,LON)
      v10(i,j)  = jf_mean(ifile,j+y1,"VA_GDS4_ISBL",sea,LAT,LON)
      qs(i,j)   = jf_mean(ifile,j+y1,"qs",sea,LAT,LON)
      q2(i,j)   = jf_mean(ifile,j+y1,"QA_GDS4_ISBL",sea,LAT,LON)
    end do
  end do
  
  VARS = [/lhf,wspd,u10,v10,qs,q2/]
  var_names = (/"lhf","wspd","u10","v10","qs","q2"/)
  NV = dimsizes(var_names)

  ;; output
  system("rm -f lhf_area_"+expr+".nc")
  o = addfile("lhf_area_"+expr+".nc","c")
  VARS[0]!0 = "run"
  VARS[0]!1 = "year"
  VARS[0]&run  = m1+ispan(0,M-1,1)
  VARS[0]&year = y1+ispan(0,Y-1,1)
  o->$var_names(0)$ = VARS[0]
  do i = 1, NV-1
    copy_VarCoords(VARS[0],VARS[i])
    o->$var_names(i)$ = VARS[i]
  end do

end
