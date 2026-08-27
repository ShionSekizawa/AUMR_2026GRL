;;; File: make_figS2_tmpfile_d4pdf.ncl
;;; Description: make a temporary (intermediate) file for creating Figure S2
;;;
;;; Input: d4PDF data, ausm_rain_HPB.nc created by make_ausm_rain.ncl
;;; Output: figS2_tmpfile_d4pdf.nc
;;;
;;; Running:
;;; ncl make_figS2_tmpfile_d4pdf.ncl
;;;

procedure dtrend_2d(x:numeric)
local x0, i, dim
begin
  x0 = dim_avg_n(x,0)
  x = dim_rmvmean_n(x,0)
  x0 = dtrend_n(x0,False,0)
  dim = dimsizes(x)
  if (dimsizes(dim).eq.2) then
    do i = 0, dim(0)-1
      x(i,:) = x(i,:) + x0
    end do
  end if
  if (dimsizes(dim).eq.4) then
    do i = 0, dim(0)-1
      x(i,:,:,:) = x(i,:,:,:) + x0
    end do
  end if
end
  
begin

  ;;; directory path referring to the d4PDF
  dir_d4pdf = "${d4PDF_dir_path}/Global/"

  ;;; get information for the HIST experiment (I=0, fixed)
  I = 0
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(I,0)
  m1 = toint(data_list(I,1))
  m2 = toint(data_list(I,2))
  y1 = toint(data_list(I,3))
  y2 = toint(data_list(I,4))
  
  setfileoption("bin","ReadByteOrder","BigEndian")
  ;; read d4PDF land-sea mask file
  land = fbindirread("${d4PDF_dir_path}/TopogRatiol_gsmuv_TL319.gd",1,(/320,640/),"float")
  land@_FillValue =  -9.99E33
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
  lat@units = "degrees_N"
  lon@units = "degrees_E"        

  ;; read d4PDF data
  M = m2-m1+1
  Y = y2-y1+1
  dum  = new((/M,Y,320,640/),float)
  dum0 = new((/M,Y,320,640/),float)
  dum1 = new((/M,Y,320,640/),float)
  dum2 = new((/M,Y,320,640/),float)
  dum3 = new((/M,Y,320,640/),float)
  do i = 0, M-1
    mxxx = "m"+sprinti("%0.3i",i+m1)
    do j = 0, Y-1
      yyyy = sprinti("%0.4i",j+y1)
      system("echo read m="+mxxx+", y="+yyyy)
      in = addfile(dir_d4pdf+expr+"_"+mxxx+"/"+yyyy+"01/sfc_avr_mon_"+expr+"_"+mxxx+"_"+yyyy+"01.grib","r")
      r1 = in->PRECIPI_GDS4_ISBL
      a1 = in->TGEF_GDS4_ISBL
      u1 = in->UA_GDS4_ISBL
      v1 = in->VA_GDS4_ISBL
      e1 = in->FLLH_GDS4_ISBL
      in = addfile(dir_d4pdf+expr+"_"+mxxx+"/"+yyyy+"02/sfc_avr_mon_"+expr+"_"+mxxx+"_"+yyyy+"02.grib","r")
      r2 = in->PRECIPI_GDS4_ISBL
      a2 = in->TGEF_GDS4_ISBL
      u2 = in->UA_GDS4_ISBL
      v2 = in->VA_GDS4_ISBL
      e2 = in->FLLH_GDS4_ISBL
      d1 = days_in_month(j+y1,1)
      d2 = days_in_month(j+y1,2)
      dum(i,j,:,:)  = 60*60*24.*(r1*d1+r2*d2)/(d1+d2)
      dum0(i,j,:,:) = where(land.ge.0.5,dum0@_FillValue,(a1*d1+a2*d2)/(d1+d2))
      dum1(i,j,:,:) = (u1*d1+u2*d2)/(d1+d2)
      dum2(i,j,:,:) = (v1*d1+v2*d2)/(d1+d2)
      dum3(i,j,:,:) = (e1*d1+e2*d2)/(d1+d2)
      delete([/a1,a2,u1,u2,v1,v2,e1,e2/])
    end do
  end do
  dtrend_2d(dum)
  dtrend_2d(dum0)
  dtrend_2d(dum1)
  dtrend_2d(dum2)
  dtrend_2d(dum3)
  
;; AUMR
  in = addfile("ausm_rain_HPB.nc","r")
  ausm = in->ausm
  dtrend_2d(ausm)
  ausm = ausm/stddev(ausm)

  d4pdf_std = dim_stddev_n(reshape(dum,(/M*Y,320,640/)),0)
  d4pdf_cor = escorc_n(ndtooned(ausm),reshape(dum0,(/M*Y,320,640/)),0,0)
  d4pdf_u   = regCoef_n(ndtooned(ausm),reshape(dum1,(/M*Y,320,640/)),0,0)
  d4pdf_v   = regCoef_n(ndtooned(ausm),reshape(dum2,(/M*Y,320,640/)),0,0)
  d4pdf_lf  = regCoef_n(ndtooned(ausm),reshape(dum3,(/M*Y,320,640/)),0,0)
  d4pdf_cor!0 = "lat"
  d4pdf_cor!1 = "lon"
  d4pdf_cor&lat = lat
  d4pdf_cor&lon = lon
  copy_VarCoords(d4pdf_cor,d4pdf_std)
  copy_VarCoords(d4pdf_cor,d4pdf_u)
  copy_VarCoords(d4pdf_cor,d4pdf_v)
  copy_VarCoords(d4pdf_cor,d4pdf_lf)

  system("rm -f figS2_tmpfile_d4pdf.nc")
  o = addfile("figS2_tmpfile_d4pdf.nc","c")
  o->d4pdf_std = d4pdf_std
  o->d4pdf_cor = d4pdf_cor
  o->d4pdf_u = d4pdf_u
  o->d4pdf_v = d4pdf_v
  o->d4pdf_lf = d4pdf_lf

  dum = dim_rmvmean_n(dum,0)
  dum0 = dim_rmvmean_n(dum0,0)
  dum1 = dim_rmvmean_n(dum1,0)
  dum2 = dim_rmvmean_n(dum2,0)
  dum3 = dim_rmvmean_n(dum3,0)
  d4pdf_std = dim_stddev_n(reshape(dum,(/M*Y,320,640/)),0)
  d4pdf_cor = escorc_n(ndtooned(ausm),reshape(dum0,(/M*Y,320,640/)),0,0)
  d4pdf_u   = regCoef_n(ndtooned(ausm),reshape(dum1,(/M*Y,320,640/)),0,0)
  d4pdf_v   = regCoef_n(ndtooned(ausm),reshape(dum2,(/M*Y,320,640/)),0,0)
  d4pdf_lf  = regCoef_n(ndtooned(ausm),reshape(dum3,(/M*Y,320,640/)),0,0)
  d4pdf_cor!0 = "lat"
  d4pdf_cor!1 = "lon"
  d4pdf_cor&lat = lat
  d4pdf_cor&lon = lon
  copy_VarCoords(d4pdf_cor,d4pdf_std)
  copy_VarCoords(d4pdf_cor,d4pdf_u)
  copy_VarCoords(d4pdf_cor,d4pdf_v)
  copy_VarCoords(d4pdf_cor,d4pdf_lf)
  o->d4pdf_std_int = d4pdf_std
  o->d4pdf_cor_int = d4pdf_cor
  o->d4pdf_u_int   = d4pdf_u
  o->d4pdf_v_int   = d4pdf_v
  o->d4pdf_lf_int  = d4pdf_lf
  
end