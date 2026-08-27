;;; File: make_figS2_tmpfile_obs.ncl
;;; Description: make a temporary (intermediate) file for creating Figure S2
;;;
;;; Input: AGCD data, COBE-SST2 data, JRA-3Q data
;;;        ausm_rain_obs.nc created by make_ausm_rain.obs.ncl
;;; Output: figS2_tmpfile_obs.nc
;;;
;;; Running:
;;; ncl make_figS2_tmpfile_obs.ncl
;;;

begin

  y1 = 1951
  y2 = 2011
  year = ispan(y1,y2,1)
  yyyymm = yyyymm_time(y1,y2,"integer")
  
;;; Observation

;; AGCD rainfall stddev
  in = addfile("${AGCD_dir_path}/v2/land_fixed_0.05.nc","r")  ;; open AGCD land-sea mask file
  land = short2flt(in->land)
  lon = in->lon
  lat = in->lat
  dum  = new((/y2-y1+1,691,886/),float)
  do i = y1, y2
    yyyy = sprinti("%0.4i",i)
    ;; open AGCD rainfall data file
    in = addfile("${AGCD_dir_path}/v2/monthly/precip_total_r005_monthly_"+yyyy+"01.nc","r")
    a1 = in->precip
    in = addfile("${AGCD_dir_path}/v2/monthly/precip_total_r005_monthly_"+yyyy+"02.nc","r")
    a2 = in->precip
    d1 = days_in_month(i,1)
    d2 = days_in_month(i,2)
    dum(i-y1,:,:) = (a1+a2)/(d1+d2)
  end do
  dum = dtrend_n(dum,False,0)
  obs_std = dim_stddev_n(dum,0)
  obs_std = where(land.lt.0.5,obs_std@_FillValue,obs_std)
  obs_std!0 = "lat_0"
  obs_std!1 = "lon_0"
  obs_std&lat_0 = lat
  obs_std&lon_0 = lon
  delete([/dum,lat,lon,a1,a2,yyyy/])

;; AUMR
  in = addfile("ausm_rain_obs.nc","r")
  ausm = in->ausm
  ausm = dtrend_n(ausm,False,0)
  ausm = ausm/stddev(ausm)
  
;; Correlations with SST
  dum = new((/y2-y1+1,180,360/),float)
  do yyyy = y1, y2
    in = addfile("${COBE_SST2_dir_path}/mon/sst-glb."+yyyy+".dat.grb","r")
    a1 = in->WTMP_GDS0_DBSL(0,:,:)
    a2 = in->WTMP_GDS0_DBSL(1,:,:)
    d1 = days_in_month(yyyy,1)
    d2 = days_in_month(yyyy,2)
    dum(yyyy-y1,:,:) = (a1*d1+a2*d2)/(d1+d2)
  end do
  dum = dtrend_n(dum,False,0)
  obs_cor = escorc_n(ausm,dum,0,0)
  copy_VarCoords(a1,obs_cor)
  obs_cor!0 = "lat_1"
  obs_cor!1 = "lon_1"
  delete([/dum,a1,a2/])

;; Regression anomalies of wind and evaporation
  dum1 = new((/y2-y1+1,480,960/),float)
  dum2 = new((/y2-y1+1,480,960/),float)
  dum3 = new((/y2-y1+1,480,960/),float)
  do yyyy = y1, y2
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/anl_surf/anl_surf."+yyyy+"01.grb2","r")
    u1 = in->UGRD_P8_L103_GGA0_avg
    v1 = in->VGRD_P8_L103_GGA0_avg
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/fcst_phy2m/fcst_phy2m."+yyyy+"01.grb2","r")
    e1 = in->LHTFL_P8_L1_GGA0_avg
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/anl_surf/anl_surf."+yyyy+"02.grb2","r")
    u2 = in->UGRD_P8_L103_GGA0_avg
    v2 = in->VGRD_P8_L103_GGA0_avg
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/fcst_phy2m/fcst_phy2m."+yyyy+"02.grb2","r")
    e2 = in->LHTFL_P8_L1_GGA0_avg
    d1 = days_in_month(yyyy,1)
    d2 = days_in_month(yyyy,2)
    dum1(yyyy-y1,:,:) = (u1*d1+u2*d2)/(d1+d2)
    dum2(yyyy-y1,:,:) = (v1*d1+v2*d2)/(d1+d2)
    dum3(yyyy-y1,:,:) = (e1*d1+e2*d2)/(d1+d2)
  end do
  dum1 = dtrend_n(dum1,False,0)
  dum2 = dtrend_n(dum2,False,0)
  dum3 = dtrend_n(dum3,False,0)
  obs_u  = regCoef_n(ausm,dum1,0,0)
  obs_v  = regCoef_n(ausm,dum2,0,0)
  obs_lf = regCoef_n(ausm,dum3,0,0)
  copy_VarCoords(u1,obs_u)
  obs_u!0 = "lat_2"
  obs_u!1 = "lon_2"
  copy_VarCoords(obs_u,obs_v)
  copy_VarCoords(obs_u,obs_lf)
  delete([/dum1,dum2,dum3,u1,u2,v1,v2,e1,e2,ausm,land,yyyy/])

  system("rm -f figS2_tmpfile_obs.nc")
  o = addfile("figS2_tmpfile_obs.nc","c")
  o->obs_std = obs_std
  o->obs_cor = obs_cor
  o->obs_u = obs_u
  o->obs_v = obs_v
  o->obs_lf = obs_lf
  
end