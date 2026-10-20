;;; File: cov.obs.ncl
;;; Description: Calculate statistics of AUMR and U10 in observations
;;;
;;; Input: sst_pc.obs.nc created by sst_pc.obs.ncl
;;;        ausm_rain_obs.nc created by make_ausm_rain.obs.ncl
;;;        JRA-3Q data
;;; Output: 
;;;
;;; Running:
;;; ncl cov.obs.ncl
;;;

begin

  in = addfile("sst_pc.obs.nc","r")
  pc1 = in->pc(0,:)
  in = addfile("ausm_rain_obs.nc","r")
  ausm = in->ausm({1951:2010})
  ausm_mean = avg(ausm)
  ausm = dtrend(ausm,False)
  ausm = dtrend_msg(pc1,ausm,False,False)
  ausm_std = stddev(ausm)
  ausm = ausm/ausm_std

  dum = new(60,float)
  lat = latGau(480,"lat","latitude","degrees_north")
  lon = lonGlobeF(960,"lon","longitude","degrees_east")
  gwt = latGauWgt(480,"lat","gaussian weights","")
  LAT = ind(lat.ge.-20.and.lat.le.-6)
  LON = ind(lon.ge.105.and.lon.le.150)
  in = addfile("${JRA-3Q_dir_path}/Const/TL479_land.grib2","r")
  wgt = 1. - in->LAND_P0_L1_GGA0(::-1,:)
  wgt = wgt*conform_dims(dimsizes(wgt),gwt,0)
  do yyyy = 1951, 2010
    print(yyyy+"")
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/anl_surf/anl_surf."+yyyy+"01.grb2","r")
    u1 = in->UGRD_P8_L103_GGA0_avg(::-1,:)
    in = addfile("${JRA-3Q_dir_path}/Hist/Monthly/anl_surf/anl_surf."+yyyy+"02.grb2","r")
    u2 = in->UGRD_P8_L103_GGA0_avg(::-1,:)
    d1 = days_in_month(yyyy,1)
    d2 = days_in_month(yyyy,2)
    dum(yyyy-1951) = wgt_areaave2((d1*u1(LAT,LON)+d2*u2(LAT,LON))/(d1+d2),wgt(LAT,LON),0)
  end do
  u10_mean = avg(dum)
  dum = dtrend(dum,False)
  dum = dtrend_msg(pc1,dum,False,False)
  u10_reg = regCoef(ausm,dum)
    
  print("AGCD:  aumr_mean="+sprintf("%3.2f",ausm_mean)+"  aumr_stddev="+sprintf("%3.2f",ausm_std))
  print("JRA-3Q:  u10_mean="+sprintf("%3.2f",u10_mean)+"  u10_reg="+sprintf("%3.2f",u10_reg)+" ratio:"+sprintf("%4.3f",u10_reg/u10_mean))

end