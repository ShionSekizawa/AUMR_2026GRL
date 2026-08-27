;;; File: make_ausm_rain.obs.ncl
;;; Description: Calculate AUMR (JF-mean rainfall over northern Australia) for observation (AGCD)
;;;
;;; Input: AGCD data
;;; Output: ausm_rain_obs.nc
;;;
;;; Running:
;;; ncl make_ausm_rain.obs.ncl
;;;

begin

  y1 = 1951
  y2 = 2011
  yyyymm = yyyymm_time(y1,y2,"integer")
  
;;; Observation
  in = addfile("${AGCD_dir_path}/v2/land_fixed_0.05.nc","r")  ;; open AGCD land-sea mask file
  wgt = short2flt(in->land)
  wgt = wgt*conform_dims(dimsizes(wgt),cos(tofloat(in->lat)*atan(1.)/45.),0)
  wgt({:-20},:) = 0.0  ;; weight array for spatial average
  ausm = new((/y2-y1+1/),float)
  do i = y1, y2
    yyyy = sprinti("%0.4i",i)
    ;; open AGCD rainfall data file
    in = addfile("${AGCD_dir_path}/v2/monthly/precip_total_r005_monthly_"+yyyy+"01.nc","r")
    a1 = wgt_areaave2(in->precip,wgt,0)
    in = addfile("${AGCD_dir_path}/v2/monthly/precip_total_r005_monthly_"+yyyy+"02.nc","r")
    a2 = wgt_areaave2(in->precip,wgt,0)
    d1 = days_in_month(i,1)
    d2 = days_in_month(i,2)
    ausm(i-y1) = (a1+a2)/(d1+d2)
  end do
  ausm!0 = "year"
  ausm&year = ispan(y1,y2,1)

  ;; output
  system("rm -f ausm_rain_obs.nc")
  o = addfile("ausm_rain_obs.nc","c")
  o->ausm = ausm

end