;;; File: cov.ncl
;;; Description: Calculate statistics of AUMR and U10 in d4PDF historical simulations
;;;
;;; Input: ausm_rain_HPB.nc created by d4PDF/make_ausm_rain.ncl
;;;        lhf_area_HPB.nc created by d4PDF/make_lhf_area.ncl
;;; Output: 
;;;
;;; Running:
;;; ncl cov.ncl
;;;

begin
  ;; get information on the d4PDF experiments
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(0,0)
  m1 = toint(data_list(0,1))
  m2 = toint(data_list(0,2))
  y1 = toint(data_list(0,3))
  y2 = toint(data_list(0,4))

;;; d4PDF data
  in1 = addfile("./lhf_area_"+expr(0)+".nc","r")
  in2 = addfile("./ausm_rain_"+expr(0)+".nc","r")
  u10_mean = avg(in1->u10)
  ausm_mean = avg(in2->ausm)
  dum2 = ndtooned(dim_rmvmean_n(in1->u10,0))
  dum = ndtooned(dim_rmvmean_n(in2->ausm,0))
  ausm_std = stddev(dum)
  dum = dum/stddev(dum)
  u10_reg = regCoef(dum,dum2)

  print("d4PDF HIST:  aumr_mean="+sprintf("%3.2f",ausm_mean)+"  aumr_stddev="+sprintf("%3.2f",ausm_std))
  print("d4PDF HIST:  u10_mean="+sprintf("%3.2f",u10_mean)+"  u10_reg="+sprintf("%3.2f",u10_reg)+" ratio:"+sprintf("%4.3f",u10_reg/u10_mean))

end