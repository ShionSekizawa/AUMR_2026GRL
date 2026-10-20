;;; File: gwlev.ncl
;;; Description: Calculate 60-year mean warming level relative to the pre-industrial level based on GMST in CMIP6 models
;;;
;;; Input: gmst_piControl.nc created by CMIP6/gmst_pi.ncl
;;;        gmst_historical_mulens.nc and gmst_ssp585_mulens.nc created by CMIP6/mulens/gmst.ncl
;;; Output: gwlev_mulens.nc
;;;
;;; Running:
;;; ncl gwlev.ncl
;;;

begin

  M = 13
  ;; read GMST time series and agerage over the analysis period
  list0 = "../ModelMemList.txt"
  models = asciiread(list0,toint(systemfunc("wc -l "+list0)),"string")
  dum2 = new(toint(systemfunc("wc -l "+list0)),string)
  do i = 0, dimsizes(dum2)-1
    dum = str_split(models(i)," ")
    dum2(i) = dum(1)
    delete(dum)
  end do
  list0 = "ModelMemList_mulens.txt"
  models_mulens = asciiread(list0,toint(systemfunc("wc -l "+list0)),"string")
  model_id = new(M,integer)
  do i = 0, M-1
    dum1 = str_split(models_mulens(i*5)," ")
    model_id(i) = ind(dum2.eq.dum1(1))
    delete(dum1)
  end do
  print(models(model_id))
  
  in = addfile("../gmst_piControl.nc","r")
  gmst = conform_dims((/M,5/),dim_avg(in->gmst(model_id,:)),0)

  in = addfile("gmst_historical_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  gmst0 = dim_avg(in->gmst(:,:,id))-gmst  ;;; relative to the pre-industrial level
  delete([/yyyy,id/])

  in = addfile("gmst_ssp585_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  gmst4 = dim_avg(in->gmst(:,:,id))-gmst  ;;; relative to the pre-industrial level
  delete([/yyyy,id/])

  gmst0!0 = "model"
  gmst0&model = ispan(1,M,1)
  gmst0!1 = "ens"
  gmst0&ens = ispan(1,5,1)
  copy_VarCoords(gmst0,gmst4)
  
  ;;; output
  system("rm -f gwlev_mulens.nc")
  o = addfile("gwlev_mulens.nc","c")
  o->gwlev_historical = gmst0
  o->gwlev_ssp585     = gmst4
end