;;; File: gwlev.ncl
;;; Description: Calculate 60-year mean warming level relative to the pre-industrial level based on GMST in CMIP6 models
;;;
;;; Input: gmst_piControl.nc created by CMIP6/gmst_pi.ncl
;;;        gmst_historical.nc, gmst_ssp245.nc, gmst_ssp585.nc created by CMIP6/gmst.ncl
;;; Output: gwlev.nc
;;;
;;; Running:
;;; ncl gwlev.ncl
;;;

begin

  ;; read GMST time series and agerage over the analysis period
  in = addfile("gmst_piControl.nc","r")
  gmst = dim_avg(in->gmst)
  M = dimsizes(gmst)

  in = addfile("gmst_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  gmst0 = dim_avg(in->gmst(:,id))-gmst  ;;; relative to the pre-industrial level
  delete([/yyyy,id/])

  in = addfile("gmst_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  gmst2 = dim_avg(in->gmst(:,id))-gmst  ;;; relative to the pre-industrial level
  delete([/yyyy,id/])

  in = addfile("gmst_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  gmst4 = dim_avg(in->gmst(:,id))-gmst  ;;; relative to the pre-industrial level
  delete([/yyyy,id/])

  gmst0!0 = "model"
  gmst0&model = ispan(1,M,1)
  copy_VarCoords(gmst0,gmst2)
  copy_VarCoords(gmst0,gmst4)
  
  ;;; output
  system("rm -f gwlev.nc")
  o = addfile("gwlev.nc","c")
  o->gwlev_historical = gmst0
  o->gwlev_ssp245     = gmst2
  o->gwlev_ssp585     = gmst4
end