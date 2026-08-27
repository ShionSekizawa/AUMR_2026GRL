;;; File: sst_pc.ncl
;;; Description: Perform EOF (PCA) analysis for SST in CMIP6 models
;;;
;;; Input: sst_historical.nc, sst_ssp245.nc, sst_ssp585.nc created by CMIP6/sst.ncl
;;; Output: sst_pc_historical.nc, sst_pc_ssp245.nc, sst_pc_ssp585.nc
;;;
;;; Running:
;;; ncl sst_pc.ncl
;;;

;;; procedure for EOF (PCA) analysis
procedure sst_eof(sst[*][*][*][*]:numeric,ofile[1]:string)
local dim, wgt, j, i, opt, neof, o, dum1, dum2, LON
begin
  neof = 5
  dim = dimsizes(sst)
  wgt = sqrt(cos(sst&$sst!2$*atan(1.)/45.))
  LON = ind(sst&$sst!3$.ge.190.and.sst&$sst!3$.le.270)
  sst = dim_rmvmean_n(sst,1)
  sst = dtrend_n(sst,False,1)
  do j = 0, dimsizes(wgt)-1
    sst(:,:,j,:) = sst(:,:,j,:) * wgt(j)
  end do
  opt   = True
  opt@jopt = 0
  eof = new((/dim(0),neof,dim(2),dim(3)/),typeof(sst))
  pc  = new((/dim(0),neof,dim(1)/),typeof(sst))
  pcvar = new((/dim(0),neof/),typeof(sst))
  do i = 0, dim(0)-1
    dum1 = eofunc_n(sst(i,:,:,:),neof,opt,0)
    dum2 = eofunc_ts_n(sst(i,:,:,:),dum1,False,0)
    pcvar(i,:) = dum1@pcvar
    do j = 0, neof-1
      dum1(j,:,:) = dum1(j,:,:)*sqrt(dum1@eval(j))
      dum2(j,:)  = dum2(j,:)/sqrt(dum1@eval(j))
      if (avg(dum1(j,:,LON)).lt.0) then
        dum1(j,:,:) = -dum1(j,:,:)
        dum2(j,:) = -dum2(j,:)
      end if
    end do
    do j = 0, dimsizes(wgt)-1
      dum1(:,j,:) = dum1(:,j,:) / wgt(j)
    end do
    eof(i,:,:,:) = dum1
    pc(i,:,:) = dum2
    delete([/dum1,dum2/])
  end do
  ;; set coordinates
  eof!0 = "model"
  eof!1 = "eof_mode"
  eof!2 = sst!2
  eof!3 = sst!3
  eof&model = ispan(1,dim(0),1)
  eof&eof_mode = ispan(1,neof,1)
  eof&$eof!2$ = sst&$sst!2$
  eof&$eof!3$ = sst&$sst!3$
  pc!0 = "model"
  pc!1 = "eof_mode"
  pc!2 = sst!1
  pc&model = ispan(1,dim(0),1)
  pc&eof_mode = ispan(1,neof,1)
  pc&$pc!2$ = sst&$sst!1$
  copy_VarCoords(pc(:,:,0),pcvar)
  ;; output
  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->eof = eof
  o->pc  = pc
  o->pcvar = pcvar
end

begin
  ;; number of models
  M = 26

  ;; read historical simulation data
  in = addfile("sst_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  Y  = dimsizes(id) 
  lat = in->lat
  lon = in->lon
  ;; EOF (PCA) domain
  LAT = ind(lat.ge.-15.and.lat.le.15)
  LON = ind(lon.ge.0.and.lon.le.360)
  lat := lat(LAT)
  lon := lon(LON)
  NY = dimsizes(lat)
  NX = dimsizes(lon)
  sst = in->sst(:,id,LAT,LON)
  ;; EOF (PCA) and output
  sst_eof(sst,"sst_pc_historical.nc")
  delete([/yyyy,id,sst/])

  ;; read SSP2-4.5 simulation data
  in = addfile("sst_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  sst = in->sst(:,id,LAT,LON)
  ;; EOF (PCA) and output
  sst_eof(sst,"sst_pc_ssp245.nc")
  delete([/yyyy,id,sst/])

  ;; read SSP5-8.5 simulation data
  in = addfile("sst_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  sst = in->sst(:,id,LAT,LON)
  ;; EOF (PCA) and output
  sst_eof(sst,"sst_pc_ssp585.nc")
  delete([/yyyy,id,sst/])
  
end