;;; File: sst_pc.ncl
;;; Description: Perform EOF (PCA) analysis for SST in CMIP6 models
;;;
;;; Input: sst_historical_mulens.nc and sst_ssp585_mulens.nc created by CMIP6/mulens/sst.ncl
;;; Output: sst_pc_historical_mulens.nc, sst_pc_ssp585_mulens.nc
;;;
;;; Running:
;;; ncl sst_pc.ncl
;;;

;;; procedure for EOF (PCA) analysis
procedure sst_eof(sst[*][*][*][*][*]:numeric,ofile[1]:string)
local dim, wgt, j, i, opt, neof, o, dum1, dum2, LON, sstin, trend, x
begin
  neof = 5
  dim = dimsizes(sst)
  wgt = sqrt(cos(sst&$sst!3$*atan(1.)/45.))
  LON = ind(sst&$sst!4$.ge.190.and.sst&$sst!4$.le.270)
  x = ispan(1,dim(2),1)
  do i = 0, dim(0)-1
    sst(i,:,:,:,:) = dim_rmvmean_n(sst(i,:,:,:,:),(/0,1/))
    trend = conform_dims(dim(2:4),regCoef_n(x,dim_avg_n(sst(i,:,:,:,:),0),0,0),(/1,2/))
    do j = 0, dim(2)-1
      trend(j,:,:) = trend(j,:,:)*(x(j)-avg(x))
    end do
    do j = 0, 4
      sst(i,j,:,:,:) = sst(i,j,:,:,:) - trend
    end do
    delete(trend)
  end do
  do j = 0, dimsizes(wgt)-1
    sst(:,:,:,j,:) = sst(:,:,:,j,:) * wgt(j)
  end do
  opt   = True
  opt@jopt = 0
  eof = new((/dim(0),neof,dim(3),dim(4)/),typeof(sst))
  pc  = new((/dim(0),neof,dim(1),dim(2)/),typeof(sst))
  pcvar = new((/dim(0),neof/),typeof(sst))
  do i = 0, dim(0)-1
    sstin = reshape(sst(i,:,:,:,:),(/dim(1)*dim(2),dim(3),dim(4)/))
    dum1 = eofunc_n(sstin,neof,opt,0)
    dum2 = eofunc_ts_n(sstin,dum1,False,0)
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
    pc(i,:,:,:) = reshape(dum2,(/neof,dim(1),dim(2)/))
    delete([/dum1,dum2,sstin/])
  end do
  ;; set coordinates
  eof!0 = "model"
  eof!1 = "eof_mode"
  eof!2 = sst!3
  eof!3 = sst!4
  eof&model = ispan(1,dim(0),1)
  eof&eof_mode = ispan(1,neof,1)
  eof&$eof!2$ = sst&$sst!3$
  eof&$eof!3$ = sst&$sst!4$
  pc!0 = "model"
  pc!1 = "eof_mode"
  pc!2 = sst!1
  pc!3 = sst!2
  pc&model = ispan(1,dim(0),1)
  pc&eof_mode = ispan(1,neof,1)
  pc&$pc!2$ = sst&$sst!1$
  pc&$pc!3$ = sst&$sst!2$
  copy_VarCoords(pc(:,:,0,0),pcvar)
  ;; output
  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->eof = eof
  o->pc  = pc
  o->pcvar = pcvar
end

begin
  ;; number of models
  M = 13

  ;; read historical simulation data
  in = addfile("sst_historical_mulens.nc","r")
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
  sst = in->sst(:,:,id,LAT,LON)
  ;; EOF (PCA) and output
  sst_eof(sst,"sst_pc_historical_mulens.nc")
  delete([/yyyy,id,sst/])

  ;; read SSP5-8.5 simulation data
  in = addfile("sst_ssp585_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  sst = in->sst(:,:,id,LAT,LON)
  ;; EOF (PCA) and output
  sst_eof(sst,"sst_pc_ssp585_mulens.nc")
  delete([/yyyy,id,sst/])
  
end