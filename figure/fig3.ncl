;;; File: fig3.ncl
;;; Description: create Figure 3
;;;
;;; Input: ausm_rain_${experiment_type_name}.nc created by d4PDF/make_ausm_rain.ncl
;;;        lhf_area_${experiment_type_name}.nc created by d4PDF/make_lhf_area.ncl
;;;        ausm_variance.nc created by d4PDF/derive_ausm_variance.ncl
;;;        sst_pc_historical.nc, sst_pc_ssp245.nc and sst_pc_ssp585.nc created by CMIP6/sst_pc.ncl
;;;        ausm_historical.nc, ausm_ssp245.nc and ausm_ssp585.nc created by CMIP6/ausm.ncl
;;;        u10_historical.nc, u10_ssp245.nc and u10_ssp585.nc created by CMIP6/u10.ncl
;;;        u10_cov.ncl created by CMIP6/cov.ncl
;;;        gwlev.nc created by CMIP6/gwlev.ncl
;;;        cmip6/ModelMemList.txt
;;; Output: fig3.eps
;;;
;;; Running:
;;; ncl fig3.ncl
;;;

;;; function for converting 1d array into diagonal 2d array for marker plotting
function diag(x[*]:numeric)
local N, i
begin
  N = dimsizes(x)
  y = new((/N,N/),typeof(x))
  do i = 0, N-1
    y(i,i) = x(i)
  end do
  return(y)
end

;;; function for calculating anomalies (for 2d array)
function anom2(x[*][*]:numeric,clim_opt[1]:integer)
local dim, n
begin
  dim = dimsizes(x)
  n   = dim(0)
  if (clim_opt.eq.6) then
    do i = 0, clim_opt-1
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:) = \
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:) - avg(x(i*n/clim_opt:(i+1)*n/clim_opt-1,:))
    end do
  else
    x(:,:) = x(:,:) - avg(x(:,:))
  end if
  x = dim_rmvmean_n(x,0)
  return(ndtooned(x))
end

;;; function for calculating anomalies (for 3d array)
function anom3(x[*][*][*]:numeric,clim_opt[1]:integer)
local dim, n
begin
  dim = dimsizes(x)
  n   = dim(0)
  if (clim_opt.eq.6) then
    do i = 0, clim_opt-1
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:,:) = dim_rmvmean_n(x(i*n/clim_opt:(i+1)*n/clim_opt-1,:,:),(/0,1/))
    end do
  else
    x = dim_rmvmean_n(x,(/0,1/))
  end if
  x = dim_rmvmean_n(x,0)
  return(reshape(x,(/dim(0)*dim(1),dim(2)/)))
end

;;; function for calculating dSST-induced correlation (d4PDF)
function dSST_cor(xin[12]:numeric,yin[12]:numeric)
local i, x, y, cor, cor2, cor4, reg, reg2, reg4, id, cor_tmp, reg_tmp, stat, ret
begin
  x = xin
  y = yin
  x(:5) = dim_rmvmean(x(:5))
  x(6:) = dim_rmvmean(x(6:))
  y(:5) = dim_rmvmean(y(:5))
  y(6:) = dim_rmvmean(y(6:))
  cor = escorc(x,y)
  reg = regCoef(x,y)
  cor2 = escorc(x(:5),y(:5))
  reg2 = regCoef(x(:5),y(:5))
  cor4 = escorc(x(6:),y(6:))
  reg4 = regCoef(x(6:),y(6:))
  print("dSST_cor stats")
  id = ispan(0,5,1)
  cor_tmp = new(6,float)
  reg_tmp = new(6,float)
  stat = False
  do i = 0, 5
    cor_tmp(i) = escorc(x(ind(id.ne.i)),y(ind(id.ne.i)))
    reg_tmp(i) = regCoef(x(ind(id.ne.i)),y(ind(id.ne.i)))
  end do
  print("   2K: "+cor2+" ("+min(cor_tmp)+" - "+max(cor_tmp)+")")
  print("       "+reg2+" ("+min(reg_tmp)+" - "+max(reg_tmp)+")")
  print("     (R2 = "+cor2*cor2+")")
  if (min(cor_tmp)*max(cor_tmp).gt.0) then
    stat = True
  end if
  do i = 0, 5
    cor_tmp(i) = escorc(x(ind(id.ne.i)+6),y(ind(id.ne.i)+6))
    reg_tmp(i) = regCoef(x(ind(id.ne.i)+6),y(ind(id.ne.i)+6))
  end do
  print("   4K: "+cor4+" ("+min(cor_tmp)+" - "+max(cor_tmp)+")")
  print("       "+reg4+" ("+min(reg_tmp)+" - "+max(reg_tmp)+")")
  print("     (R2 = "+cor4*cor4+")")
  if (stat.and.min(cor_tmp)*max(cor_tmp).gt.0) then
    stat = True
  else
    stat = False
  end if
  print("   2K/4K: "+cor)
  print("          "+reg)
  print("       (R2 = "+cor*cor+")")
  print("   robustness test: "+stat)
  ret = (/cor,reg/)
  ret@cor2 = cor2
  ret@cor4 = cor4
  ret@reg2 = reg2
  ret@reg4 = reg4
  ret@stat = stat
  return(ret)
end

;;; function for calculating dSST-induced correlation (CMIP6)
function dSST_cor_c6(x1[*]:numeric,x2[*]:numeric,y1[*]:numeric,y2[*]:numeric)
local x, y, cor, cor2, cor4, reg, reg2, reg4, ret, N, t, r, z, se, zlow, zhi, low, hi
begin
  x = array_append_record(dim_rmvmean(x1),dim_rmvmean(x2),0)
  y = array_append_record(dim_rmvmean(y1),dim_rmvmean(y2),0)
  N = dimsizes(x)
  t = cdft_t(0.975,N-2)
  r = t/sqrt(N-2+t^2)
  cor  = escorc(x,y)
  reg  = regCoef(x,y)
  ret  = (/cor,reg/)
  cor2 = escorc(x1,y1)
  cor4 = escorc(x2,y2)
  reg2 = regCoef(x1,y1)
  reg4 = regCoef(x2,y2)
  ret@cor2 = cor2
  ret@cor4 = cor4
  ret@reg2 = reg2
  ret@reg4 = reg4
  ret@stat = abs(cor).gt.r
  print("dSST_cor_c6 stats")
  z  = 0.5*log((1+cor2)/(1-cor2))
  se = 1.0/sqrt(N/2-3)
  zlow = z - 1.96*se
  zhi  = z + 1.96*se
  low = (exp(2*zlow)-1)/(exp(2*zlow)+1)
  hi  = (exp(2*zhi)-1)/(exp(2*zhi)+1)  
  print("   245: "+cor2+" ("+low+" - "+hi+")")
  t = cdft_t(0.975,N/2-2)
  low = reg2 - reg2@rstd*t
  hi  = reg2 + reg2@rstd*t
  print("        "+reg2+" ("+low+" - "+hi+")")
  z  = 0.5*log((1+cor4)/(1-cor4))
  se = 1.0/sqrt(N/2-3)
  zlow = z - 1.96*se
  zhi  = z + 1.96*se
  low = (exp(2*zlow)-1)/(exp(2*zlow)+1)
  hi  = (exp(2*zhi)-1)/(exp(2*zhi)+1)
  print("         r2 = "+cor2*cor2)
  print("   585: "+cor4+" ("+low+" - "+hi+")")
  t = cdft_t(0.975,N/2-2)
  low = reg4 - reg4@rstd*t
  hi  = reg4 + reg4@rstd*t
  print("        "+reg4+" ("+low+" - "+hi+")")
  z  = 0.5*log((1+cor)/(1-cor))
  se = 1.0/sqrt(N-3)
  zlow = z - 1.96*se
  zhi  = z + 1.96*se
  low = (exp(2*zlow)-1)/(exp(2*zlow)+1)
  hi  = (exp(2*zhi)-1)/(exp(2*zhi)+1)
  print("         r2 = "+cor4*cor4)
  print("   245/585: "+cor+" ("+low+" - "+hi+")")
  t = cdft_t(0.975,N-2)
  low = reg - reg@rstd*t
  hi  = reg + reg@rstd*t
  print("            "+reg+" ("+low+" - "+hi+")")
  print("         r2 = "+cor*cor)
  print("  stat = "+ret@stat)
  return(ret)
end

function sst_index(sst)
local lat, mc, wp, wgt, ss, w2, e2, w3, e3, ret
begin
  lat = sst&lat
  lat!0 = "lat"
  lat&lat = lat
  w2 = 105
  e2 = 135
  w3 = 135
  e3 = 180
  ss  = 10.
  wgt = cos(lat({-ss:ss})*atan(1.)/45)
  if (dimsizes(dimsizes(sst)).eq.3) then
    mc  = wgt_areaave(sst(:,{-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst(:,{-ss:ss},{w3:e3}),wgt,1.,0)
  end if
  if (dimsizes(dimsizes(sst)).eq.2) then
    mc  = wgt_areaave(sst({-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst({-ss:ss},{w3:e3}),wgt,1.,0)
  end if
  ret = wp-mc
  return(ret)
end

begin
  ;; get information on the d4PDF experiments
  data_list = str_split_csv(asciiread("../d4PDF/data_list",14,"string")," ",0)
  expr = data_list(:,0)
  m1 = toint(data_list(:,1))
  m2 = toint(data_list(:,2))
  y1 = toint(data_list(:,3))
  y2 = toint(data_list(:,4))

;;; d4PDF data
  u10_mean = new(15,float)
  wspd_reg  = new(15,float)
  sst_mean = new(15,float)
  in1 = addfile("../d4PDF/lhf_area_"+expr(0)+".nc","r")
  in2 = addfile("../d4PDF/ausm_rain_"+expr(0)+".nc","r")
  in3 = addfile("../d4PDF/sst_clim_"+expr(0)+".nc","r")
  lat = in3->lat
  u10_mean(0) = avg(in1->u10)
  wspd_reg(0)  = regCoef_n(anom2(in1->u10,-1),anom2(in1->wspd,-1),0,0)
  sst_mean(0) = sst_index(in3->sst)
  u10_hist_ens = dim_avg_n(in1->u10,0)
  var_hist_int = dim_variance_n(in2->ausm,0)
  delete([/in1,in2,in3/])
  in1 = addfiles("../d4PDF/lhf_area_"+expr(2:7)+".nc","r")
  in2 = addfiles("../d4PDF/ausm_rain_"+expr(2:7)+".nc","r")
  in3 = addfiles("../d4PDF/sst_clim_"+expr(2:7)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in2,"cat")
  ListSetType(in3,"join")
  u10_mean(1) = avg(in1[:]->u10)
  wspd_reg(1) = regCoef_n(anom2(in1[:]->u10,6),anom2(in1[:]->wspd,6),0,0)
  sst_mean(1) = sst_index(dim_avg_n_Wrap(in3[:]->sst,0))
  in1 = addfiles("../d4PDF/lhf_area_"+expr(8:13)+".nc","r")
  in2 = addfiles("../d4PDF/ausm_rain_"+expr(8:13)+".nc","r")
  in3 = addfiles("../d4PDF/sst_clim_"+expr(8:13)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in2,"cat")
  ListSetType(in3,"join")
  u10_mean(2) = avg(in1[:]->u10)
  wspd_reg(2) = regCoef_n(anom2(in1[:]->u10,6),anom2(in1[:]->wspd,6),0,0)
  sst_mean(2) = sst_index(dim_avg_n_Wrap(in3[:]->sst,0))
  delete([/in1,in2,in3/])
  do i = 2, 13
    in1 = addfile("../d4PDF/lhf_area_"+expr(i)+".nc","r")
    in2 = addfile("../d4PDF/ausm_rain_"+expr(i)+".nc","r")
    in3 = addfile("../d4PDF/sst_clim_"+expr(i)+".nc","r")
    u10_mean(i+1) = avg(in1->u10)
    wspd_reg(i+1)  = regCoef_n(anom2(in1->u10,-1),anom2(in1->wspd,-1),0,0)
    sst_mean(i+1) = sst_index(in3->sst)
  end do
  u10_2K  = new((/6,(m2(2)-m1(2)+1)*(y2(2)-y1(2)+1)/),float)
  wspd_2K = new(dimsizes(u10_2K),float)
  do i = 2, 7
    in1 = addfile("../d4PDF/lhf_area_"+expr(i)+".nc","r")
    u10_2K(i-2,:) = ndtooned(in1->u10)
    wspd_2K(i-2,:) = ndtooned(in1->wspd)
  end do
  u10_4K  = new((/6,(m2(8)-m1(8)+1)*(y2(8)-y1(8)+1)/),float)
  wspd_4K = new(dimsizes(u10_4K),float)
  do i = 8, 13
    in1 = addfile("../d4PDF/lhf_area_"+expr(i)+".nc","r")
    u10_4K(i-8,:) = ndtooned(in1->u10)
    wspd_4K(i-8,:) = ndtooned(in1->wspd)
  end do
  delete([/in1/])

  in = addfile("../d4PDF/ausm_variance.nc","r")
  var = in->var(0,1,:)

;;; CMIP6 data
  dtrend_type = "linear"
  in = addfile("../CMIP6/sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  M = dimsizes(pc1(:,0))
  in = addfile("../CMIP6/ausm_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  dum = dim_rmvmean(in->ausm(:,id))
  in = addfile("../CMIP6/gmst_historical.nc","r")
  gmst = in->gmst(:,id)
  gmst_hist = dim_avg(gmst)
  if (dtrend_type.eq."linear") then
    dum = dtrend(dum,False)
  end if
  if (dtrend_type.eq."gmst") then
    do i = 0, M-1
      dum(i,:) = dtrend_msg(gmst(i,:),dum(i,:),False,False)
    end do
  end if
  if (dtrend_type.eq."gmst_mme") then
    dum = dtrend_msg(dim_avg_n(gmst,0),dum,False,False)
  end if
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  cov = asciiread("../CMIP6/u10_cov.txt",M,"float")
  M1 = ind(cov.lt.0.12)
  M2 = ind(cov.ge.0.12)
  print(dimsizes(M2))
  var_c6_hist = dim_variance(dum)
  in = addfile("../CMIP6/u10_historical.nc","r")
  u10_c6_hist = dim_avg(in->u10(:,id))
  delete([/yyyy,id,dum,pc1,gmst/])

  in = addfile("../CMIP6/sst_pc_ssp245.nc","r")
  pc1 = in->pc(:,0,:)
  in = addfile("../CMIP6/ausm_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  dum = in->ausm(:,id)
  dum = dim_rmvmean(in->ausm(:,id))
  in = addfile("../CMIP6/gmst_ssp245.nc","r")
  gmst = in->gmst(:,id)
  gmst_245 = dim_avg(gmst)
  if (dtrend_type.eq."linear") then
    dum = dtrend(dum,False)
  end if
  if (dtrend_type.eq."gmst") then
    do i = 0, M-1
      dum(i,:) = dtrend_msg(gmst(i,:),dum(i,:),False,False)
    end do
  end if
  if (dtrend_type.eq."gmst_mme") then
    dum = dtrend_msg(dim_avg_n(gmst,0),dum,False,False)
  end if
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  var_c6_245 = dim_variance(dum)
  in = addfile("../CMIP6/u10_ssp245.nc","r")
  u10_c6_245 = dim_avg(in->u10(:,id))
  delete([/yyyy,id,dum,pc1,gmst/])

  in = addfile("../CMIP6/sst_pc_ssp585.nc","r")
  pc1 = in->pc(:,0,:)
  in = addfile("../CMIP6/ausm_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  dum = in->ausm(:,id)
  dum = dim_rmvmean(in->ausm(:,id))
  in = addfile("../CMIP6/gmst_ssp585.nc","r")
  gmst = in->gmst(:,id)
  gmst_585 = dim_avg(gmst)
  if (dtrend_type.eq."linear") then
    dum = dtrend(dum,False)
  end if
  if (dtrend_type.eq."gmst") then
    do i = 0, M-1
      dum(i,:) = dtrend_msg(gmst(i,:),dum(i,:),False,False)
    end do
  end if
  if (dtrend_type.eq."gmst_mme") then
    dum = dtrend_msg(dim_avg_n(gmst,0),dum,False,False)
  end if
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  var_c6_585 = dim_variance(dum)
  in = addfile("../CMIP6/u10_ssp585.nc","r")
  u10_c6_585 = dim_avg(in->u10(:,id))
  delete([/yyyy,id,dum,pc1,gmst/])

  ;;; scaled by global warming level
  in = addfile("../CMIP6/gwlev.nc","r")
  gwlev_hist = in->gwlev_historical
  gwlev_245  = in->gwlev_ssp245
  gwlev_585  = in->gwlev_ssp585
  var_c6_245 = (var_c6_245-var_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  var_c6_585 = (var_c6_585-var_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  u10_c6_245 = (u10_c6_245-u10_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  u10_c6_585 = (u10_c6_585-u10_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
 
  m_tmp = asciiread("../CMIP6/ModelMemList.txt",M,"string")
  m_name = new(M,string)
  do i = 0, M-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do

  cmp = read_colormap_file("NCV_banded")
  NC = dimsizes(cmp(:,0))
  m_colors = cmp((ispan(1,2*M-1,2)*NC)/(2*M)-1,:)
  m_markers_2K = new(M,integer)
  m_markers_4K = new(M,integer)
  m_markers_2K(M1) = 2
  m_markers_4K(M1) = 3
  m_markers_2K(M2) = 4
  m_markers_4K(M2) = 6
  m_sizes = new(M,float)
  m_sizes(M1) = 0.005
  m_sizes(M2) = 0.007

;;; open workstation
  wks = gsn_open_wks("pdf","fig3")
;;; Panels (a-d)
  res1 = True
  res1@gsnDraw  = False
  res1@gsnFrame = False
  res1@tmBorderThicknessF  = 1.25
  res1@tmYLMajorThicknessF = 1.25
  res1@tmXBMajorThicknessF = 1.25
  res1@tmXBMajorLengthF = 0.005
  res1@tmYLMajorLengthF = 0.005
  res1@tmXBMinorLengthF = 0.002
  res1@tmYLMinorLengthF = 0.002
  res1@tmXBMajorOutwardLengthF = 0
  res1@tmYLMajorOutwardLengthF = 0
  res1@tmYLMinorOutwardLengthF = 0
  res1@tmXBMinorOutwardLengthF = 0
  res1@tmXBLabelFontHeightF = 0.008
  res1@tmYLLabelFontHeightF = 0.008
  res1@tiYAxisFontHeightF = 0.008
  res1@tiXAxisFontHeightF = 0.008
  res1@tmYROn = False
  res1@tmXTOn = False
  res1@xyMarkLineMode = "Markers"
  res1@xyMarkerThicknessF = 1.0
  res1@tmYLMode = "Manual"
  res1@tmXBMode = "Manual"
  colors = (/"black","orange","firebrick","lightgoldenrod2","pink3"/)
  res1@gsnYRefLine = 0
  res1@gsnXRefLine = 0
  
  res1@vpYF = 0.92
  x1 = 0.045
  x2 = 0.99
  a  = 0.20
  dx = (x2-x1)/(4+3*a)
  res1@vpWidthF  = dx
  res1@vpHeightF = dx
  res1@tiXAxisString = "10-m zonal wind [m/s]"
  res1@tiYAxisString = "10-m scalar wind speed [m/s]"
  res1@tmYLPrecision = 2
  res1@trXMinF = -2.0
  res1@trXMaxF =  6.0
  res1@tmXBTickSpacingF = 1.0
  res1@trYMinF = 2.0
  res1@trYMaxF = 8.0
  res1@tmYLTickSpacingF = 1.0
  res1@vpXF = x1+0*(1+a)*dx
  res1@xyMarkerColors := colors(3)
  res1@xyMarkers := 16
  res1@xyMarkerSizes := 0.001
  plot11 = gsn_csm_xy(wks,ndtooned(u10_2K((/1,2,3,4,5/),:)),ndtooned(wspd_2K((/1,2,3,4,5/),:)),res1)
  res1@xyMarkerColors := colors(1)
  res1@xyMarkers := 16
  res1@xyMarkerSizes := 0.003
  plot12 = gsn_csm_xy(wks,u10_2K(0,:),wspd_2K(0,:),res1)
  rc = regline(u10_2K(0,:),wspd_2K(0,:))
  xr = fspan(-1.5,5.5,10)
  yr = rc*xr+rc@yintercept
  resl = True
  resl@gsLineThicknessF = 1.5
  resl@gsLineDashPattern = 2
  line1 = gsn_add_polyline(wks,plot12,xr,yr,resl)
  restx = True
  restx@txFontHeightF = 0.009
  rc1 = gsn_add_text(wks,plot12,"~F33~b~F~ = "+sprintf("%3.2f",rc)+" ~F34~1~F~ "+sprintf("%3.2f",rc@rstd*cdft_t(0.975,rc@nptxy-2)),3.5,3.0,restx)
  overlay(plot11,plot12)
  draw(plot11)
  res1@tiYAxisString = ""
  res1@vpXF = x1+1*(1+a)*dx
  res1@xyMarkerColors := colors(3)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.001
  plot13 = gsn_csm_xy(wks,ndtooned(u10_2K((/0,1,2,4,5/),:)),ndtooned(wspd_2K((/0,1,2,4,5/),:)),res1)
  res1@xyMarkerColors := colors(1)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.003
  plot14 = gsn_csm_xy(wks,u10_2K(3,:),wspd_2K(3,:),res1)
  rc = regline(u10_2K(3,:),wspd_2K(3,:))
  yr = rc*xr+rc@yintercept
  line2 = gsn_add_polyline(wks,plot14,xr,yr,resl)
  rc2 = gsn_add_text(wks,plot14,"~F33~b~F~ = "+sprintf("%3.2f",rc)+" ~F34~1~F~ "+sprintf("%3.2f",rc@rstd*cdft_t(0.975,rc@nptxy-2)),3.5,3.0,restx)
  overlay(plot13,plot14)
  draw(plot13)
  res1@vpXF = x1+2*(1+a)*dx
  res1@xyMarkerColors := colors(4)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.001
  plot15 = gsn_csm_xy(wks,ndtooned(u10_4K((/1,2,3,4,5/),:)),ndtooned(wspd_4K((/1,2,3,4,5/),:)),res1)
  res1@xyMarkerColors := colors(2)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.003
  plot16 = gsn_csm_xy(wks,u10_4K(0,:),wspd_4K(0,:),res1)
  rc = regline(u10_4K(0,:),wspd_4K(0,:))
  yr = rc*xr+rc@yintercept
  line3 = gsn_add_polyline(wks,plot16,xr,yr,resl)
  rc3 = gsn_add_text(wks,plot16,"~F33~b~F~ = "+sprintf("%3.2f",rc)+" ~F34~1~F~ "+sprintf("%3.2f",rc@rstd*cdft_t(0.975,rc@nptxy-2)),3.5,3.0,restx)
  overlay(plot15,plot16)
  draw(plot15)
  res1@vpXF = x1+3*(1+a)*dx
  res1@xyMarkerColors := colors(4)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.001
  plot17 = gsn_csm_xy(wks,ndtooned(u10_4K((/0,1,2,4,5/),:)),ndtooned(wspd_4K((/0,1,2,4,5/),:)),res1)
  res1@xyMarkerColors := colors(2)
  res1@xyMarkers := 16
  res1@xyMarkerSizes = 0.003
  plot18 = gsn_csm_xy(wks,u10_4K(3,:),wspd_4K(3,:),res1)
  rc = regline(u10_4K(3,:),wspd_4K(3,:))
  yr = rc*xr+rc@yintercept
  line4 = gsn_add_polyline(wks,plot18,xr,yr,resl)
  rc4 = gsn_add_text(wks,plot18,"~F33~b~F~ = "+sprintf("%3.2f",rc)+" ~F34~1~F~ "+sprintf("%3.2f",rc@rstd*cdft_t(0.975,rc@nptxy-2)),3.5,3.0,restx)
  overlay(plot17,plot18)
  draw(plot17)
  
;;; Panels (e-g)
  res1@vpYF = 0.59
  res1@tiXAxisString = "10-m zonal wind clim. [m/s]"
  res1@tiYAxisString = "AUMR variance [(mm/day)~S~2~N~]"
  res1@gsnLeftString = ""
  res1@vpXF = x1+0*(1+1.4*a)*dx
  res1@tmYLPrecision = 2
  res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
  res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
  res1@xyMarkerSizes := 0.009+(res1@xyMarkers*0)
  res1@trXMinF = 0.0
  res1@trXMaxF = 3.0
  res1@tmXBTickSpacingF = 0.5
  res1@trYMinF = 1.2
  res1@trYMaxF = 2.7
  res1@tmYLTickSpacingF = 0.3
  res1@gsnRightStringOrthogonalPosF = -0.94
  res1@gsnRightStringParallelPosF = 0.96
  res1@gsnRightStringFontHeightF = 0.010
  temp = dSST_cor(u10_mean(3:14),var(3:14))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
  plot21 = gsn_csm_xy(wks,diag(u10_mean(3:14)),diag(var(3:14)),res1)
  temp = dSST_cor(sst_mean(3:14),var(3:14))
  res1@xyMarkerColors := colors((/0,1,2/))
  res1@xyMarkers := (/16,16,16/)
  res1@xyMarkerSizes := 0.012+(res1@xyMarkers*0)
  res1@gsnRightString = ""
  plot22 = gsn_csm_xy(wks,diag(u10_mean(0:2)),diag(var(0:2)),res1)
  overlay(plot21,plot22)
  draw(plot21)
  res1@tiXAxisString = "Ensemble mean 10-m zonal wind [m/s]"
  res1@tiYAxisString = "Inter-member AUMR variance [(mm/day)~S~2~N~]"
  res1@vpXF = x1+1*(1+1.4*a)*dx
  res1@xyMarkerColors := colors(0)
  res1@xyMarkers := (/16/)
  res1@xyMarkerSizes := 0.005
  res1@trXMinF = 0.0
  res1@trXMaxF = 4.0
  res1@tmXBTickSpacingF = 1.0
  res1@trYMinF = 0.6
  res1@trYMaxF = 2.7
  res1@tmYLTickSpacingF = 0.3
  res1@gsnRightString = "r = "+sprintf("%3.2f",escorc(u10_hist_ens,var_hist_int))+"*"
  plot23 = gsn_csm_xy(wks,u10_hist_ens,var_hist_int,res1)
  draw(plot23)
  res1@tiYAxisString = "AUMR variance change [(mm/day)~S~2~N~]"
  res1@tiXAxisString = "10-m zonal wind clim. change [m/s]"
  res1@vpXF = x1+2*(1+1.4*a)*dx
  res1@trXMinF = -1.5
  res1@trXMaxF =  1.5
  res1@tmXBTickSpacingF = 0.5
  res1@trYMinF = -3.0
  res1@trYMaxF =  3.0
  res1@tmYLTickSpacingF = 1.0
  res1@xyMarkerColors := m_colors
  res1@xyMarkerSizes := m_sizes
  res1@xyMarkers := m_markers_2K
  temp = dSST_cor_c6(u10_c6_245(M2),u10_c6_585(M2),var_c6_245(M2),var_c6_585(M2))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+dimsizes(M2)+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot24 = gsn_csm_xy(wks,diag(u10_c6_245),diag(var_c6_245),res1)
  res1@xyMarkers := m_markers_4K
  res1@gsnRightString = ""
  plot25 = gsn_csm_xy(wks,diag(u10_c6_585),diag(var_c6_585),res1)
  overlay(plot24,plot25)
  draw(plot24)

;;; Panels (h-i)
  res1@vpYF = 0.26
  d = "~F34~6~F~"
  res1@tiXAxisString = d+"W/"+d+"u~B~s~N~ [(m/s)/(m/s)]"
  res1@tiYAxisString = "AUMR variance [(mm/day)~S~2~N~]"
  res1@gsnLeftString = ""
  res1@vpXF = x1+0*(1+1.4*a)*dx
  res1@tmYLPrecision = 2
  res1@tmXBPrecision = 1
  res1@trXMinF = 0.4
  res1@trXMaxF = 0.8
  res1@tmXBTickSpacingF = 0.1
  res1@trYMinF = 1.2
  res1@trYMaxF = 2.7
  res1@tmYLTickSpacingF = 0.3
  res1@gsnRightStringOrthogonalPosF = -0.94
  res1@gsnRightStringParallelPosF = 0.96
  res1@gsnRightStringFontHeightF = 0.010
  res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
  res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
  res1@xyMarkerSizes := 0.009+(res1@xyMarkers*0)
  temp = dSST_cor(wspd_reg(3:14),var(3:14))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
  plot31 = gsn_csm_xy(wks,diag(wspd_reg(3:14)),diag(var(3:14)),res1)
  res1@xyMarkerColors := colors((/0,1,2/))
  res1@xyMarkers := (/16,16,16/)
  res1@xyMarkerSizes := 0.012+(res1@xyMarkers*0)
  res1@gsnRightString = ""
  plot32 = gsn_csm_xy(wks,diag(wspd_reg(0:2)),diag(var(0:2)),res1)
  overlay(plot31,plot32)
  draw(plot31)
  res1@tiXAxisString = "10-m zonal wind clim. [m/s]"
  res1@tiYAxisString = d+"W/"+d+"u~B~s~N~ [(m/s)/(m/s)]"
  res1@vpXF = x1+1*(1+1.4*a)*dx
  res1@tmXBPrecision = 2
  res1@trXMinF = 0.0
  res1@trXMaxF = 3.0
  res1@tmXBTickSpacingF = 0.5
  res1@tmYLPrecision = 1
  res1@trYMinF = 0.4
  res1@trYMaxF = 0.8
  res1@tmYLTickSpacingF = 0.1
  res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
  res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
  res1@xyMarkerSizes := 0.009+(res1@xyMarkers*0)
  temp = dSST_cor(u10_mean(3:14),wspd_reg(3:14))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")+ \
                        " ~C~ "
  plot33 = gsn_csm_xy(wks,diag(u10_mean(3:14)),diag(wspd_reg(3:14)),res1)
  res1@xyMarkerColors := colors((/0,1,2/))
  res1@xyMarkers := (/16,16,16/)
  res1@xyMarkerSizes := 0.012+(res1@xyMarkers*0)
  res1@gsnRightString = ""
  plot34 = gsn_csm_xy(wks,diag(u10_mean(0:2)),diag(wspd_reg(0:2)),res1)
  overlay(plot33,plot34)
  draw(plot33)
  res1@tiXAxisString = "Relative SST index (WP-MC) [K]"
  res1@tiYAxisString = "10-m zonal wind clim. [m/s]"
  res1@vpXF = x1+2*(1+1.4*a)*dx
  res1@trXMinF = -0.5
  res1@trXMaxF = 1.5
  res1@tmXBTickSpacingF = 0.5
  res1@tmYLPrecision = 2
  res1@trYMinF = 0.0
  res1@trYMaxF = 3.0
  res1@tmYLTickSpacingF = 0.5
  res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
  res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
  res1@xyMarkerSizes := 0.009+(res1@xyMarkers*0)
  temp = dSST_cor(sst_mean(3:14),u10_mean(3:14))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
  plot35 = gsn_csm_xy(wks,diag(sst_mean(3:14)),diag(u10_mean(3:14)),res1)
  res1@xyMarkerColors := colors((/0,1,2/))
  res1@xyMarkers := (/16,16,16/)
  res1@xyMarkerSizes := 0.012+(res1@xyMarkers*0)
  res1@gsnRightString = ""
  plot36 = gsn_csm_xy(wks,diag(sst_mean(0:2)),diag(u10_mean(0:2)),res1)
  overlay(plot35,plot36)
  draw(plot35)
  
  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) d4PDF 2K_CC",x1+0*(1+a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(b) d4PDF 2K_MI",x1+1*(1+a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(c) d4PDF 4K_CC",x1+2*(1+a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(d) d4PDF 4K_MI",x1+3*(1+a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(e) d4PDF HIST/2K/4K",x1+0*(1+a)*dx,0.60,restx)
  gsn_text_ndc(wks,"(f) d4PDF HIST (each year)",x1+1*(1+1.4*a)*dx,0.60,restx)
  gsn_text_ndc(wks,"(g) CMIP6 (26 models)",x1+2*(1+1.4*a)*dx,0.60,restx)
  gsn_text_ndc(wks,"(h) d4PDF HIST/2K/4K",x1+0*(1+a)*dx,0.27,restx)
  gsn_text_ndc(wks,"(i) d4PDF HIST/2K/4K",x1+1*(1+1.4*a)*dx,0.27,restx)
  gsn_text_ndc(wks,"(j) d4PDF HIST/2K/4K",x1+2*(1+1.4*a)*dx,0.27,restx)

  restx@txFontHeightF = 0.016
  gsn_text_ndc(wks,"JF-mean 10-m zonal wind & scalar wind speed for individual years",x1,0.96,restx)
  gsn_text_ndc(wks,"Background monsoon westerlies & AUMR variance",x1,0.63,restx)
  gsn_text_ndc(wks,"AUMR variance, windspeed sensitivity, background monsoon westerlies & SST gradient",x1,0.30,restx)

  resm = True
  restx@txFontHeightF = 0.007
  restx@txJust = "CenterLeft"
  do i = 0, M-1
    ym = 0.608 - 0.018*(i%15)
    xm = 0.791 + (i/15)*0.115
    resm@gsMarkerColor = m_colors(i,:)
    resm@gsMarkerSizeF = m_sizes(i)
    resm@gsMarkerIndex = m_markers_2K(i)
    gsn_polymarker_ndc(wks,xm,ym,resm)
    resm@gsMarkerIndex = m_markers_4K(i)
    gsn_polymarker_ndc(wks,xm+0.011,ym,resm)
    gsn_text_ndc(wks,m_name(i),xm+0.020,ym,restx)
  end do
  restx@txFontHeightF = 0.007
  restx@txJust = "BottomCenter"
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.7963,0.62,restx)
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.9113,0.62,restx)

  restx@txFontHeightF = 0.008
  restx@txJust = "CenterLeft"
  mtext = (/"All","CC","GF","HA","MI","MP","MR"/)
  do i = 0, 13
    resm@gsMarkerIndex := where(i%7.eq.0,16,i%7+3)
    resm@gsMarkerColor := colors((i/7)+1)
    resm@gsMarkerSizeF := where(i%7.eq.0,0.012,0.009)
    xm = 0.835+where(i%7.le.3,i%7,i%7-3)*0.04
    ym = 0.23-where(i%7.le.3,0,1)*0.02-(i/7)*0.05
    gsn_polymarker_ndc(wks,xm,ym,resm)
    gsn_text_ndc(wks,mtext(i%7),xm+0.012,ym,restx)
    if (i%7.eq.0) then
      gsn_text_ndc(wks,where(i/7.eq.0,"2K :","4K :"),xm-0.034,ym,restx)
    end if
  end do
  xm = 0.835
  ym = 0.26
  gsn_text_ndc(wks,"HIST :",xm-0.045,ym,restx)
  resm@gsMarkerIndex = 16
  resm@gsMarkerColor = colors(0)
  resm@gsMarkerSizeF = 0.012
  gsn_polymarker_ndc(wks,xm,ym,resm)
  
;  drawNDCGrid(wks)
  frame(wks)
  
end