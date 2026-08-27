;;; File: fig3.ncl
;;; Description: create Figure 3
;;;
;;; Input: ausm_rain_${experiment_type_name}.nc created by d4PDF/make_ausm_rain.ncl
;;;        lhf_area_${experiment_type_name}.nc created by d4PDF/make_lhf_area.ncl
;;;        ausm_variance.nc created by d4PDF/derive_ausm_variance.ncl
;;;        sst_pc_historical.nc, sst_pc_ssp245.nc and sst_pc_ssp585.nc created by CMIP6/sst_pc.ncl
;;;        ausm_historical.nc, ausm_ssp245.nc and ausm_ssp585.nc created by CMIP6/ausm.ncl
;;;        u10_historical.nc, u10_ssp245.nc and u10_ssp585.nc created by CMIP6/u10.ncl
;;;        gmst_historical.nc, gmst_ssp245.nc and gmst_ssp585.nc created by CMIP6/gmst.ncl
;;;        gwlev.nc created by CMIP6/gwlev.ncl
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

;;; function for calculating dSST-induced correlation (d4PDF)
function dSST_cor(xin[12]:numeric,yin[12]:numeric)
local x, y
begin
  x = xin
  y = yin
  x(:5) = dim_rmvmean(x(:5))
  x(6:) = dim_rmvmean(x(6:))
  y(:5) = dim_rmvmean(y(:5))
  y(6:) = dim_rmvmean(y(6:))
  return(escorc(x,y))
end

;;; function for calculating dSST-induced correlation (CMIP6)
function dSST_cor_c6(x1[*]:numeric,x2[*]:numeric,y1[*]:numeric,y2[*]:numeric)
local x, y
begin
  x = array_append_record(dim_rmvmean(x1),dim_rmvmean(x2),0)
  y = array_append_record(dim_rmvmean(y1),dim_rmvmean(y2),0)
  return(escorc(x,y))
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
  in1 = addfile("../d4PDF/lhf_area_"+expr(0)+".nc","r")
  in2 = addfile("../d4PDF/ausm_rain_"+expr(0)+".nc","r")
  u10_mean(0) = avg(in1->u10)
  u10_hist_ens = dim_avg_n(in1->u10,0)
  var_hist_int = dim_variance_n(in2->ausm,0)
  delete([/in1,in2/])
  in1 = addfiles("../d4PDF/lhf_area_"+expr(2:7)+".nc","r")
  u10_mean(1) = avg(in1[:]->u10)
  ListSetType(in1,"cat")
  in1 = addfiles("../d4PDF/lhf_area_"+expr(8:13)+".nc","r")
  u10_mean(2) = avg(in1[:]->u10)
  ListSetType(in1,"cat")
  delete([/in1/])
  do i = 2, 13
    in1 = addfile("../d4PDF/lhf_area_"+expr(i)+".nc","r")
    u10_mean(i+1) = avg(in1->u10)
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
  cor = escorc(dum,pc1)
  M1 = ind(abs(cor).gt.0.4)
  M2 = ind(abs(cor).lt.0.4)
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
  wks = gsn_open_wks("eps","fig3")
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
  colors = (/"black","orange","firebrick","lightgoldenrod2","pink3"/)
  res1@gsnYRefLine = 0
  res1@gsnXRefLine = 0
  
  res1@vpYF = 0.90
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
  res1@vpYF = 0.57
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
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",dSST_cor(u10_mean(3:14),var(3:14)))
  plot21 = gsn_csm_xy(wks,diag(u10_mean(3:14)),diag(var(3:14)),res1)
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
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",escorc(u10_hist_ens,var_hist_int))
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
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",dSST_cor_c6(u10_c6_245(M2),u10_c6_585(M2),var_c6_245(M2),var_c6_585(M2)))
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot24 = gsn_csm_xy(wks,diag(u10_c6_245),diag(var_c6_245),res1)
  res1@xyMarkers := m_markers_4K
  res1@gsnRightString = ""
  plot25 = gsn_csm_xy(wks,diag(u10_c6_585),diag(var_c6_585),res1)
  overlay(plot24,plot25)
  draw(plot24)

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) 2K_CC",x1+0*(1+a)*dx,0.91,restx)
  gsn_text_ndc(wks,"(b) 2K_MI",x1+1*(1+a)*dx,0.91,restx)
  gsn_text_ndc(wks,"(c) 4K_CC",x1+2*(1+a)*dx,0.91,restx)
  gsn_text_ndc(wks,"(d) 4K_MI",x1+3*(1+a)*dx,0.91,restx)
  gsn_text_ndc(wks,"(e) d4PDF HIST/2K/4K",x1+0*(1+a)*dx,0.58,restx)
  gsn_text_ndc(wks,"(f) d4PDF HIST (each year)",x1+1*(1+1.4*a)*dx,0.58,restx)
  gsn_text_ndc(wks,"(g) CMIP6 (26 models)",x1+2*(1+1.4*a)*dx,0.58,restx)

  restx@txFontHeightF = 0.016
  gsn_text_ndc(wks,"Relationship b/w JF-mean 10-m zonal wind & scalar wind speed",x1,0.94,restx)
  gsn_text_ndc(wks,"Relationship b/w background monsoon westerlies & AUMR variability",x1,0.61,restx)

  resm = True
  restx@txFontHeightF = 0.007
  restx@txJust = "CenterLeft"
  do i = 0, M-1
    ym = 0.558 - 0.018*(i%15)
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
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.7963,0.57,restx)
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.9113,0.57,restx)

  restx@txFontHeightF = 0.008
  restx@txJust = "CenterLeft"
  mtext = (/"All","CC","GF","HA","MI","MP","MR"/)
  do i = 0, 13
    resm@gsMarkerIndex := where(i%7.eq.0,16,i%7+3)
    resm@gsMarkerColor := colors((i/7)+1)
    resm@gsMarkerSizeF := where(i%7.eq.0,0.012,0.009)
    xm = 0.17+(i%7)*0.04
    ym = 0.30-(i/7)*0.02
    gsn_polymarker_ndc(wks,xm,ym,resm)
    gsn_text_ndc(wks,mtext(i%7),xm+0.012,ym,restx)
    if (i%7.eq.0) then
      gsn_text_ndc(wks,where(i/7.eq.0,"2K :","4K :"),xm-0.036,ym,restx)
    end if
  end do
  xm = 0.1
  gsn_text_ndc(wks,"HIST :",xm-0.05,ym+0.02,restx)
  resm@gsMarkerIndex = 16
  resm@gsMarkerColor = colors(0)
  resm@gsMarkerSizeF = 0.012
  gsn_polymarker_ndc(wks,xm,ym+0.02,resm)
  
;  drawNDCGrid(wks)
  frame(wks)
  
end