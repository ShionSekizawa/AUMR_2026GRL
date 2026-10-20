;;; File: figS6_1.ncl
;;; Description: create Figure S6 (upper part)
;;;
;;; Input: sst_pc_historical.nc, sst_pc_ssp245.nc and sst_pc_ssp585.nc created by CMIP6/sst_pc.ncl
;;;        ausm_historical.nc, ausm_ssp245.nc and ausm_ssp585.nc created by CMIP6/ausm.ncl
;;;        u10_historical.nc, u10_ssp245.nc and u10_ssp585.nc created by CMIP6/u10.ncl
;;;        sst_historical.nc, sst_ssp245.nc and sst_ssp585.nc created by CMIP6/sst.ncl
;;;        u10_cov.txt created by CMIP6/cov.ncl
;;;        gwlev.nc created by CMIP6/gwlev.ncl
;;;        CMIP6/ModelMemList.txt
;;; Output: figS6_1.eps
;;;
;;; Running:
;;; ncl figS6_1.ncl
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

  ;;; CMIP6 data
  in = addfile("../CMIP6/sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  M = dimsizes(pc1(:,0))
  in = addfile("../CMIP6/ausm_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  dum = dim_rmvmean(in->ausm(:,id))
  dum = dtrend(dum,False)
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  cov = asciiread("../CMIP6/u10_cov.txt",26,"float")
  M1 = ind(cov.lt.0.12)
  M2 = ind(cov.ge.0.12)
  print(dimsizes(M2))
  var_c6_hist = dim_variance(dum)
  in = addfile("../CMIP6/u10_historical.nc","r")
  u10_c6_hist = dim_avg(in->u10(:,id))
  in = addfile("../CMIP6/sst_historical.nc","r")
  sst_c6_hist = sst_index(dim_avg_n_Wrap(in->sst(:,id,:,:),1))
  delete([/yyyy,id,dum,pc1/])

  in = addfile("../CMIP6/sst_pc_ssp245.nc","r")
  pc1 = in->pc(:,0,:)
  in = addfile("../CMIP6/ausm_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  dum = in->ausm(:,id)
  dum = dim_rmvmean(in->ausm(:,id))
  dum = dtrend(dum,False)
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  var_c6_245 = dim_variance(dum)
  in = addfile("../CMIP6/u10_ssp245.nc","r")
  u10_c6_245 = dim_avg(in->u10(:,id))
  in = addfile("../CMIP6/sst_ssp245.nc","r")
  sst_c6_245 = sst_index(dim_avg_n_Wrap(in->sst(:,id,:,:),1))
  delete([/yyyy,id,dum,pc1/])

  in = addfile("../CMIP6/sst_pc_ssp585.nc","r")
  pc1 = in->pc(:,0,:)
  in = addfile("../CMIP6/ausm_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  dum = in->ausm(:,id)
  dum = dim_rmvmean(in->ausm(:,id))
  dum = dtrend(dum,False)
  do i = 0, M-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  var_c6_585 = dim_variance(dum)
  in = addfile("../CMIP6/u10_ssp585.nc","r")
  u10_c6_585 = dim_avg(in->u10(:,id))
  in = addfile("../CMIP6/sst_ssp585.nc","r")
  sst_c6_585 = sst_index(dim_avg_n_Wrap(in->sst(:,id,:,:),1))  
  delete([/yyyy,id,dum,pc1/])

  ;;; scaled by global warming level
  in = addfile("../CMIP6/gwlev.nc","r")
  gwlev_hist = in->gwlev_historical
  gwlev_245  = in->gwlev_ssp245
  gwlev_585  = in->gwlev_ssp585
  var_c6_245 = (var_c6_245-var_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  var_c6_585 = (var_c6_585-var_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  u10_c6_245 = (u10_c6_245-u10_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  u10_c6_585 = (u10_c6_585-u10_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  sst_c6_245 = (sst_c6_245-sst_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  sst_c6_585 = (sst_c6_585-sst_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
 
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
  wks = gsn_open_wks("eps","figS6_1")
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
  res1@gsnYRefLine = 0
  res1@gsnXRefLine = 0

  res1@vpYF = 0.92
  x1 = 0.045
  x2 = 0.99
  a  = 0.20
  dx = (x2-x1)/(4+3*a)
  res1@vpWidthF  = dx
  res1@vpHeightF = dx
  
  res1@tiYAxisString = "AUMR variance change [(mm/day)~S~2~N~]"
  res1@tiXAxisString = "10-m zonal wind clim. change [m/s]"
  res1@vpXF = x1+0*(1+1.4*a)*dx
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
  res1@gsnRightStringOrthogonalPosF = -0.94
  res1@gsnRightStringParallelPosF = 0.96
  res1@gsnRightStringFontHeightF = 0.010
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+dimsizes(M2)+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot11 = gsn_csm_xy(wks,diag(u10_c6_245),diag(var_c6_245),res1)
  res1@xyMarkers := m_markers_4K
  res1@gsnRightString = ""
  plot12 = gsn_csm_xy(wks,diag(u10_c6_585),diag(var_c6_585),res1)
  overlay(plot11,plot12)
  draw(plot11)

  res1@tiYAxisString = "AUMR variance change [(mm/day)~S~2~N~]"
  res1@tiXAxisString = "Relative SST change (WP-MC) [K]"
  res1@vpXF = x1+1*(1+1.4*a)*dx
  res1@trXMinF = -0.9
  res1@trXMaxF =  0.6
  res1@tmXBTickSpacingF = 0.3
  res1@tmXBPrecision = 1
  res1@trYMinF = -3.0
  res1@trYMaxF =  3.0
  res1@tmYLTickSpacingF = 1.0
  res1@xyMarkerColors := m_colors
  res1@xyMarkerSizes := m_sizes
  res1@xyMarkers := m_markers_2K
  temp = dSST_cor_c6(sst_c6_245(M2),sst_c6_585(M2),var_c6_245(M2),var_c6_585(M2))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+dimsizes(M2)+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot13 = gsn_csm_xy(wks,diag(sst_c6_245),diag(var_c6_245),res1)
  res1@xyMarkers := m_markers_4K
  res1@gsnRightString = ""
  plot14 = gsn_csm_xy(wks,diag(sst_c6_585),diag(var_c6_585),res1)
  overlay(plot13,plot14)
  draw(plot13)

  res1@tiYAxisString = "10-m zonal wind clim. change [m/s]"
  res1@tiXAxisString = "Relative SST change (WP-MC) [K]"
  res1@vpXF = x1+2*(1+1.4*a)*dx
  res1@trXMinF = -0.9
  res1@trXMaxF =  0.6
  res1@tmXBTickSpacingF = 0.3
  res1@trYMinF = -1.5
  res1@trYMaxF =  1.5
  res1@tmYLTickSpacingF = 0.5
  res1@xyMarkerColors := m_colors
  res1@xyMarkerSizes := m_sizes
  res1@xyMarkers := m_markers_2K
  temp = dSST_cor_c6(sst_c6_245(:),sst_c6_585(:),u10_c6_245(:),u10_c6_585(:))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+M+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot15 = gsn_csm_xy(wks,diag(sst_c6_245),diag(u10_c6_245),res1)
  res1@xyMarkers := m_markers_4K
  res1@gsnRightString = ""
  plot16 = gsn_csm_xy(wks,diag(sst_c6_585),diag(u10_c6_585),res1)
  overlay(plot15,plot16)
  draw(plot15)

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) CMIP6 (26 models)",x1+0*(1+1.4*a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(b) CMIP6 (26 models)",x1+1*(1+1.4*a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(c) CMIP6 (26 models)",x1+2*(1+1.4*a)*dx,0.93,restx)

  restx@txFontHeightF = 0.016
  gsn_text_ndc(wks,"AUMR variance, background monsoon westerlies & relative SST changes",x1,0.96,restx)

  resm = True
  restx@txFontHeightF = 0.007
  restx@txJust = "CenterLeft"
  do i = 0, M-1
    ym = 0.910 - 0.018*(i%15)
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
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.7963,0.92,restx)
  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.9113,0.92,restx)

  ;  drawNDCGrid(wks)
  frame(wks)
end