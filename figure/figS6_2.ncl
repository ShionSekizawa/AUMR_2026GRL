;;; File: figS6_2.ncl
;;; Description: create Figure S6 (lower part)
;;;
;;; Input: sst_pc_historical_mulens.nc and sst_pc_ssp585_mulens.nc created by CMIP6/mulens/sst_pc_mulens.ncl
;;;        ausm_historical_mulens.nc and ausm_ssp585_mulens.nc created by CMIP6/mulens/ausm_mulens.ncl
;;;        u10_historical_mulens.nc and u10_ssp585_mulens.nc created by CMIP6/mulens/u10_mulens.ncl
;;;        sst_historical_mulens.nc and sst_ssp585_mulens.nc created by CMIP6/mulens/sst_mulens.ncl
;;;        CMIP6/mulens/u10_cov_mulens.txt created by CMIP6/mulens/cov_mulens.ncl
;;;        CMIP6/gwlev_mulens.nc created by CMIP6/mulens/gwlev.ncl
;;;        CMIP6/mulens/ModelMemList_mulens.txt
;;; Output: figS6_2.eps
;;;
;;; Running:
;;; ncl figS6_2.ncl
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
  if (dimsizes(dimsizes(sst)).eq.4) then
    mc  = wgt_areaave(sst(:,:,{-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst(:,:,{-ss:ss},{w3:e3}),wgt,1.,0)
  end if
  if (dimsizes(dimsizes(sst)).eq.3) then
    mc  = wgt_areaave(sst(:,{-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst(:,{-ss:ss},{w3:e3}),wgt,1.,0)
  end if
  ret = wp-mc
  return(ret)
end

begin

  ;;; CMIP6 data
  M = 13
  in = addfile("../CMIP6/mulens/sst_pc_historical_mulens.nc","r")
  pc1  = reshape(in->pc(:,:,0,:),(/13*5,60/))
  pc1m = reshape(in->pc(:,:,0,:),(/13,5*60/))
  in = addfile("../CMIP6/mulens/ausm_historical_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  dum  = reshape(dim_rmvmean(in->ausm(:,:,id)),(/13*5,60/))
  dumm = reshape(dim_rmvmean_n(in->ausm(:,:,id),(/1,2/)),(/13,5*60/))
  dum  = dtrend(dum,False)
  dumm = dtrend_msg(ndtooned(conform_dims((/5,60/),yyyy(id),1)),dumm,False,False)
  do i = 0, 13*5-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  do i = 0, 13-1
    dumm(i,:) = dtrend_msg(pc1m(i,:),dumm(i,:),False,False)
  end do
  cov = asciiread("../CMIP6/mulens/u10_cov_mulens.txt",13,"float")
  M1 = ind(cov.lt.0.12)
  M2 = ind(cov.ge.0.12)
  var_c6_hist  = dim_variance(dum)
  var_c6m_hist = dim_variance(dumm)
  in = addfile("../CMIP6/mulens/u10_historical_mulens.nc","r")
  u10_c6_hist  = ndtooned(dim_avg(in->u10(:,:,id)))
  u10_c6m_hist = dim_avg_n(in->u10(:,:,id),(/1,2/))
  in = addfile("../CMIP6/mulens/sst_historical_mulens.nc","r")
  sst_c6_hist = ndtooned(sst_index(dim_avg_n_Wrap(in->sst(:,:,id,:,:),2)))
  sst_c6m_hist = sst_index(dim_avg_n_Wrap(in->sst(:,:,id,:,:),(/1,2/)))
  delete([/yyyy,id,dum,dumm,pc1,pc1m/])

  in = addfile("../CMIP6/mulens/sst_pc_ssp585_mulens.nc","r")
  pc1  = reshape(in->pc(:,:,0,:),(/13*5,60/))
  pc1m = reshape(in->pc(:,:,0,:),(/13,5*60/))
  in = addfile("../CMIP6/mulens/ausm_ssp585_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  dum  = reshape(dim_rmvmean(in->ausm(:,:,id)),(/13*5,60/))
  dumm = reshape(dim_rmvmean_n(in->ausm(:,:,id),(/1,2/)),(/13,5*60/))
  dum  = dtrend(dum,False)
  dumm = dtrend_msg(ndtooned(conform_dims((/5,60/),yyyy(id),1)),dumm,False,False)
  do i = 0, 13*5-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  do i = 0, 13-1
    dumm(i,:) = dtrend_msg(pc1m(i,:),dumm(i,:),False,False)
  end do
  var_c6_585  = dim_variance(dum)
  var_c6m_585 = dim_variance(dumm)
  in = addfile("../CMIP6/mulens/u10_ssp585_mulens.nc","r")
  u10_c6_585  = ndtooned(dim_avg(in->u10(:,:,id)))
  u10_c6m_585 = dim_avg_n(in->u10(:,:,id),(/1,2/))
  in = addfile("../CMIP6/mulens/sst_ssp585_mulens.nc","r")
  sst_c6_585 = ndtooned(sst_index(dim_avg_n_Wrap(in->sst(:,:,id,:,:),2)))
  sst_c6m_585 = sst_index(dim_avg_n_Wrap(in->sst(:,:,id,:,:),(/1,2/)))
  delete([/yyyy,id,dum,dumm,pc1,pc1m/])

  ;;; scaled by global warming level
  in = addfile("../CMIP6/mulens/gwlev_mulens.nc","r")
  gwlev_hist = ndtooned(in->gwlev_historical)
  gwlev_585  = ndtooned(in->gwlev_ssp585)
  var_c6_585 = (var_c6_585-var_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  u10_c6_585 = (u10_c6_585-u10_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  sst_c6_585 = (sst_c6_585-sst_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  delete([/gwlev_hist,gwlev_585/])
  gwlev_hist = dim_avg(in->gwlev_historical)
  gwlev_585  = dim_avg(in->gwlev_ssp585)
  var_c6m_585 = (var_c6m_585-var_c6m_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  u10_c6m_585 = (u10_c6m_585-u10_c6m_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  sst_c6m_585 = (sst_c6m_585-sst_c6m_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  delete([/gwlev_hist,gwlev_585/])
 
  m_tmp  = asciiread("../CMIP6/ModelMemList.txt",26,"string")
  m_name = new(26,string)
  do i = 0, 26-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do
  mm_tmp = asciiread("../CMIP6/mulens/ModelMemList_mulens.txt",13*5,"string")
  m_name_index = new(13,integer)
  do i = 0, 13-1
    tmp = str_split(mm_tmp(i*5)," ")
    m_name_index(i) = ind(m_name.eq.tmp(1))
    delete(tmp)
  end do
  m_name := m_name(m_name_index)
  print(m_name)

  cmp = read_colormap_file("NCV_banded")
  NC = dimsizes(cmp(:,0))
  m_colors = cmp((ispan(1,2*26-1,2)*NC)/(2*26)-1,:)
  mm_colors = m_colors(m_name_index,:)
  m_colors := reshape(conform_dims((/13,5,4/),m_colors(m_name_index,:),(/0,2/)),(/13*5,4/))
  m_markers  = new(13*5,integer)
  mm_markers = new(13,integer)
  m_markers  = 16
  mm_markers(M1) = 3
  mm_markers(M2) = 6
  m_sizes  = new(13*5,float)
  mm_sizes = new(13,float)
  m_sizes  = 0.0012
  mm_sizes(M1) = 0.006
  mm_sizes(M2) = 0.007

;;; open workstation
  wks = gsn_open_wks("eps","figS6_2")
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
  res1@xyMarkers := m_markers
  temp = dSST_cor_c6(u10_c6m_585(M2),u10_c6m_585(M2),var_c6m_585(M2),var_c6m_585(M2))
  res1@gsnRightStringOrthogonalPosF = -0.94
  res1@gsnRightStringParallelPosF = 0.96
  res1@gsnRightStringFontHeightF = 0.010
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+dimsizes(M2)+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot11 = gsn_csm_xy(wks,diag(u10_c6_585),diag(var_c6_585),res1)
  res1@xyMarkerColors := mm_colors
  res1@xyMarkerSizes := mm_sizes
  res1@xyMarkers := mm_markers
  res1@gsnRightString = ""
  plot12 = gsn_csm_xy(wks,diag(u10_c6m_585),diag(var_c6m_585),res1)
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
  res1@xyMarkers := m_markers
  temp = dSST_cor_c6(sst_c6m_585(M2),sst_c6m_585(M2),var_c6m_585(M2),var_c6m_585(M2))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+dimsizes(M2)+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot13 = gsn_csm_xy(wks,diag(sst_c6_585),diag(var_c6_585),res1)
  res1@xyMarkerColors := mm_colors
  res1@xyMarkerSizes := mm_sizes
  res1@xyMarkers := mm_markers
  res1@gsnRightString = ""
  plot14 = gsn_csm_xy(wks,diag(sst_c6m_585),diag(var_c6m_585),res1)
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
  res1@xyMarkers := m_markers
  temp = dSST_cor_c6(sst_c6m_585(M2),sst_c6m_585(M2),u10_c6m_585(M2),u10_c6m_585(M2))
  res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")+" ~C~ ~S~             ("+M+" models)"
  res1@gsnRightStringOrthogonalPosF = -0.94
  plot15 = gsn_csm_xy(wks,diag(sst_c6_585),diag(u10_c6_585),res1)
  res1@xyMarkerColors := mm_colors
  res1@xyMarkerSizes := mm_sizes
  res1@xyMarkers := mm_markers
  res1@gsnRightString = ""
  plot16 = gsn_csm_xy(wks,diag(sst_c6m_585),diag(u10_c6m_585),res1)
  overlay(plot15,plot16)
  draw(plot15)

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(d) CMIP6 (13 models x 5 mem)",x1+0*(1+1.4*a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(e) CMIP6 (13 models x 5 mem)",x1+1*(1+1.4*a)*dx,0.93,restx)
  gsn_text_ndc(wks,"(f) CMIP6 (13 models x 5 mem)",x1+2*(1+1.4*a)*dx,0.93,restx)

  restx@txFontHeightF = 0.016
;  gsn_text_ndc(wks,"AUMR variance, background monsoon westerlies & relative SST gradient changes",x1,0.96,restx)

  resm = True
  restx@txFontHeightF = 0.007
  restx@txJust = "CenterLeft"
  do i = 0, 13-1
    ym = 0.910 - 0.018*(i%15)
    xm = 0.821 + (i/15)*0.115
    resm@gsMarkerColor = mm_colors(i,:)
    resm@gsMarkerSizeF = mm_sizes(i)
    resm@gsMarkerIndex = mm_markers(i)
    gsn_polymarker_ndc(wks,xm,ym,resm)
    resm@gsMarkerColor = m_colors(i*5,:)
    resm@gsMarkerSizeF = m_sizes(i*5)
    resm@gsMarkerIndex = m_markers(i*5)
    gsn_polymarker_ndc(wks,xm+0.011,ym,resm)
    gsn_text_ndc(wks,m_name(i),xm+0.020,ym,restx)
  end do
  restx@txFontHeightF = 0.007
  restx@txJust = "BottomCenter"
;  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.7963,0.92,restx)
;  gsn_text_ndc(wks,"2K~S~ ~N~4K",0.9113,0.92,restx)

  ;  drawNDCGrid(wks)
  frame(wks)
end