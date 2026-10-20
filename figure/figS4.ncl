;;; File: figS4.ncl
;;; Description: create Figure S4
;;;
;;; Input: lhf_area_${experiment_type_name}.nc created by d4PDF/make_lhf_area.ncl
;;;        sst_clim_${experiment_type_name}.nc created by d4PDF/make_sst_clim.ncl
;;;        sst_pc_historical.nc, sst_pc_ssp245.nc and sst_pc_ssp585.nc created by CMIP6/sst_pc.ncl
;;;        sst_historical.nc, sst_ssp245.nc, sst_ssp585.nc created by CMIP6/sst.ncl
;;;        u10_historical.nc, u10_ssp245.nc, u10_ssp585.nc created by CMIP6/u10.ncl
;;;        gwlev.nc created by CMIP6/gwlev.ncl
;;;        CMIP6/ModelMemList.txt
;;; Output: figS4.eps
;;;
;;; Running:
;;; ncl figS4.ncl
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
local x, y, cor, reg, ret, N, t, r
begin
  x = array_append_record(dim_rmvmean(x1),dim_rmvmean(x2),0)
  y = array_append_record(dim_rmvmean(y1),dim_rmvmean(y2),0)
  N = dimsizes(x1)
  t = cdft_t(0.975,N)
  r = t/sqrt(N+t^2)
  cor  = escorc(x,y)
  reg  = regCoef(x,y)
  ret  = (/cor,reg/)
  ret@cor2 = escorc(x1,y1)
  ret@cor4 = escorc(x2,y2)
  ret@reg2 = regCoef(x1,y1)
  ret@reg4 = regCoef(x2,y2)
  ret@stat = abs(cor).gt.r
  print("dSST_cor_c6 stats")
  print("   245: "+ret@cor2+", 585: "+ret@cor4+", 245/585: "+cor+", sig_cor="+r)
  return(ret)
end

function sst_index(sst,lat)
local trp, mc, wp, ep, io, wgt, ss, w1, w2, w3, w4, e1, e2, e3, e4, ret
begin
  w1 =  60
  e1 = 105
  w2 = 105
  e2 = 135
  w3 = 135
  e3 = 180
  w4 = 225
  e4 = 270
  ss  = 10.
  wgt = cos(lat({-ss:ss})*atan(1.)/45)
  if (dimsizes(dimsizes(sst)).eq.3) then
    trp = wgt_areaave(sst(:,{-ss:ss},:),wgt,1.,0)
    io  = wgt_areaave(sst(:,{-ss:ss},{w1:e1}),wgt,1.,0)
    mc  = wgt_areaave(sst(:,{-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst(:,{-ss:ss},{w3:e3}),wgt,1.,0)
    ep  = wgt_areaave(sst(:,{-ss:ss},{w4:e4}),wgt,1.,0)
  end if
  if (dimsizes(dimsizes(sst)).eq.2) then
    trp = wgt_areaave(sst({-ss:ss},:),wgt,1.,0)
    io  = wgt_areaave(sst({-ss:ss},{w1:e1}),wgt,1.,0)
    mc  = wgt_areaave(sst({-ss:ss},{w2:e2}),wgt,1.,0)
    wp  = wgt_areaave(sst({-ss:ss},{w3:e3}),wgt,1.,0)
    ep  = wgt_areaave(sst({-ss:ss},{w4:e4}),wgt,1.,0)
  end if
  ret = (/mc-trp,io-mc,wp-mc,wp-ep/)
  return(ret)
end

;;; procedure for display SST_index correlation matrix
procedure cor_mat(x)
local str, i, j, N, cor, cor_inv
begin
  N = dimsizes(x(:,0))
  print("SST_index correlation matrix")
  do i = 0, N-1
    str = sprintf("%5.4f",escorc(x(i,:),x(0,:)))
    do j = 1, N-1
      str = str + "   " +sprintf("%5.4f",escorc(x(i,:),x(j,:)))
    end do
    print("  "+str)
  end do
  cor = new((/N-1,N-1/),float)
  do i = 1, N-1
    do j = 1, N-1
      cor(i-1,j-1) = escorc(x(i,:),x(j,:))
    end do
  end do
  cor_inv = inverse_matrix(cor)
  do i = 0, N-2
    print(" vif:"+cor_inv(i,i))
  end do
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
  in1 = addfile("../d4PDF/lhf_area_"+expr(0)+".nc","r")
  in2 = addfile("../d4PDF/sst_clim_"+expr(0)+".nc","r")
  lat = in2->lat
  u10_hist = avg(in1->u10)
  sst_hist = sst_index(in2->sst,lat)
  delete([/in1,in2/])
  in1 = addfiles("../d4PDF/lhf_area_"+expr(2:7)+".nc","r")
  in2 = addfiles("../d4PDF/sst_clim_"+expr(2:7)+".nc","r")
  ListSetType(in1,"join")
  ListSetType(in2,"join")
  u10_2K = dim_avg_n(in1[:]->u10(:,:,:),(/1,2/))
  sst_2K = sst_index(in2[:]->sst(:,:,:),lat)
  in1 = addfiles("../d4PDF/lhf_area_"+expr(8:13)+".nc","r")
  in2 = addfiles("../d4PDF/sst_clim_"+expr(8:13)+".nc","r")
  ListSetType(in1,"join")
  ListSetType(in2,"join")
  u10_4K = dim_avg_n(in1[:]->u10(:,:,:),(/1,2/))
  sst_4K = sst_index(in2[:]->sst(:,:,:),lat)
  delete([/in1,in2,lat/])
  u10_2K4K = array_append_record(u10_2K,u10_4K,0)
  sst_2K4K = transpose(array_append_record(transpose(sst_2K),transpose(sst_4K),0))
  
;;; CMIP6 data (single member)
  in = addfile("../CMIP6/sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  M = dimsizes(pc1(:,0))
  in = addfile("../CMIP6/sst_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  lat = in->lat
  lon = in->lon
  M2 = ispan(0,M-1,1)
  sst_c6_hist = sst_index(dim_avg_n_Wrap(in->sst(M2,id,:,:),1),lat)
  in = addfile("../CMIP6/u10_historical.nc","r")
  u10_c6_hist = dim_avg(in->u10(M2,id))
  delete([/yyyy,id/])

  in = addfile("../CMIP6/sst_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  sst_c6_245 = sst_index(dim_avg_n_Wrap(in->sst(M2,id,:,:),1),lat)
  in = addfile("../CMIP6/u10_ssp245.nc","r")
  u10_c6_245 = dim_avg(in->u10(M2,id))
  delete([/yyyy,id/])

  in = addfile("../CMIP6/sst_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  sst_c6_585 = sst_index(dim_avg_n_Wrap(in->sst(M2,id,:,:),1),lat)
  in = addfile("../CMIP6/u10_ssp585.nc","r")
  u10_c6_585 = dim_avg(in->u10(M2,id))
  delete([/yyyy,id/])

  in = addfile("../CMIP6/gwlev.nc","r")
  gwlev_hist = in->gwlev_historical(M2)
  gwlev_245 = in->gwlev_ssp245(M2)
  gwlev_585 = in->gwlev_ssp585(M2)
  sst_c6_245 = sst_c6_245-sst_c6_hist
  sst_c6_585 = sst_c6_585-sst_c6_hist
  do i = 0, dimsizes(M2)-1
    sst_c6_245(:,i) = sst_c6_245(:,i)*(2.-gwlev_hist(i))/(gwlev_245(i)-gwlev_hist(i))
    sst_c6_585(:,i) = sst_c6_585(:,i)*(4.-gwlev_hist(i))/(gwlev_585(i)-gwlev_hist(i))
  end do
  u10_c6_245 = (u10_c6_245-u10_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  u10_c6_585 = (u10_c6_585-u10_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  delete([/lat,lon,gwlev_hist,gwlev_245,gwlev_585/])
  u10_c6 = array_append_record(u10_c6_245,u10_c6_585,0)
  sst_c6 = transpose(array_append_record(transpose(sst_c6_245),transpose(sst_c6_585),0))

  cor_mat(sst_c6_245)
  cor_mat(sst_c6_585)
  cor_mat(transpose(array_append_record(transpose(dim_rmvmean(sst_c6_245)),transpose(dim_rmvmean(sst_c6_585)),0)))

;;; multi-linear regression
  opt = True
  opt@print_anova = True 
  opt@print_data  = False
  xid = (/1,2,3/)
  print("reg_multlin_stats: CMIP6 245")
  mulin_y = dim_standardize(dim_rmvmean(u10_c6_245),0)
  mulin_x = transpose(dim_standardize(dim_rmvmean(sst_c6_245(xid,:)),0))
  b = reg_multlin_stats(mulin_y,mulin_x,opt)
  print(b)
  print(b@tval/sqrt(b@tval^2+b@MSE_dof))
  print("reg_multlin_stats: CMIP6 585")
  mulin_y = dim_standardize(dim_rmvmean(u10_c6_585),0)
  mulin_x = transpose(dim_standardize(dim_rmvmean(sst_c6_585(xid,:)),0))
  b = reg_multlin_stats(mulin_y,mulin_x,opt)
  print(b)
  print(b@tval/sqrt(b@tval^2+b@MSE_dof))
  print("reg_multlin_stats: CMIP6 245+585")
  mulin_y = dim_standardize(dim_rmvmean(u10_c6_245),0)
  mulin_x = transpose(dim_standardize(dim_rmvmean(sst_c6_245(xid,:)),0))
  mulin_y := array_append_record(mulin_y,dim_standardize(dim_rmvmean(u10_c6_585),0),0)
  mulin_x := array_append_record(mulin_x,transpose(dim_standardize(dim_rmvmean(sst_c6_585(xid,:)),0)),0)
  b = reg_multlin_stats(mulin_y,mulin_x,opt) 
  print(b)
  print(b@tval/sqrt(b@tval^2+b@MSE_dof))

  m_tmp = asciiread("../CMIP6/ModelMemList.txt",M,"string")
  m_name = new(M,string)
  do i = 0, M-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do

  cmp = read_colormap_file("NCV_banded")
  NC = dimsizes(cmp(:,0))
  m_colors  = cmp((ispan(1,2*M-1,2)*NC)/(2*M)-1,:)
  m_colors := m_colors(M2,:)
  m_markers_2K = new(dimsizes(M2),integer)
  m_markers_4K = new(dimsizes(M2),integer)
  m_markers_2K = 4
  m_markers_4K = 6
  m_sizes = new(dimsizes(M2),float)
  m_sizes = 0.007

;;; open workstation
  wks = gsn_open_wks("eps","figS4")
  plot1 = new(4,graphic)
  plot2 = new(4,graphic)
  plot3 = new(4,graphic)
  plot4 = new(4,graphic)
  plot5 = new(4,graphic)
  
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
  
  res1@vpYF = 0.88
  x1 = 0.045
  x2 = 0.99
  a  = 0.20
  dx = (x2-x1)/(4+3*a)
  res1@vpWidthF  = dx
  res1@vpHeightF = dx

  printMinMax(sst_2K4K,0)
  printMinMax(u10_2K4K,0)
  printMinMax(sst_c6,0)
  printMinMax(u10_c6,0)

  XMin = (/0.,-0.4,-0.4,1.2/)
  XMax = (/1.,0.4,1.2,3.0/)
  XInt = (/0.2,0.2,0.4,0.3/)
  XPre = (/2,1,2,2/)
  res1@tiXAxisString = "Relative SST [K]"
  res1@tiYAxisString = "10-m zonal wind clim. [m/s]"
  res1@gsnLeftString = ""
  do i = 0, 3
    res1@vpXF = x1+i*(1+a)*dx
    res1@tmYLPrecision = 2
    res1@tmXBPrecision = XPre(i)
    res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
    res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
    res1@xyMarkerSizes := 0.009+(res1@xyMarkers*0)
    res1@trXMinF = XMin(i)
    res1@trXMaxF = XMax(i)
    res1@tmXBTickSpacingF = XInt(i)
    res1@trYMinF = 0.0
    res1@trYMaxF = 3.0
    res1@tmYLTickSpacingF = 0.5
    res1@gsnRightStringOrthogonalPosF = -0.94
    res1@gsnRightStringParallelPosF = 0.96
    res1@gsnRightStringFontHeightF = 0.010
    temp = dSST_cor(sst_2K4K(i,:),u10_2K4K)
    res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
    plot1(i) = gsn_csm_xy(wks,diag(sst_2K4K(i,:)),diag(u10_2K4K),res1)
    res1@xyMarkerColors := colors((/0,1,2/))
    res1@xyMarkers := (/16,16,16/)
    res1@xyMarkerSizes := 0.012+(res1@xyMarkers*0)
    res1@gsnRightString = ""
    plot2(i) = gsn_csm_xy(wks,diag((/sst_hist(i),avg(sst_2K(i,:)),avg(sst_4K(i,:))/)),diag((/u10_hist,avg(u10_2K),avg(u10_4K)/)),res1)
    overlay(plot1(i),plot2(i))
    draw(plot1(i))
    res1@tiYAxisString = ""
  end do

  XMin = (/-0.4,-0.8,-0.8,-1.2/)
  XMax = (/0.4,0.8,0.8,0.4/)
  XInt = (/0.2,0.4,0.4,0.4/)
  XPre = (/1,1,1,2/)
  res1@vpYF = 0.50
  res1@tiXAxisString = "Relative SST change [K]"
  res1@tiYAxisString = "10-m zonal wind clim. change [m/s]"
  do i = 0, 3
    res1@vpXF = x1+i*(1+a)*dx
    res1@trXMinF = XMin(i)
    res1@trXMaxF = XMax(i)
    res1@tmXBTickSpacingF = XInt(i)
    res1@tmXBPrecision = XPre(i)
    res1@trYMinF = -1.5
    res1@trYMaxF =  1.5
    res1@tmYLTickSpacingF = 0.5
    res1@xyMarkerColors := m_colors
    res1@xyMarkerSizes := m_sizes
    res1@xyMarkers := m_markers_2K
    temp = dSST_cor_c6(sst_c6_245(i,:),sst_c6_585(i,:),u10_c6_245,u10_c6_585)
    res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"*","")
    res1@gsnRightStringOrthogonalPosF = -0.94
    plot3(i) = gsn_csm_xy(wks,diag(sst_c6_245(i,:)),diag(u10_c6_245),res1)
    res1@xyMarkers := m_markers_4K
    res1@gsnRightString = ""
    plot4(i) = gsn_csm_xy(wks,diag(sst_c6_585(i,:)),diag(u10_c6_585),res1)
    overlay(plot3(i),plot4(i))
    draw(plot3(i))
    res1@tiYAxisString = ""
  end do

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) MC",x1+0*(1+a)*dx,0.89,restx)
  gsn_text_ndc(wks,"(b) IO - MC",x1+1*(1+a)*dx,0.89,restx)
  gsn_text_ndc(wks,"(c) WP - MC",x1+2*(1+a)*dx,0.89,restx)
  gsn_text_ndc(wks,"(d) WP - EP",x1+3*(1+a)*dx,0.89,restx)
  gsn_text_ndc(wks,"(e) MC",x1+0*(1+a)*dx,0.51,restx)
  gsn_text_ndc(wks,"(f) IO - MC",x1+1*(1+a)*dx,0.51,restx)
  gsn_text_ndc(wks,"(g) WP - MC",x1+2*(1+a)*dx,0.51,restx)
  gsn_text_ndc(wks,"(h) WP - EP",x1+3*(1+a)*dx,0.51,restx)

  restx@txFontHeightF = 0.016
  gsn_text_ndc(wks,"Relationship b/w background monsoon westerlies & area-averaged relative SST indices",x1,0.96,restx)
  restx@txFontHeightF = 0.014
  gsn_text_ndc(wks,"d4PDF 2K/4K",x1,0.92,restx)
  gsn_text_ndc(wks,"CMIP6 ("+M+" models, single member)",x1,0.54,restx)

  resm = True
  restx@txFontHeightF = 0.007
  restx@txJust = "CenterLeft"
  do i = 0, dimsizes(M2)-1
    ym = 0.220 - 0.018*(i%4)
    xm = 0.080 + (i/4)*0.115
    resm@gsMarkerColor = m_colors(i,:)
    resm@gsMarkerSizeF = m_sizes(i)
    resm@gsMarkerIndex = m_markers_2K(i)
    gsn_polymarker_ndc(wks,xm,ym,resm)
    resm@gsMarkerIndex = m_markers_4K(i)
    gsn_polymarker_ndc(wks,xm+0.011,ym,resm)
    gsn_text_ndc(wks,m_name(M2(i)),xm+0.020,ym,restx)
  end do
  restx@txFontHeightF = 0.007
  restx@txJust = "BottomCenter"
  do i = 0, dimsizes(M2)/4
    gsn_text_ndc(wks,"2K~S~ ~N~4K",0.085+0.115*i,0.23,restx)
  end do
    
  restx@txFontHeightF = 0.008
  restx@txJust = "CenterLeft"
  mtext = (/"All","CC","GF","HA","MI","MP","MR"/)
  do i = 0, 13
    resm@gsMarkerIndex := where(i%7.eq.0,16,i%7+3)
    resm@gsMarkerColor := colors((i/7)+1)
    resm@gsMarkerSizeF := where(i%7.eq.0,0.012,0.009)
    xm = 0.170+(i%7)*0.04
    ym = 0.615-(i/7)*0.02
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