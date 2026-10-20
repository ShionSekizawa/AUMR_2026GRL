;;; File: fig2.ncl
;;; Description: create Figure 2
;;;
;;; Input: ausm_rain_${experiment_type_name}.nc created by d4PDF/make_ausm_rain.ncl
;;;        MSEbudget_${experiment_type_name}_advect.nc created by d4PDF/derive_MSEbudget_advect.ncl
;;;        MSEbudget_${experiment_type_name}_forcing.nc created by d4PDF/derive_MSEbudget_forcing.ncl
;;;        MSEbudget_area_${experiment_type_name}.nc created by d4PDF/make_MSEbudget_area.ncl
;;;        lhf_area_${experiment_type_name}.nc created by d4PDF/make_lhf_area.ncl
;;;        ausm_variance.nc created by d4PDF/derive_ausm_variance.ncl
;;; Output: fig2.eps
;;;
;;; Running:
;;; ncl fig2.ncl
;;;

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

;;; function for calculating anomalies (for 4d array)
function anom4(x[*][*][*][*]:numeric,clim_opt[1]:integer)
local dim, n
begin
  dim = dimsizes(x)
  n   = dim(0)
  if (clim_opt.eq.6) then
    do i = 0, clim_opt-1
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:,:,:) = dim_rmvmean_n(x(i*n/clim_opt:(i+1)*n/clim_opt-1,:,:,:),(/0,1/))
    end do
  else
    x = dim_rmvmean_n(x,(/0,1/))
  end if
  x = dim_rmvmean_n(x,0)
  return(reshape(x,(/dim(0)*dim(1),dim(2),dim(3)/)))
end

;;; function for converting 1d array into diagonal 2d array for marker plotting
function diag(x[*]:numeric)
local N, i, y
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

begin
  ;; constants
  L = 2.507e6 
  scale = L/(60*60*24)
  g = 9.8

  ;; get information on the d4PDF experiments
  data_list = str_split_csv(asciiread("../d4PDF/data_list",14,"string")," ",0)
  expr = data_list(:,0)
  m1 = toint(data_list(:,1))
  m2 = toint(data_list(:,2))
  y1 = toint(data_list(:,3))
  y2 = toint(data_list(:,4))

  ;; array for the panels (a-d)
  sfc_map  = new((/49,97/),float)
  rad_map  = new((/49,97/),float)
  hadv_map = new((/49,97/),float)
  vadv_map = new((/49,97/),float)
  sfc_p  = new((/49,97/),float)
  rad_p  = new((/49,97/),float)
  hadv_p = new((/49,97/),float)
  vadv_p = new((/49,97/),float)

  ;; array for the panels (e-h)
  sfc_area = new((/15,3/),float)
  rad_area = new((/15,3/),float)
  hadv_area = new((/15,3/),float)
  vadv_area = new((/15,3/),float)

  ;; array for the panels (i-k)
  dEdW_area = new(15,float)
  dWdu_area = new(15,float)
  u10_area = new(15,float)

  ;; HIST experiment
  in1 = addfile("../d4PDF/ausm_rain_"+expr(0)+".nc","r")
  in2 = addfile("../d4PDF/MSEbudget_"+expr(0)+"_forcing.nc","r")
  in3 = addfile("../d4PDF/MSEbudget_"+expr(0)+"_advect.nc","r")
  in4 = addfile("../d4PDF/MSEbudget_area_"+expr(0)+".nc","r")
  in5 = addfile("../d4PDF/lhf_area_"+expr(0)+".nc","r")
  lat = in2->lat
  lon = in2->lon
  sfc_map = regCoef_n(anom2(in1->ausm,-1),anom4(in2->lhf+in2->shf,-1),0,0)
  sfc_p  = reshape(cdft_p(sfc_map@tval,sfc_map@nptxy-1),(/49,97/))
  rad_map = regCoef_n(anom2(in1->ausm,-1),anom4(in2->lwr+in2->swr,-1),0,0)
  rad_p = reshape(cdft_p(rad_map@tval,rad_map@nptxy-1),(/49,97/))
  hadv_map = regCoef_n(anom2(in1->ausm,-1),anom4(in3->shadv+L*in3->qhadv,-1),0,0)
  hadv_p = reshape(cdft_p(hadv_map@tval,hadv_map@nptxy-1),(/49,97/))
  vadv_map = regCoef_n(anom2(in1->ausm,-1),anom4(in3->svadv+L*in3->qvadv,-1),0,0)
  vadv_p = reshape(cdft_p(vadv_map@tval,vadv_map@nptxy-1),(/49,97/))
  sfc_area(0,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->lhf+in4->shf,-1),0,0)
  rad_area(0,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->lwr+in4->swr,-1),0,0)
  hadv_area(0,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->shadv+L*in4->qhadv,-1),0,0)
  vadv_area(0,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->svadv+L*in4->qvadv,-1),0,0)
  dEdW_area(0) = regCoef_n(anom2(in5->wspd,-1),anom2(in5->lhf,-1)/scale,0,0)
  dWdu_area(0) = regCoef_n(anom2(in5->u10,-1),anom2(in5->wspd,-1),0,0)
  u10_area(0) = regCoef_n(anom2(in1->ausm,-1),anom2(in5->u10,-1),0,0)
  delete([/in1,in2,in3,in4,in5/])

  ;; 2K experiments as a whole
  in1 = addfiles("../d4PDF/ausm_rain_"+expr(2:7)+".nc","r")
  in4 = addfiles("../d4PDF/MSEbudget_area_"+expr(2:7)+".nc","r")
  in5 = addfiles("../d4PDF/lhf_area_"+expr(2:7)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in4,"cat")
  ListSetType(in5,"cat")
  sfc_area(1,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->lhf+in4[:]->shf,6),0,0)
  rad_area(1,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->lwr+in4[:]->swr,6),0,0)
  hadv_area(1,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->shadv+L*in4[:]->qhadv,6),0,0)
  vadv_area(1,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->svadv+L*in4[:]->qvadv,6),0,0)
  dEdW_area(1) = regCoef_n(anom2(in5[:]->wspd,6),anom2(in5[:]->lhf,6)/scale,0,0)
  dWdu_area(1) = regCoef_n(anom2(in5[:]->u10,6),anom2(in5[:]->wspd,6),0,0)
  u10_area(1) = regCoef_n(anom2(in1[:]->ausm,6),anom2(in5[:]->u10,6),0,0)

  ;; 4K experiments as a whole
  in1 = addfiles("../d4PDF/ausm_rain_"+expr(8:13)+".nc","r")
  in4 = addfiles("../d4PDF/MSEbudget_area_"+expr(8:13)+".nc","r")
  in5 = addfiles("../d4PDF/lhf_area_"+expr(8:13)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in4,"cat")
  ListSetType(in5,"cat")
  sfc_area(2,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->lhf+in4[:]->shf,6),0,0)
  rad_area(2,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->lwr+in4[:]->swr,6),0,0)
  hadv_area(2,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->shadv+L*in4[:]->qhadv,6),0,0)
  vadv_area(2,:) = regCoef_n(anom2(in1[:]->ausm,6),anom3(in4[:]->svadv+L*in4[:]->qvadv,6),0,0)
  dEdW_area(2) = regCoef_n(anom2(in5[:]->wspd,6),anom2(in5[:]->lhf,6)/scale,0,0)
  dWdu_area(2) = regCoef_n(anom2(in5[:]->u10,6),anom2(in5[:]->wspd,6),0,0)
  u10_area(2) = regCoef_n(anom2(in1[:]->ausm,6),anom2(in5[:]->u10,6),0,0)
  delete([/in1,in4,in5/])

  do i = 2, 13
    ;; 2K (i=2-7), 4K (i=8-13) experiments
    in1 = addfile("../d4PDF/ausm_rain_"+expr(i)+".nc","r")
    in4 = addfile("../d4PDF/MSEbudget_area_"+expr(i)+".nc","r")
    in5 = addfile("../d4PDF/lhf_area_"+expr(i)+".nc","r")
    sfc_area(i+1,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->lhf+in4->shf,-1),0,0)
    rad_area(i+1,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->lwr+in4->swr,-1),0,0)
    hadv_area(i+1,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->shadv+L*in4->qhadv,-1),0,0)
    vadv_area(i+1,:) = regCoef_n(anom2(in1->ausm,-1),anom3(in4->svadv+L*in4->qvadv,-1),0,0)
    dEdW_area(i+1) = regCoef_n(anom2(in5->wspd,-1),anom2(in5->lhf,-1)/scale,0,0)
    dWdu_area(i+1) = regCoef_n(anom2(in5->u10,-1),anom2(in5->wspd,-1),0,0)
    u10_area(i+1) = regCoef_n(anom2(in1->ausm,-1),anom2(in5->u10,-1),0,0)
  end do
  delete([/in1,in4,in5/])

  in = addfile("../d4PDF/ausm_variance.nc","r")
  var = in->var(0,1,:)

  sfc_map!0 = "lat"
  sfc_map!1 = "lon"
  sfc_map&lat = lat
  sfc_map&lon = lon
  copy_VarCoords(sfc_map,rad_map)
  copy_VarCoords(sfc_map,hadv_map)
  copy_VarCoords(sfc_map,vadv_map)
  copy_VarCoords(sfc_map,sfc_p)
  copy_VarCoords(sfc_map,rad_p)
  copy_VarCoords(sfc_map,hadv_p)
  copy_VarCoords(sfc_map,vadv_p)

;;; open workstation
  wks = gsn_open_wks("pdf","fig2")

;;; Panels (a-d)
  plot1 = new(4,graphic)
  plot2 = new(4,graphic)
  box = new(4,graphic)

  res1                 = True
  res1@gsnDraw         = False
  res1@gsnFrame        = False
  res1@cnLinesOn       = False
  res1@cnFillOn        = True
  res1@cnInfoLabelOn   = False
  res1@lbLabelBarOn    = False
  res1@cnLineLabelsOn  = False
  res1@cnFillPalette   = "ncl_default"
  res1@mpMaxLatF       =   0.
  res1@mpMinLatF       = -45.
  res1@mpMaxLonF       = 165.
  res1@mpMinLonF       =  90.
  res1@mpCenterLonF    = 135.
  res1@gsnAddCyclic    = False
  res1@tmBorderThicknessF  = 1.25
  res1@tmYLMajorThicknessF = 1.25
  res1@tmXBMajorThicknessF = 1.25
  res1@tmXBMajorLengthF = 0.004
  res1@tmYLMajorLengthF = 0.004
  res1@tmXBMajorOutwardLengthF = 0
  res1@tmYLMajorOutwardLengthF = 0
  res1@tmYLMinorOutwardLengthF = 0
  res1@tmXBMinorOutwardLengthF = 0
  res1@tmXTOn = False
  res1@tmYROn = False
  res1@tmYLLabelFontHeightF = 0.008
  res1@tmXBLabelFontHeightF = 0.008
  res1@lbBoxEndCapStyle     = "TriangleBothEnds"
  res1@lbLabelFontHeightF   = 0.009
  res1@lbTitleOn            = True
  res1@lbTitleOffsetF       = -0.01
  res1@lbTitleString        = "[(W/m~S~2~N~)/(mm/day)]"
  res1@lbTitlePosition      = "Right"
  res1@lbTitleDirection     = "Across" 
  res1@lbTitleFontHeightF   = 0.009
  res1@pmLabelBarWidthF  = 0.88
  res1@pmLabelBarHeightF = 0.03
  res1@pmLabelBarOrthogonalPosF = 0.2

  res2                 = True
  res2@gsnDraw         = False
  res2@gsnFrame        = False
  res2@cnLinesOn       = False
  res2@cnLineLabelsOn  = False
  res2@cnFillOn        = False
  res2@cnInfoLabelOn   = False
  res2@cnFillColor     = "gray70"
  res2@cnFillScaleF    = 0.4
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLevels = (/0.0249,0.025,0.0251,0.9749,0.975,0.9751/)

  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMinLevelValF  = -30.0
  res1@cnMaxLevelValF  =  30.0
  res1@cnLevelSpacingF =  3.

  resbx = True
  resbx@gsLineColor = "magenta"
  resbx@gsLineThicknessF = 2
  resbx@tfPolyDrawOrder = "PostDraw"
  s1 = -20
  s2 =  -6
  s3 = 105
  s4 = 150

  res1@vpYF = 0.9
  x1 = 0.05
  x2 = 0.985
  a  = 0.20
  dx = (x2-x1)/(4+3*a)
  res1@vpHeightF = dx*(res1@mpMaxLatF-res1@mpMinLatF)/ \
                      (res1@mpMaxLonF-res1@mpMinLonF)
  VAR_map = [/sfc_map,rad_map,hadv_map,vadv_map/]
  VAR_p   = [/  sfc_p,  rad_p,  hadv_p,  vadv_p/]
  do i = 0, 3
    if (i.eq.3) then
      res1@lbLabelBarOn = True
      res1@pmLabelBarParallelPosF = -1-1.5*a
    end if
    res1@vpXF = x1+i*(1+a)*dx
    res1@vpWidthF = dx
    plot1(i) = gsn_csm_contour_map_ce(wks,VAR_map[i],res1)
    plot2(i) = gsn_csm_contour(wks,VAR_p[i],res2)
    plot2(i) = ShadeGeLeContour(plot2(i),0.025,0.975,6)
    overlay(plot1(i),plot2(i))
    box(i) = gsn_add_polyline(wks,plot1(i),(/s3,s4,s4,s3,s3/),(/s1,s1,s2,s2,s1/),resbx)
    draw(plot1(i))
  end do

;;; Panels (e-h)
  plot3 = new((/4,2/),graphic) 
  res3 = True
  res3@gsnDraw  = False
  res3@gsnFrame = False
  res3@tmBorderThicknessF  = 1.25
  res3@tmYLMajorThicknessF = 1.25
  res3@tmXBMajorThicknessF = 1.25
  res3@tmXBMajorLengthF = 0.005
  res3@tmYLMajorLengthF = 0.005
  res3@tmXBMinorLengthF = 0.002
  res3@tmYLMinorLengthF = 0.002
  res3@tmXBMajorOutwardLengthF = 0.
  res3@tmYLMajorOutwardLengthF = 0.
  res3@tmYLMinorOutwardLengthF = 0.
  res3@tmXBMinorOutwardLengthF = 0.
  res3@tmXBLabelFontHeightF = 0.008
  res3@tmYLLabelFontHeightF = 0.008
  res3@tiYAxisFontHeightF = 0.008
  res3@tiXAxisFontHeightF = 0.008
  res3@tmYROn = False
  res3@tmXTOn = False
  res3@xyMarkLineMode = "Markers"
  res3@xyMarkerThicknessF = 1.0
  res3@tmYLMode = "Manual"
  res3@gsnYRefLine = 0
  res3@trXMinF = 1.2
  res3@trXMaxF = 2.7
  colors = (/"black","orange","firebrick"/)
  res3@vpYF = 0.68
  res3@vpWidthF  = dx
  res3@vpHeightF = dx
  VAR_area3 = [/sfc_area,rad_area,hadv_area,vadv_area/]
  res3@tiXAxisString = "AUMR variance [(mm/day)~S~2~N~]"
  res3@gsnLeftString = ""
  res3@gsnRightStringOrthogonalPosF = -0.12
  res3@gsnRightStringParallelPosF = 0.96
  res3@gsnRightStringFontHeightF = 0.010
  ypre = (/2,2,2,2/)
  aid = (/2,2,2,2/)
  y1 := (/ 0., 0.,-2.5,-5./)
  y2 := (/ 5., 5., 2.5, 0./)
  res3@tmYLTickSpacingF = 1.
  res3@tiYAxisString = "MSE tendency [(W/m~S~2~N~)/(mm/day)]"
  do i = 0, 3
    res3@vpXF = x1+i*(1+a)*dx
    res3@trYMinF = y1(i)
    res3@trYMaxF = y2(i)
    res3@tmYLPrecision = ypre(i)
    res3@tmYLTickStartF = ceil(y1(i))
    res3@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
    res3@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
    res3@xyMarkerSizeF = 0.009
    printMinMax(VAR_area3[i](3:14,aid(i)),0)
    temp = dSST_cor(var(3:14),VAR_area3[i](3:14,aid(i)))
    res3@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
    plot3(i,0) = gsn_csm_xy(wks,diag(var(3:14)),diag(VAR_area3[i](3:14,aid(i))),res3)
    res3@xyMarkerColors := colors((/0,1,2/))
    res3@xyMarkers := (/16,16,16/)
    res3@xyMarkerSizeF = 0.012
    res3@gsnRightString = ""
    res3@tiYAxisString = ""
    plot3(i,1) = gsn_csm_xy(wks,diag(var(0:2)),diag(VAR_area3[i](0:2,aid(i))),res3)
    overlay(plot3(i,0),plot3(i,1))
    draw(plot3(i,0))
  end do
  delete(res3@tmYLTickStartF)

;;; Panels (i-k)
  plot4 = new((/3,2/),graphic)
  res4 = res3
  res4@vpYF = 0.35
  ypre = (/1,1,1,1/)
  y1 := (/ 0.50, 0.45, 0.20/)
  y2 := (/ 0.80, 0.75, 0.50/)
  Spacing = (/0.10,0.10,0.10/)
  VAR_area4 = [/dEdW_area,dWdu_area,u10_area/]
  d = "~F34~6~F~"
  YString = (/d+"E/"+d+"W [(mm/day)/(m/s)]",d+"W/"+d+"u~B~s~N~ [(m/s)/(m/s)]","10-m zonal wind anom. [(m/s)/(mm/day)]"/)
  res4@gsnRightStringOrthogonalPosF = 0.02
  res4@gsnRightStringParallelPosF = 0.98
  res5 = True
  res5@tmBorderThicknessF  = 1.25
  res5@tmYRMajorThicknessF = 1.25
  res5@tmYRMajorLengthF = 0.005
  res5@tmYRMinorLengthF = 0.002
  res5@tmYRMajorOutwardLengthF = 0.
  res5@tmYRMinorOutwardLengthF = 0.
  res5@tmYRLabelFontHeightF = 0.008
  res5@tiYAxisFontHeightF = 0.008
  res5@tiYAxisString = "Ratio relative to HIST [%]"
  res5@tmYRMode = "Explicit"
  res5@tmYRPrecision = 2
  YRLabels = ispan(0,200,5)
  res5@tmYRValues = YRLabels(ind(YRLabels%20.eq.0))
  res5@tmYRMinorValues = YRLabels(ind(YRLabels%20.ne.0))
  res5@tmYRLabels = res5@tmYRValues
  do i = 0, 2
    res4@vpXF = x1+i*(1+3.*a)*dx
    res4@trYMinF = y1(i)
    res4@trYMaxF = y2(i)
    res4@tmYLTickSpacingF = Spacing(i)
    res4@tiYAxisString = YString(i)
    res4@tmYLPrecision = ypre(i)
    res4@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
    res4@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
    res4@xyMarkerSizeF = 0.009
    printMinMax(VAR_area4[i](3:14),0)
    dummy = diag(VAR_area4[i](3:14))
    res5@trYMinF = 100*y1(i)/VAR_area4[i](0)
    res5@trYMaxF = 100*y2(i)/VAR_area4[i](0)
    temp = dSST_cor(var(3:14),VAR_area4[i](3:14))
    res4@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",temp(0))+where(temp@stat,"~F18S~M","")
    plot4(i,0) = gsn_csm_xy2(wks,diag(var(3:14)),diag(VAR_area4[i](3:14)),dummy,res4,res5)
    res4@xyMarkerColors := colors((/0,1,2/))
    res4@xyMarkers := (/16,16,16/)
    res4@xyMarkerSizeF = 0.012
    res4@gsnRightString = ""
    plot4(i,1) = gsn_csm_xy(wks,diag(var(0:2)),diag(VAR_area4[i](0:2)),res4)
    overlay(plot4(i,0),plot4(i,1))
    draw(plot4(i,0))
  end do

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) Surface heat flux",x1+0*(1+a)*dx,res1@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(b) Radiation",x1+1*(1+a)*dx,res1@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(c) Horizontal advection",x1+2*(1+a)*dx,res1@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(d) Vertical advection",x1+3*(1+a)*dx,res1@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(e) Surface heat flux",x1+0*(1+a)*dx,res3@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(f) Radiation",x1+1*(1+a)*dx,res3@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(g) Horizontal advection",x1+2*(1+a)*dx,res3@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(h) Vertical advection",x1+3*(1+a)*dx,res3@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(i) "+d+"E/"+d+"W",x1+0*(1+3.*a)*dx,res4@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(j) "+d+"W/"+d+"u~B~s",x1+1*(1+3.*a)*dx,res4@vpYF+0.01,restx)
  gsn_text_ndc(wks,"(k) u~B~s~N~'",x1+2*(1+3.*a)*dx,res4@vpYF+0.01,restx)

  restx@txFontHeightF = 0.016
  gsn_text_ndc(wks,"MSE budget for JF-mean regression anomalies onto AUMR",x1,res1@vpYF+0.04,restx)
  gsn_text_ndc(wks,"Decomposition of surface evaporation anomalies",x1,res4@vpYF+0.04,restx)

  resm = True
  restx@txFontHeightF = 0.009
  restx@txJust = "CenterLeft"
  mtext = (/"All","CC","GF","HA","MI","MP","MR"/)
  do i = 0, 13
    resm@gsMarkerIndex = where(i%7.eq.0,16,i%7+3)
    resm@gsMarkerColor = colors((i/7)+1)
    resm@gsMarkerSizeF = where(i%7.eq.0,0.012,0.009)
    xm = 0.18+(i%7)*0.042+(i/7)*0.36
    ym = 0.08
    gsn_polymarker_ndc(wks,xm,ym,resm)
    gsn_text_ndc(wks,mtext(i%7),xm+0.012,ym,restx)
    if (i%7.eq.0) then
      gsn_text_ndc(wks,where(i/7.eq.0,"2K :","4K :"),xm-0.036,ym,restx)
    end if
  end do
  xm = 0.1
  gsn_text_ndc(wks,"HIST :",xm-0.05,ym,restx)
  resm@gsMarkerIndex = 16
  resm@gsMarkerColor = colors(0)
  resm@gsMarkerSizeF = 0.012
  gsn_polymarker_ndc(wks,xm,ym,resm)

;  drawNDCGrid(wks)
  frame(wks)
  
end