;;; File: figS3.ncl
;;; Description: create Figure S3
;;;
;;; Input: ausm_rain_${experiment_type_name}.nc created by d4PDF/make_ausm_rain.ncl
;;;        vadv_area_${experiment_type_name}.nc created by d4PDF/make_vadv_area.ncl
;;; Output: figS3.eps
;;;
;;; Running:
;;; ncl figS3.ncl
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

;;; function for calculating dSST-induced correlation
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

begin
  ;; constants
  L = 2.507e6
  scale = L/(60*60*24)
  g = 9.8

  ;; get information on the d4PDF experiments
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(:,0)
  m1 = toint(data_list(:,1))
  m2 = toint(data_list(:,2))
  y1 = toint(data_list(:,3))
  y2 = toint(data_list(:,4))

  dhc_area = new(15,float)  ;; array for the panel (a)
  wa_area  = new(15,float)  ;; array for the panel (b)
  dha_area = new(15,float)  ;; array for the panel (x)
  wc_area  = new(15,float)  ;; array for the panel (x)

;;; d4PDF data
  ;; HIST experiment
  in1 = addfile("../d4PDF/ausm_rain_"+expr(0)+".nc","r")
  in2 = addfile("../d4PDF/vadv_area_"+expr(0)+".nc","r")
  dhc_area(0) = avg((in2->su+L*in2->qu)-(in2->sl+L*in2->ql))
  wa_area(0) = regCoef_n(anom2(in1->ausm,-1),anom2(in2->wm,-1),0,0)
  dha_area(0) = regCoef_n(anom2(in1->ausm,-1),anom2((in2->su+L*in2->qu)-(in2->sl+L*in2->ql),-1),0,0)
  wc_area(0) = avg(in2->wm)
  delete([/in1,in2/])

  ;;; 2K experiments as a whole
  in1 = addfiles("../d4PDF/ausm_rain_"+expr(2:7)+".nc","r")
  in2 = addfiles("../d4PDF/vadv_area_"+expr(2:7)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in2,"cat")
  dhc_area(1) = avg((in2[:]->su+L*in2[:]->qu)-(in2[:]->sl+L*in2[:]->ql))
  wa_area(1) = regCoef_n(anom2(in1[:]->ausm,6),anom2(in2[:]->wm,6),0,0)
  dha_area(1) = regCoef_n(anom2(in1[:]->ausm,6),anom2((in2[:]->su+L*in2[:]->qu)-(in2[:]->sl+L*in2[:]->ql),6),0,0)
  wc_area(1) = avg(in2[:]->wm)

  ;; 4K experiments as a whole
  in1 = addfiles("../d4PDF/ausm_rain_"+expr(8:13)+".nc","r")
  in2 = addfiles("../d4PDF/vadv_area_"+expr(8:13)+".nc","r")
  ListSetType(in1,"cat")
  ListSetType(in2,"cat")
  dhc_area(2) = avg((in2[:]->su+L*in2[:]->qu)-(in2[:]->sl+L*in2[:]->ql))
  wa_area(2) = regCoef_n(anom2(in1[:]->ausm,6),anom2(in2[:]->wm,6),0,0)
  dha_area(2) = regCoef_n(anom2(in1[:]->ausm,6),anom2((in2[:]->su+L*in2[:]->qu)-(in2[:]->sl+L*in2[:]->ql),6),0,0)
  wc_area(2) = avg(in2[:]->wm)
  delete([/in1,in2/])
  
  do i = 2, 13
    ;; 2K (i=2-7), 4K (i=8-13) experiments
    in1 = addfile("../d4PDF/ausm_rain_"+expr(i)+".nc","r")
    in2 = addfile("../d4PDF/vadv_area_"+expr(i)+".nc","r")
    dhc_area(i+1) = avg((in2->su+L*in2->qu)-(in2->sl+L*in2->ql))
    wa_area(i+1) = regCoef_n(anom2(in1->ausm,-1),anom2(in2->wm,-1),0,0)
    dha_area(i+1) = regCoef_n(anom2(in1->ausm,-1),anom2((in2->su+L*in2->qu)-(in2->sl+L*in2->ql),-1),0,0)
    wc_area(i+1) = avg(in2->wm)
  end do
  delete([/in1,in2/])

  dhc_area = dhc_area/1000
  wa_area = -wa_area*100
  dha_area = -dha_area/1000
  wc_area = -wc_area*100

  in = addfile("../d4PDF/ausm_variance.nc","r")
  var = in->var(0,1,:)

;;; open workstation
  wks = gsn_open_wks("eps","figS3")
  x1 = 0.05
  x2 = 0.985
  a  = 0.2
  dx = (x2-x1)/(4+3*a)

  plot1 = new((/4,2/),graphic) 
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
  res1@tmXBMajorOutwardLengthF = 0.
  res1@tmYLMajorOutwardLengthF = 0.
  res1@tmYLMinorOutwardLengthF = 0.
  res1@tmXBMinorOutwardLengthF = 0.
  res1@tmXBLabelFontHeightF = 0.008
  res1@tmYLLabelFontHeightF = 0.008
  res1@tiYAxisFontHeightF = 0.008
  res1@tiXAxisFontHeightF = 0.008
  res1@tiXAxisString = "AUMR variance [(mm/day)~S~2~N~]"
  res1@tmYROn = False
  res1@tmXTOn = False
  res1@xyMarkLineMode = "Markers"
  res1@xyMarkerThicknessF = 1.0
  res1@tmYLMode = "Manual"
  res1@gsnYRefLine = 0
  res1@trXMinF = 1.2
  res1@trXMaxF = 2.7
  colors = (/"black","orange","firebrick"/)
  res1@vpWidthF  = dx
  res1@vpHeightF = dx

  ypre = (/2,2,2,2/)
  y1 := (/ 1., 0.6, 0., 10./)
  y2 := (/ 6., 1.0, 1., 40./)
  Spacing = (/1.0,0.1,0.20,5.0/)
  VAR_area = [/dhc_area,wa_area,dha_area,wc_area/]
  YString = (/"Clim. MSE diff. [kJ/m~S~2~N~]",\
              "-~F33~w~F~ anom. [(hPa/s)/(mm/day)]","",""/)
  res1@gsnRightStringOrthogonalPosF = 0.02
  res1@gsnRightStringParallelPosF = 0.98
  
  res2 = True
  res2@tmBorderThicknessF  = 1.25
  res2@tmYRMajorThicknessF = 1.25
  res2@tmYRMajorLengthF = 0.005
  res2@tmYRMinorLengthF = 0.002
  res2@tmYRMajorOutwardLengthF = 0.
  res2@tmYRMinorOutwardLengthF = 0.
  res2@tmYRLabelFontHeightF = 0.008
  res2@tiYAxisFontHeightF = 0.008
  res2@tiYAxisString = "Ratio relative to HIST [%]"
  res2@tmYRMode = "Explicit"
  res2@tmYRPrecision = 2
  YRLabels = ispan(0,200,5)
  res2@tmYRValues = YRLabels(ind(YRLabels%20.eq.0))
  res2@tmYRMinorValues = YRLabels(ind(YRLabels%20.ne.0))
  res2@tmYRLabels = res2@tmYRValues
  
  do i = 0, 1
    res1@vpXF = x1+i%2*(1+4*a)*dx
    res1@vpYF = 0.90-dx*1.25*(i/2)
    res1@trYMinF = y1(i)
    res1@trYMaxF = y2(i)
    res1@tmYLTickSpacingF = Spacing(i)
    res1@tiYAxisString = YString(i)
    res1@tmYLPrecision = ypre(i)
    res1@xyMarkerColors := colors((/1,1,1,1,1,1,2,2,2,2,2,2/))
    res1@xyMarkers := (/4,5,6,7,8,9,4,5,6,7,8,9/)
    res1@xyMarkerSizeF = 0.009
    printMinMax(VAR_area[i],0)
    dummy = diag(VAR_area[i](3:14))
    res2@trYMinF = 100*y1(i)/VAR_area[i](0)
    res2@trYMaxF = 100*y2(i)/VAR_area[i](0)
;    res1@gsnRightString = "r~BF33~D~F~SST~N~ = "+sprintf("%3.2f",dSST_cor(var(3:14),VAR_area[i](3:14)))
    plot1(i,0) = gsn_csm_xy2(wks,diag(var(3:14)),diag(VAR_area[i](3:14)),dummy,res1,res2)
    res1@xyMarkerColors := colors((/0,1,2/))
    res1@xyMarkers := (/16,16,16/)
    res1@xyMarkerSizeF = 0.012
    res1@gsnRightString = ""
    plot1(i,1) = gsn_csm_xy(wks,diag(var(0:2)),diag(VAR_area[i](0:2)),res1)
    overlay(plot1(i,0),plot1(i,1))
    draw(plot1(i,0))
  end do

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.012
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) Clim. MSE diff. (150hPa - 925hPa)",x1+0*(1+4*a)*dx,0.91-0*dx*1.25,restx)
  gsn_text_ndc(wks,"(b) -~F33~w~F~ anom. (500hPa)",x1+1*(1+4*a)*dx,0.91-0*dx*1.25,restx)

  restx@txFontHeightF = 0.015
  gsn_text_ndc(wks,"Decomposition of column-integrated vertical MSE advection anomalies",x1,0.94,restx)

  resm = True
  restx@txFontHeightF = 0.009
  restx@txJust = "CenterLeft"
  mtext = (/"All","CC","GF","HA","MI","MP","MR"/)
  do i = 0, 13
    resm@gsMarkerIndex = where(i%7.eq.0,16,i%7+3)
    resm@gsMarkerColor = colors((i/7)+1)
    resm@gsMarkerSizeF = where(i%7.eq.0,0.012,0.009)
    xm = 0.18+(i%7)*0.042
    ym = 0.64-(i/7)*0.02
    gsn_polymarker_ndc(wks,xm,ym,resm)
    gsn_text_ndc(wks,mtext(i%7),xm+0.012,ym,restx)
    if (i%7.eq.0) then
      gsn_text_ndc(wks,where(i/7.eq.0,"2K :","4K :"),xm-0.036,ym,restx)
    end if
  end do
  xm = 0.1
  ym = 0.64
  gsn_text_ndc(wks,"HIST :",xm-0.05,ym,restx)
  resm@gsMarkerIndex = 16
  resm@gsMarkerColor = colors(0)
  resm@gsMarkerSizeF = 0.012
  gsn_polymarker_ndc(wks,xm,ym,resm)

;  drawNDCGrid(wks)
  frame(wks)
  
end