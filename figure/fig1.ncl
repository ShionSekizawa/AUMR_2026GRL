;;; File: fig1.ncl
;;; Description: create Figure 1
;;;
;;; Input: clim_seasonality.nc created by d4PDF/clim_seasonality.ncl
;;;        AGCD data (land-sea mask file)
;;;        ausm_variance.nc created by d4PDF/derive_ausm_variance.ncl
;;; Output: fig1.eps
;;;
;;; Running:
;;; ncl fig1.ncl
;;;

;;; function for calculating dSST-induced correlation
function dSST_cor(x[12]:numeric,y[12]:numeric)
begin
  x(:5) = dim_rmvmean(x(:5))
  x(6:) = dim_rmvmean(x(6:))
  y(:5) = dim_rmvmean(y(:5))
  y(6:) = dim_rmvmean(y(6:))
  return(escorc(x,y))
end

begin
;;; open workstation
  wks = gsn_open_wks("eps","fig1")
  
;;; Panel (a): Rainfall climatology and Australia map
  in = addfile("../d4PDF/clim_seasonality.nc","r")
  mn = (/7,8,9,10,11,0,1,2,3,4,5,6/)
  mn_label = (/"Aug","Sep","Oct","Nov","Dec","Jan", \
               "Feb","Mar","Apr","May","Jun","Jul"/)
  obs = in->rain_obs(mn)
  l_obs = obs - in->err_obs(mn)
  u_obs = obs + in->err_obs(mn)
  d4 = in->rain_d4pdf(mn)
  l_d4 = d4 - in->err_d4pdf(mn)
  u_d4 = d4 + in->err_d4pdf(mn)

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@tmBorderThicknessF  = 1.5
  res@tmYLMajorThicknessF = 1.5
  res@tmXBMajorThicknessF = 1.5
  res@tmYLMinorThicknessF = 1.0
  res@tmXBMinorThicknessF = 1.0
  res@tmXBMajorLengthF = 0.006
  res@tmYLMajorLengthF = 0.006
  res@tmXBMinorLengthF = 0.003
  res@tmYLMinorLengthF = 0.003
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMajorOutwardLengthF = 0
  res@tmXBMinorOutwardLengthF = 0
  res@tmYLMinorOutwardLengthF = 0
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.011
  x = fspan(0,11,12)
  dx = 0.033
  res@trXMinF     =  0.-dx
  res@trXMaxF     = 11.+dx
  res@trYMinF     =  0.
  res@trYMaxF     = 10.5
  res@tmYROn      = False
  res@tmXTOn      = False
  res@tmXBMode    = "Explicit"
  res@tmXBValues  = x
  res@tmXBLabels  = mn_label
  res@tmYLMode    = "Manual"
  res@tmYLTickSpacingF  = 3
  res@tmYLTickStartF    = 0
  res@tmYLMinorPerMajor = 2
  res@xyLineThicknessF = 2.5
  res@xyDashPattern = 0
  res@tiYAxisFontHeightF = 0.011

  res@vpXF = 0.06
  res@vpYF = 0.95
  res@vpWidthF  = 0.66
  res@vpHeightF = 0.20
  res@xyLineColor = "black"
  res@xyDashPattern = 2
  res@tiYAxisString = "Rainfall [mm/day]"
  plot11 = gsn_csm_xy(wks,x-dx,obs,res)
  resjf = True
  resjf@tfPolyDrawOrder = "PreDraw"
  resjf@gsFillColor = (/1.00,0.90,0.97/)
  shade = gsn_add_polygon(wks,plot11,(/4.5,6.5,6.5,4.5/),(/0,0,12,12/),resjf)
  
  res@xyLineColor = "red"
  res@xyDashPattern = 0
  plot12 = gsn_csm_xy(wks,x+dx,d4,res)

  res_err = True
  res_err@tfPolyDrawOrder = "Draw"
  res_err@gsLineThicknessF = 1.5
  err_obs = new((/3,12/),graphic)
  err_d4 = new((/3,12/),graphic)
  dx_err = 0.033
  res_err@gsLineColor = "black"
  res_err@gsLineThicknessF = 1.5
  do i = 0, 11
    x_obs = x(i)-dx
    err_obs(0,i) = gsn_add_polyline(wks,plot11,(/x_obs,x_obs/),(/l_obs(i),u_obs(i)/)>0,res_err)
    err_obs(1,i) = gsn_add_polyline(wks,plot11,(/x_obs-dx_err,x_obs+dx_err/),(/l_obs(i),l_obs(i)/)>0,res_err)
    err_obs(2,i) = gsn_add_polyline(wks,plot11,(/x_obs-dx_err,x_obs+dx_err/),(/u_obs(i),u_obs(i)/)>0,res_err)
  end do
  overlay(plot11,plot12)
  res_err@gsLineColor = "red"
  do i = 0, 11
    x_d4 = x(i)+dx
    err_d4(0,i) = gsn_add_polyline(wks,plot11,(/x_d4,x_d4/),(/l_d4(i),u_d4(i)/),res_err)
    err_d4(1,i) = gsn_add_polyline(wks,plot11,(/x_d4-dx_err,x_d4+dx_err/),(/l_d4(i),l_d4(i)/),res_err)
    err_d4(2,i) = gsn_add_polyline(wks,plot11,(/x_d4-dx_err,x_d4+dx_err/),(/u_d4(i),u_d4(i)/),res_err)
  end do
  draw(plot11)

  ;; read AGCD land-sea mask file
  in = addfile("${AGCD_dir_path}/v2/land_fixed_0.05.nc","r")
  land = short2flt(in->land)
  land({-20:},:) = land({-20:},:)*2.0
  land = where(land.eq.1,land@_FillValue,land)
  res_map = True
  res_map@gsnFrame = False
  res_map@gsnDraw  = False
  res_map@cnLinesOn      = False
  res_map@cnLineLabelsOn = False
  res_map@lbLabelBarOn   = False
  res_map@cnFillOn       = True
  res_map@cnFillMode = "CellFill"
  res_map@cnInfoLabelOn  = False
  res_map@mpMinLatF = -45.
  res_map@mpMaxLatF = -8.
  res_map@mpMinLonF = 109.
  res_map@mpMaxLonF = 156.
  res_map@gsnAddCyclic = False
  res_map@vpYF = 0.95
  res_map@vpHeightF = 0.20
  res_map@vpWidthF = res_map@vpHeightF*(res_map@mpMaxLonF-res_map@mpMinLonF)/(res_map@mpMaxLatF-res_map@mpMinLatF)
  res_map@vpXF = 0.985-res_map@vpWidthF
  res_map@gsnLeftString  = ""
  res_map@tmBorderThicknessF = 1.5
  res_map@tmXBOn = False
  res_map@tmXTOn = False
  res_map@tmYLOn = False
  res_map@tmYROn = False
  res_map@cnLevelSelectionMode = "ExplicitLevels"
  res_map@cnLevels = (/1.5/)
  res_map@cnFillPalette = (/"white","magenta"/)
  plot13 = gsn_csm_contour_map_ce(wks,land,res_map)
  draw(plot13)
  delete(res)

;;; Panel (b, c): AUMR variance and mean
  in = addfile("../d4PDF/ausm_variance.nc","r")
  var = in->var
  p = new((/4,15/),float)
  do i = 1, 14
    p(0,i) = ttest(var(0,0,i),var(0,3,i),var(3,0,i),var(0,0,0),var(0,3,0),var(3,0,0),True,False)>0.0001
    if (var(0,0,i).lt.var(0,0,0)) then
      p(0,i) = -p(0,i)
    end if
    do j = 1, 3
      p(j,i) = ftest(var(0,j,i),var(3,j,i),var(0,j,0),var(3,j,0),0)>0.0001
      if (var(0,j,i).lt.var(0,j,0)) then
        p(j,i) = -p(j,i)
      end if
    end do
  end do
  
  x := fspan(0,14,15)
  x(3:8)  = x(3:8)  + 0.5
  x(9:14) = x(9:14) + 1.0
  dx = 0.25
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@tmBorderThicknessF  = 1.5
  res@tmYLMajorThicknessF = 1.5
  res@tmXBMajorThicknessF = 1.5
  res@tmYLMinorThicknessF = 1.0
  res@tmXBMinorThicknessF = 1.0
  res@tmXBMajorLengthF = 0.006
  res@tmYLMajorLengthF = 0.006
  res@tmXBMinorLengthF = 0.003
  res@tmYLMinorLengthF = 0.003
  res@tmXBMajorOutwardLengthF = 0
  res@tmYLMajorOutwardLengthF = 0
  res@tmXBMinorOutwardLengthF = 0
  res@tmYLMinorOutwardLengthF = 0
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.011
  res@trXMinF     = min(x)-0.75
  res@trXMaxF     = max(x)+0.75
  res@trYMinF     =  0.
  res@trYMaxF     =  3.6
  res@tmYROn      = False
  res@tmXTOn      = False
  res@tmXBMode    = "Explicit"
  res@tmXBValues  = x
  res@tmXBLabels  = (/"HIST","2K","4K", \
                      "CC","GF","HA","MI","MP","MR", \
                      "CC","GF","HA","MI","MP","MR" /)
  res@gsnXYBarChart         = True
  res@gsnXYBarChartBarWidth = 0.8*dx
  res@tmYLMode    = "Manual"
  res@tmYLTickSpacingF  = 1
  res@tmYLTickStartF    = 0
  res@tmYLMinorPerMajor = 4
  res@tiYAxisString = "AUMR variance [(mm/day)~S~2~N~]"
  res@tiYAxisFontHeightF = 0.011
  colors = (/"gray50","orange","firebrick"/)
  res@gsnXYBarChartColors = colors((/0,1,2,1,1,1,1,1,1,2,2,2,2,2,2/))
  res@vpXF = 0.06
  res@vpYF = 0.65
  res@vpWidthF  = 0.925
  res@vpHeightF = 0.225
  res@gsnXYBarChartFillScaleF = 0.4
  res@gsnXYBarChartPatterns = 0
  plot21 = gsn_csm_xy(wks,x-dx,var(0,3,:),res)
  res@gsnXYBarChartPatterns = 11
  plot22 = gsn_csm_xy(wks,x,var(0,1,:),res)
  res@gsnXYBarChartPatterns = 10
  plot23 = gsn_csm_xy(wks,x+dx,var(0,2,:),res)
  overlay(plot21,plot22)
  overlay(plot21,plot23)
  err = new((/3,3,15/),graphic)
  dx_err = 0.025
  res_err@gsLineColor = "black"
  do i = 0, 14
    xx = x(i)-dx
    err(0,0,i) = gsn_add_polyline(wks,plot21,(/xx,xx/),(/var(1,3,i),var(2,3,i)/),res_err)
    err(1,0,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(1,3,i),var(1,3,i)/),res_err)
    err(2,0,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(2,3,i),var(2,3,i)/),res_err)
    xx = x(i)
    err(0,1,i) = gsn_add_polyline(wks,plot21,(/xx,xx/),(/var(1,1,i),var(2,1,i)/),res_err)
    err(1,1,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(1,1,i),var(1,1,i)/),res_err)
    err(2,1,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(2,1,i),var(2,1,i)/),res_err)
    xx = x(i)+dx
    err(0,2,i) = gsn_add_polyline(wks,plot21,(/xx,xx/),(/var(1,2,i),var(2,2,i)/),res_err)
    err(1,2,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(1,2,i),var(1,2,i)/),res_err)
    err(2,2,i) = gsn_add_polyline(wks,plot21,(/xx-dx_err,xx+dx_err/),(/var(2,2,i),var(2,2,i)/),res_err)
  end do
  resl = True
  resl@gsLineThicknessF = 1.5
  x23 = 0.5*(x(2)+x(3))
  x89 = 0.5*(x(8)+x(9))
  line1 = gsn_add_polyline(wks,plot21,(/x23,x23/),(/0,3.6/),resl)
  line2 = gsn_add_polyline(wks,plot21,(/x89,x89/),(/0,3.6/),resl)
  resm = True
  resm@gsMarkerThicknessF = 1.0
  resm@gsMarkerSizeF = 0.005
  resm@gsMarkerIndex = 3
  sig = new((/3,15/),graphic)
  do j = 0, 2
    if (j.eq.0) then
      mi = 3
      mx = x-dx
    end if
    if (j.eq.1) then
      mi = 1
      mx = x
    end if
    if (j.eq.2) then
      mi = 2
      mx = x+dx
    end if
    do i = 1, 14
      if (abs(p(mi,i)).le.0.05) then
        if (p(mi,i).gt.0) then
          resm@gsMarkerColor = "magenta"
        else
          resm@gsMarkerColor = "skyblue"
        end if
        sig(j,i) = gsn_add_polymarker(wks,plot21,mx(i),res@trYMaxF*0.95,resm)
      end if
    end do
  end do
  draw(plot21)

  res@vpYF = 0.31
  res@trYMinF =  0.
  res@trYMaxF =  8.
  res@tmYLTickSpacingF  = 2
  res@tmYLTickStartF    = 0
  res@tmYLMinorPerMajor = 3
  res@tiYAxisString = "AUMR clim. [mm/day]"
  res@gsnXYBarChartBarWidth = res@gsnXYBarChartBarWidth*1.5
  res@gsnXYBarChartPatterns = 0
  plot31 = gsn_csm_xy(wks,x,var(0,0,:),res)
  dx_err = dx_err*1.5
  do i = 0, 14
    xx = x(i)
    err(0,1,i) = gsn_add_polyline(wks,plot31,(/xx,xx/),(/var(1,0,i),var(2,0,i)/),res_err)
    err(1,1,i) = gsn_add_polyline(wks,plot31,(/xx-dx_err,xx+dx_err/),(/var(1,0,i),var(1,0,i)/),res_err)
    err(2,1,i) = gsn_add_polyline(wks,plot31,(/xx-dx_err,xx+dx_err/),(/var(2,0,i),var(2,0,i)/),res_err)
  end do
  line1 = gsn_add_polyline(wks,plot31,(/x23,x23/),(/0,8/),resl)
  line2 = gsn_add_polyline(wks,plot31,(/x89,x89/),(/0,8/),resl)
  mi = 0
  mx = x
  do i = 1, 14
    if (abs(p(mi,i)).le.0.05) then
      if (p(mi,i).gt.0) then
        resm@gsMarkerColor = "magenta"
      else
        resm@gsMarkerColor = "skyblue"
      end if
      sig(0,i) = gsn_add_polymarker(wks,plot31,mx(i),res@trYMaxF*0.95,resm)
    end if
  end do
  draw(plot31)

  print("Total vs Inter-member variance: r="+dSST_cor(var(0,3,3:14),var(0,1,3:14)))
  print("Total vs SST-forced variance: r="+dSST_cor(var(0,3,3:14),var(0,2,3:14)))
  print("Mean vs Inter-member variance: r="+dSST_cor(var(0,0,3:14),var(0,1,3:14)))
  print("Mean vs SST-forced variance: r="+dSST_cor(var(0,0,3:14),var(0,2,3:14)))
  print("Mean vs Total variance: r="+dSST_cor(var(0,0,3:14),var(0,3,3:14)))

  restx = True
  restx@txJust = "BottomLeft"
  restx@txFontHeightF = 0.018
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(a) Climatological mean rainfall over northern Australia",0.06,0.96,restx)
  restx@txFontHeightF = 0.013
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"Observation (AGCD)",0.53,0.915,restx)
  restx@txFontColor   = "red"
  gsn_text_ndc(wks,"d4PDF HIST",0.53,0.89,restx)
  restx@txFontHeightF = 0.018
  restx@txFontColor   = "black"
  gsn_text_ndc(wks,"(b) AUMR variance in the d4PDF",0.06,0.66,restx)
  gsn_text_ndc(wks,"(c) Climatological-mean AUMR in the d4PDF",0.06,0.32,restx)
  restx@txFontHeightF = 0.014
  xb = (/0.50,0.61,0.78/)
  gsn_text_ndc(wks,"Total",xb(0)+restx@txFontHeightF*1.5,0.66,restx)
  gsn_text_ndc(wks,"Inter-member",xb(1)+restx@txFontHeightF*1.5,0.66,restx)
  gsn_text_ndc(wks,"SST-forced",xb(2)+restx@txFontHeightF*1.5,0.66,restx)
  resb = True
  resb@gsFillIndex = 0
  resb@gsFillScaleF = 0.4
  resb@gsEdgesOn = True
  gsn_polygon_ndc(wks,xb(0)+(/0,restx@txFontHeightF,restx@txFontHeightF,0,0/),0.66+(/0,0,restx@txFontHeightF,restx@txFontHeightF,0/),resb)
  resb@gsFillIndex = 11
  gsn_polygon_ndc(wks,xb(1)+(/0,restx@txFontHeightF,restx@txFontHeightF,0,0/),0.66+(/0,0,restx@txFontHeightF,restx@txFontHeightF,0/),resb)
  resb@gsFillIndex = 10
  gsn_polygon_ndc(wks,xb(2)+(/0,restx@txFontHeightF,restx@txFontHeightF,0,0/),0.66+(/0,0,restx@txFontHeightF,restx@txFontHeightF,0/),resb)
  restx@txFontHeightF = 0.012
  restx@txJust = "TopCenter"
  gsn_text_ndc(wks,"2K",0.435,0.385,restx)
  gsn_text_ndc(wks,"4K",0.800,0.385,restx)
  gsn_text_ndc(wks,"2K",0.435,0.045,restx)
  gsn_text_ndc(wks,"4K",0.800,0.045,restx)
  gsn_polyline_ndc(wks,(/0.280,0.600/),(/0.395,0.395/),resl)
  gsn_polyline_ndc(wks,(/0.645,0.965/),(/0.395,0.395/),resl)
  gsn_polyline_ndc(wks,(/0.280,0.600/),(/0.055,0.055/),resl)
  gsn_polyline_ndc(wks,(/0.645,0.965/),(/0.055,0.055/),resl)
  
;  drawNDCGrid(wks)
  frame(wks)

end