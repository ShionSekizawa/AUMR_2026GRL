;;; File: figS2.ncl
;;; Description: create Figure S2
;;;
;;; Input: d4PDF/figS2_tmpfile_d4pdf.nc created by d4PDF/make_figS2_tmpfile_d4pdf.ncl
;;;        d4PDF/figS2_tmpfile_obs.nc created by obs/make_figS2_tmpfile_obs.ncl
;;; Output: figS2.eps
;;;
;;; Running:
;;; ncl figS2.ncl
;;;

begin

  n = 2010-1951+1-2
  t = cdft_t(0.975,n)
  r = t/sqrt(n+t^2)
  
  in = addfile("../obs/figS2_tmpfile_obs.nc","r")
  obs_std = in->obs_std
  obs_cor = in->obs_cor
  obs_u   = in->obs_u
  obs_v   = in->obs_v
  obs_lf  = in->obs_lf
  obs_lfx = cdft_p(reshape(obs_lf@tval,dimsizes(obs_lf)), \
                   reshape(obs_lf@nptxy,dimsizes(obs_lf))-2)
  copy_VarCoords(obs_lf,obs_lfx)
                   
  in = addfile("../d4PDF/figS2_tmpfile_d4pdf.nc","r")
  d4pdf_std = in->d4pdf_std
  d4pdf_cor = in->d4pdf_cor
  d4pdf_u   = in->d4pdf_u
  d4pdf_v   = in->d4pdf_v
  d4pdf_lf  = in->d4pdf_lf
  copy_VarCoords(d4pdf_cor,d4pdf_std)
  d4pdf_lfx = cdft_p(reshape(d4pdf_lf@tval,dimsizes(d4pdf_lf)), \
                     reshape(d4pdf_lf@nptxy,dimsizes(d4pdf_lf))-2)
  copy_VarCoords(d4pdf_lf,d4pdf_lfx)
  d4pdf_std_int = in->d4pdf_std_int
  d4pdf_cor_int = in->d4pdf_cor_int
  d4pdf_u_int   = in->d4pdf_u_int
  d4pdf_v_int   = in->d4pdf_v_int
  d4pdf_lf_int  = in->d4pdf_lf_int
  copy_VarCoords(d4pdf_cor,d4pdf_std_int)
  copy_VarCoords(d4pdf_cor,d4pdf_cor_int)
  copy_VarCoords(d4pdf_cor,d4pdf_u_int)
  copy_VarCoords(d4pdf_cor,d4pdf_v_int)
  copy_VarCoords(d4pdf_cor,d4pdf_lf_int)
  d4pdf_lfx_int = cdft_p(reshape(d4pdf_lf_int@tval,dimsizes(d4pdf_lf_int)), \
                         reshape(d4pdf_lf_int@nptxy,dimsizes(d4pdf_lf_int))-2)
  copy_VarCoords(d4pdf_lf_int,d4pdf_lfx_int)
  
  wks = gsn_open_wks("eps","figS2")
  plot1 = new(9,graphic)
  plot2 = new(3,graphic)
  plot3 = new(3,graphic)
  plot4 = new(3,graphic)

  res1                 = True
  res1@gsnDraw         = False
  res1@gsnFrame        = False
  res1@cnLinesOn       = False
  res1@cnFillOn        = True
  res1@cnInfoLabelOn   = False
  res1@lbLabelBarOn    = False
  res1@cnLineLabelsOn  = False
  ncl_default = read_colormap_file("ncl_default")
  res1@tmBorderThicknessF  = 1.25
  res1@tmXTOn = False
  res1@tmYROn = False
  res1@tmXBOn = False
  res1@tmYLOn = False
  res1@gsnLeftStringFontHeightF = 0.0293
  res1@cnLevelSelectionMode = "ManualLevels"

  res2                 = True
  res2@gsnDraw         = False
  res2@gsnFrame        = False
  res2@cnLinesOn       = False
  res2@cnLineLabelsOn  = False
  res2@cnFillOn        = False
  res2@cnInfoLabelOn   = False
  res2@cnFillColor     = "gray70"
  res2@cnFillScaleF    = 0.45
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLevels = (/0.0249,0.025,0.0251,0.9749,0.975,0.9751/)

  res4 = res2
  res4@cnLevels = r*(/-1.01,-1.,-.99,.99,1.,1.01/)

  res3  = True
  res3@gsnDraw                   = False
  res3@gsnFrame                  = False
  res3@vcGlyphStyle              = "FillArrow"
  res3@vcFillArrowWidthF         = 0.04
  res3@vcFillArrowHeadXF         = 0.4
  res3@vcFillArrowHeadYF         = 0.12
  res3@vcFillArrowHeadInteriorXF = 0.2
  res3@vcRefLengthF              = 0.04
  res3@vcRefMagnitudeF           = 1.
  obs_u = where(sqrt(obs_u*obs_u+obs_v*obs_v).lt.res3@vcRefMagnitudeF/5.,obs_u@_FillValue,obs_u)
  d4pdf_u = where(sqrt(d4pdf_u*d4pdf_u+d4pdf_v*d4pdf_v).lt.res3@vcRefMagnitudeF/5.,d4pdf_u@_FillValue,d4pdf_u)
  d4pdf_u_int = where(sqrt(d4pdf_u_int*d4pdf_u_int+d4pdf_v_int*d4pdf_v_int).lt.res3@vcRefMagnitudeF/5.,d4pdf_u_int@_FillValue,d4pdf_u_int)
  res3@vcRefAnnoOn             = False
  res3@vcRefAnnoFontHeightF    = 0.018
  res3@vcRefAnnoString1        = res3@vcRefMagnitudeF+" [m/s]"
  res3@vcRefAnnoString2On      = False
  res3@vcRefAnnoPerimOn        = False
  res3@vcRefAnnoOrthogonalPosF = 0.12
  res3@vcMinDistanceF          = 0.030
  res3@vcVectorDrawOrder       = "PostDraw"

  res1@mpMaxLatF       =  -8.
  res1@mpMinLatF       = -45.
  res1@mpMaxLonF       = 160.
  res1@mpMinLonF       = 110.
  res1@mpCenterLonF    = 135.
  res1@gsnAddCyclic    = False
  res1@cnFillPalette  := ncl_default(127:,:)
  res1@cnMinLevelValF  =  0.6
  res1@cnMaxLevelValF  =  4.5
  res1@cnLevelSpacingF =  0.3
  res1@mpLandFillColor = -1
  res1@mpOceanFillColor = "white"
  res1@mpFillDrawOrder  = "PostDraw"
  res1@gsnLeftString = "(a) AGCD"
  plot1(0) = gsn_csm_contour_map_ce(wks,obs_std,res1)
  res1@gsnLeftString = "(b) d4PDF HIST (Total)"
  plot1(1) = gsn_csm_contour_map_ce(wks,d4pdf_std,res1)
  res1@gsnLeftString = "(c) d4PDF HIST (Inter-member)"
  plot1(2) = gsn_csm_contour_map_ce(wks,d4pdf_std_int,res1)
  delete([/res1@mpLandFillColor,res1@mpOceanFillColor,res1@mpFillDrawOrder/])
  
  res1@mpMaxLatF       =  30.
  res1@mpMinLatF       = -30.
  res1@mpMaxLonF       = 360.
  res1@mpMinLonF       =   0.
  res1@mpCenterLonF    = 180.
  res1@gsnAddCyclic    = True
  res1@lbTitleString   = ""
  res1@cnFillPalette  := ncl_default(:,:)
  res1@cnMinLevelValF  = -0.6
  res1@cnMaxLevelValF  =  0.6
  res1@cnLevelSpacingF =  0.1
  res1@gsnLeftStringFontHeightF = 0.034
  res1@gsnLeftString = "(d) AGCD"
  plot1(3) = gsn_csm_contour_map_ce(wks,obs_cor,res1)
  plot4(0) = gsn_csm_contour(wks,obs_cor,res4)
  plot4(0) = ShadeGeLeContour(plot4(0),-r,r,6)
  overlay(plot1(3),plot4(0))
  res1@gsnLeftString = "(e) d4PDF HIST (Total)"
  plot1(4) = gsn_csm_contour_map_ce(wks,d4pdf_cor,res1)
  plot4(1) = gsn_csm_contour(wks,d4pdf_cor,res4)
  plot4(1) = ShadeGeLeContour(plot4(1),-r,r,6)
  overlay(plot1(4),plot4(1))
  res1@gsnLeftString = "(f) d4PDF HIST (Inter-member)"
  plot1(5) = gsn_csm_contour_map_ce(wks,d4pdf_cor_int,res1)
  plot4(2) = gsn_csm_contour(wks,d4pdf_cor_int,res4)
  plot4(2) = ShadeGeLeContour(plot4(2),-r,r,6)
  overlay(plot1(5),plot4(2))

  res1@mpMaxLatF       =   0.
  res1@mpMinLatF       = -45.
  res1@mpMaxLonF       = 180.
  res1@mpMinLonF       =  90.
  res1@mpCenterLonF    = 135.
  res1@gsnAddCyclic    = False
  res1@cnMinLevelValF  = -16.
  res1@cnMaxLevelValF  =  16.
  res1@cnLevelSpacingF =  2.
  res1@gsnLeftString = "(g) JRA-3Q"
  plot1(6) = gsn_csm_contour_map_ce(wks,obs_lf,res1)
  plot2(0) = gsn_csm_contour(wks,obs_lfx,res2)
  plot2(0) = ShadeGeLeContour(plot2(0),0.025,0.975,6)
  overlay(plot1(6),plot2(0))
  plot3(0) = gsn_csm_vector(wks,obs_u,obs_v,res3)
  overlay(plot1(6),plot3(0))
  res1@gsnLeftString = "(h) d4PDF HIST (Total)"
  plot1(7) = gsn_csm_contour_map_ce(wks,d4pdf_lf,res1)
  plot2(1) = gsn_csm_contour(wks,d4pdf_lfx,res2)
  plot2(1) = ShadeGeLeContour(plot2(1),0.025,0.975,6)
  overlay(plot1(7),plot2(1))
  plot3(1) = gsn_csm_vector(wks,d4pdf_u,d4pdf_v,res3)
  overlay(plot1(7),plot3(1))
  res1@gsnLeftString = "(i) d4PDF HIST (Inter-member)"
  plot1(8) = gsn_csm_contour_map_ce(wks,d4pdf_lf_int,res1)
  plot2(2) = gsn_csm_contour(wks,d4pdf_lfx_int,res2)
  plot2(2) = ShadeGeLeContour(plot2(2),0.025,0.975,6)
  overlay(plot1(8),plot2(2))
  res3@vcRefAnnoOn  = True
  plot3(2) = gsn_csm_vector(wks,d4pdf_u_int,d4pdf_v_int,res3)
  overlay(plot1(8),plot3(2))
  
  resp = True
  resp@gsnFrame = False
  resp@gsnPanelLabelBar  = True
  resp@lbBoxEndCapStyle  = "TriangleBothEnds"
  resp@lbLabelFontHeightF = 0.01
  resp@lbTitleOn          = True
  resp@lbTitleOffsetF     = -0.01
  resp@lbTitlePosition    = "Right"
  resp@lbTitleDirection   = "Across" 
  resp@lbTitleFontHeightF = 0.01
  resp@pmLabelBarWidthF  = 0.80
  resp@pmLabelBarHeightF = 0.05
  resp@pmLabelBarOrthogonalPosF = -0.01
  dy = 0.35
  y = (/0.95,0.67,0.39/)
  resp@gsnPanelTop = y(0)
  resp@gsnPanelBottom = resp@gsnPanelTop - dy
  resp@lbTitleString   = "[mm/day]"
  gsn_panel(wks,plot1(0:2),(/1,3/),resp)
  resp@gsnPanelTop = y(1)
  resp@gsnPanelBottom = resp@gsnPanelTop - dy
  resp@lbTitleString   = ""
  gsn_panel(wks,plot1(3:5),(/1,3/),resp)
  resp@gsnPanelTop = y(2)
  resp@gsnPanelBottom = resp@gsnPanelTop - dy
  resp@pmLabelBarOrthogonalPosF = 0.05
  resp@lbTitleString   = "[W/m~S~2~N~]"
  gsn_panel(wks,plot1(6:8),(/1,3/),resp)

  restx = True
  restx@txFontHeightF = 0.016
  restx@txJust = "CenterLeft"
  gsn_text_ndc(wks,"Interannual standard deviation of JF-mean rainfall",0.01,0.97,restx)
  gsn_text_ndc(wks,"Correlation b/w JF-mean SST and AUMR",0.01,0.57,restx)
  gsn_text_ndc(wks,"JF-mean 10-m wind and surface latent heat flux anomalies regressed onto AUMR",0.01,0.36,restx)
;  drawNDCGrid(wks)
  frame(wks)
  
end