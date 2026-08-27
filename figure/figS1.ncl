;;; File: figS1.ncl
;;; Description: create Figure S1
;;;
;;; Input: sst_clim_${experiment_type_name}.nc created by d4PDF/make_sst_clim.ncl
;;; Output: figS1.eps
;;;
;;; Running:
;;; ncl figS1.ncl
;;;

begin
  ;; get information on the d4PDF experiments
  data_list = str_split_csv(asciiread("../d4PDF/data_list",14,"string")," ",0)
  expr = data_list(:,0)
  m1 = toint(data_list(:,1))
  m2 = toint(data_list(:,2))
  y1 = toint(data_list(:,3))
  y2 = toint(data_list(:,4))

;;; d4PDF SST climatology
  in = addfile("../d4PDF/sst_clim_"+expr(0)+".nc","r")
  lat = in->lat
  lon = in->lon
  NY  = dimsizes(lat)
  NX  = dimsizes(lon)
  sst = new((/15,NY,NX/),float)
  do i = 0, 14
    sst(i,:,:) = in->sst
  end do
  delete(in)
  in = addfiles("../d4PDF/sst_clim_"+expr(2:7)+".nc","r")
  ListSetType(in,"join")
  sst(0,:,:) = dim_avg_n(in[:]->sst,0) - sst(0,:,:)
  in = addfiles("../d4PDF/sst_clim_"+expr(8:13)+".nc","r")
  ListSetType(in,"join")
  sst(1,:,:) = dim_avg_n(in[:]->sst,0) - sst(1,:,:)
  delete(in)
  do i = 0, 5
    in = addfile("../d4PDF/sst_clim_"+expr(i+2)+".nc","r")
    sst(i+3,:,:) = in->sst - sst(i+3,:,:) - sst(0,:,:)
    in = addfile("../d4PDF/sst_clim_"+expr(i+8)+".nc","r")
    sst(i+9,:,:) = in->sst - sst(i+9,:,:) - sst(1,:,:)
  end do

  sst!1 = "lat"
  sst!2 = "lon"
  sst&lat = lat
  sst&lon = lon

;;; open workstation
  wks = gsn_open_wks("eps","figS1")
  plot1 = new(15,graphic)
  res1                 = True
  res1@gsnDraw         = False
  res1@gsnFrame        = False
  res1@cnLinesOn       = False
  res1@cnFillOn        = True
  res1@cnInfoLabelOn   = False
  res1@lbLabelBarOn    = False
  res1@cnLineLabelsOn  = False
  res1@mpMaxLatF       =  30.
  res1@mpMinLatF       = -30.
  res1@mpMaxLonF       = 360.
  res1@mpMinLonF       =   0.
  res1@mpCenterLonF    = 180.
  res1@gsnAddCyclic    = True
  res1@tmBorderThicknessF  = 1.2
  res1@tmYLMajorThicknessF = 1.2
  res1@tmXBMajorThicknessF = 1.2
  res1@tmXBMajorLengthF = 0.004
  res1@tmYLMajorLengthF = 0.004
  res1@tmXBMajorOutwardLengthF = 0
  res1@tmYLMajorOutwardLengthF = 0
  res1@tmYLMinorOutwardLengthF = 0
  res1@tmXBMinorOutwardLengthF = 0
  res1@tmYLLabelFontHeightF = 0.014
  res1@tmXBLabelFontHeightF = 0.014
  res1@cnLevelSelectionMode = "ExplicitLevels"
  res1@gsnLeftStringFontHeightF = 0.024
  res1@lbBoxEndCapStyle    = "TriangleBothEnds"
  res1@lbLabelFontHeightF  = 0.018
  res1@lbTitleOn            = True
  res1@lbTitleOffsetF       = -0.01
  res1@lbTitleString        = "[K]"
  res1@lbTitlePosition      = "Right"
  res1@lbTitleDirection     = "Across" 
  res1@lbTitleFontHeightF   = 0.02
  res1@pmLabelBarWidthF    = 1.6
  res1@pmLabelBarHeightF   = 0.06
  res1@pmLabelBarOrthogonalPosF = 0.4
  cmp = read_colormap_file("ncl_default")
  res1@cnFillPalette := cmp(127:,:)
  res1@cnLevels := fspan(0.4,4.,19)
  res1@gsnLeftString = "(a) 2K (All) - HIST"
  plot1(0) = gsn_csm_contour_map_ce(wks,sst(0,:,:),res1)
  res1@lbLabelBarOn  = True
  res1@gsnLeftString = "(b) 4K (All) - HIST"
  plot1(1) = gsn_csm_contour_map_ce(wks,sst(1,:,:),res1)
  res1@cnFillPalette := cmp(:,:)
  res1@cnLevels := fspan(-1.8,1.8,19)
  LeftString = (/"(c) 2K (CC) - 2K (All)", \
                 "(d) 2K (GF) - 2K (All)", \
                 "(e) 2K (HA) - 2K (All)", \
                 "(f) 2K (MI) - 2K (All)", \
                 "(g) 2K (MP) - 2K (All)", \
                 "(h) 2K (MR) - 2K (All)", \
                 "(i) 4K (CC) - 4K (All)", \
                 "(j) 4K (GF) - 4K (All)", \
                 "(k) 4K (HA) - 4K (All)", \
                 "(l) 4K (MI) - 4K (All)", \
                 "(m) 4K (MP) - 4K (All)", \
                 "(n) 4K (MR) - 4K (All)" /)
  do i = 3, 14
    if (i.eq.13) then
      res1@lbLabelBarOn = True
    else
      res1@lbLabelBarOn = False
    end if
    res1@gsnLeftString = LeftString(i-3)
    plot1(i) = gsn_csm_contour_map_ce(wks,sst(i,:,:),res1)
  end do
  
  resp = True
  resp@gsnFrame = False
  resp@gsnPanelTop = 1.0
  resp@gsnPanelBottom = 0.3
  gsn_panel(wks,plot1(:2),(/1,3/),resp)
  resp@gsnPanelTop = 0.75
  resp@gsnPanelBottom = 0.0
  resp@gsnPanelYWhiteSpacePercent = 10
  gsn_panel(wks,plot1(3:),(/4,3/),resp)

;  drawNDCGrid(wks)
  frame(wks)
  
end