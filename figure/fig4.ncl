;;; File: fig4.ncl
;;; Description: create Figure 4
;;;
;;; Input: lhf_area_${experiment_type_name}.nc created by d4PDF/make_lhf_area.ncl
;;;        sst_hist_ens_mean.nc created by d4PDF/make_sst_hist.ncl
;;;        sst_clim_${experiment_type_name}.nc created by d4PDF/make_sst_clim.ncl
;;;        sst_pc_historical.nc, sst_pc_ssp245.nc, sst_pc_ssp585.nc created by CMIP6/sst_pc.ncl
;;;        ausm_historical.nc created by CMIP6/ausm.ncl
;;;        sst_historical.nc, sst_ssp245.nc, sst_ssp585.nc created by CMIP6/sst.ncl
;;;        u10_historical.nc, u10_ssp245.nc, u10_ssp585.nc created by CMIP6/u10.ncl
;;;        gmst_historical.nc, gmst_ssp245.nc, gmst_ssp585.nc created by CMIP6/gmst.ncl
;;;        gwlev.nc created by CMIP6/gwlev.ncl
;;; Output: fig4.eps
;;;
;;; Running:
;;; ncl fig4.ncl
;;;

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
  in2 = addfile("../d4PDF/sst_hist_ens_mean.nc","r")
  lat = in2->lat
  lon = in2->lon
  cor_hist = escorc_n(dim_avg_n(in1->u10,0),in2->sst,0,0)
  delete([/in1,in2/])
  
  in1 = addfiles("../d4PDF/lhf_area_"+expr(2:7)+".nc","r")
  in2 = addfiles("../d4PDF/sst_clim_"+expr(2:7)+".nc","r")
  ListSetType(in1,"join")
  ListSetType(in2,"join")
  sst_2K = dim_rmvmean_n(dim_avg_n(in1[:]->u10(:,:,:),(/1,2/)),0)
  u10_2K = dim_rmvmean_n(in2[:]->sst(:,:,:),0)
  in1 = addfiles("../d4PDF/lhf_area_"+expr(8:13)+".nc","r")
  in2 = addfiles("../d4PDF/sst_clim_"+expr(8:13)+".nc","r")
  ListSetType(in1,"join")
  ListSetType(in2,"join")
  sst_4K = dim_rmvmean_n(dim_avg_n(in1[:]->u10(:,:,:),(/1,2/)),0)
  u10_4K = dim_rmvmean_n(in2[:]->sst(:,:,:),0)
  delete([/in1,in2/])
  cor_2K4K = escorc_n(array_append_record(u10_2K,u10_4K,0), \
                      array_append_record(sst_2K,sst_4K,0),0,0)

  cor_hist!0 = "lat"
  cor_hist!1 = "lon"
  cor_hist&lat = lat
  cor_hist&lon = lon
  copy_VarCoords(cor_hist,cor_2K4K)
  delete([/lat,lon/])
  
;;; CMIP6 data
  in = addfile("../CMIP6/sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  M = dimsizes(pc1(:,0))
  in = addfile("../CMIP6/ausm_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  cor = escorc(dim_rmvmean(in->ausm(:,id)),pc1)
  M1 = ind(abs(cor).ge.0.4)
  M2 = ind(abs(cor).lt.0.4)
  in = addfile("../CMIP6/sst_historical.nc","r")
  lat = in->lat
  lon = in->lon
  sst_c6_hist = dim_avg_n(in->sst(:,id,:,:),1)
  in = addfile("../CMIP6/gmst_historical.nc","r")
  gmst = in->gmst(:,id)
  gmst_hist = dim_avg(gmst)
  in = addfile("../CMIP6/u10_historical.nc","r")
  u10_c6_hist = dim_avg(in->u10(:,id))
  delete([/yyyy,id,gmst/])

  in = addfile("../CMIP6/sst_ssp245.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2021.and.yyyy.le.2080)
  sst_c6_245 = dim_avg_n(in->sst(:,id,:,:),1)
  in = addfile("../CMIP6/gmst_ssp245.nc","r")
  gmst = in->gmst(:,id)
  gmst_245 = dim_avg(gmst)
  in = addfile("../CMIP6/u10_ssp245.nc","r")
  u10_c6_245 = dim_avg(in->u10(:,id))
  delete([/yyyy,id,gmst/])

  in = addfile("../CMIP6/sst_ssp585.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.2041.and.yyyy.le.2100)
  sst_c6_585 = dim_avg_n(in->sst(:,id,:,:),1)
  in = addfile("../CMIP6/gmst_ssp585.nc","r")
  gmst = in->gmst(:,id)
  gmst_585 = dim_avg(gmst)
  in = addfile("../CMIP6/u10_ssp585.nc","r")
  u10_c6_585 = dim_avg(in->u10(:,id))
  delete([/yyyy,id,gmst/])

  in = addfile("../CMIP6/gwlev.nc","r")
  gwlev_hist = in->gwlev_historical
  gwlev_245 = in->gwlev_ssp245
  gwlev_585 = in->gwlev_ssp585
  sst_c6_245 = sst_c6_245-sst_c6_hist
  sst_c6_585 = sst_c6_585-sst_c6_hist
  do i = 0, M-1
    sst_c6_245(i,:,:) = sst_c6_245(i,:,:)*(2.-gwlev_hist(i))/(gwlev_245(i)-gwlev_hist(i))
    sst_c6_585(i,:,:) = sst_c6_585(i,:,:)*(4.-gwlev_hist(i))/(gwlev_585(i)-gwlev_hist(i))
  end do
  sst_c6_245 = dim_rmvmean_n(sst_c6_245,0)
  sst_c6_585 = dim_rmvmean_n(sst_c6_585,0)
  u10_c6_245 = (u10_c6_245-u10_c6_hist)*(2.-gwlev_hist)/(gwlev_245-gwlev_hist)
  u10_c6_585 = (u10_c6_585-u10_c6_hist)*(4.-gwlev_hist)/(gwlev_585-gwlev_hist)
  u10_c6_245 = dim_rmvmean_n(u10_c6_245,0)
  u10_c6_585 = dim_rmvmean_n(u10_c6_585,0)
  
  cor_c6 = escorc_n(array_append_record(u10_c6_245(M2),u10_c6_585(M2),0), \
                    array_append_record(sst_c6_245(M2,:,:),sst_c6_585(M2,:,:),0),0,0)
  cor_c6!0 = "lat"
  cor_c6!1 = "lon"
  cor_c6&lat = lat
  cor_c6&lon = lon
  delete([/lat,lon/])

  ;;; correlations at the 95% confidence level
  n1 = 6*2-2
  n2 = y2(0)-y1(0)+1-2
  n3 = dimsizes(M2)*2-2
  t = cdft_t(0.975,n1)
  r1 = t/sqrt(n1+t^2)
  t = cdft_t(0.975,n2)
  r2 = t/sqrt(n2+t^2)
  t = cdft_t(0.975,n3)
  r3 = t/sqrt(n3+t^2)
  print(r1+" "+r2+" "+r3)

;;; open workstation
  wks = gsn_open_wks("eps","fig4")
  plot1 = new(3,graphic)
  plot2 = new(3,graphic)
  box = new(3,graphic)
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
  res1@tmYLLabelFontHeightF = 0.008
  res1@tmXBLabelFontHeightF = 0.008
  res1@cnLevelSelectionMode = "ExplicitLevels"
  res1@cnLevels = fspan(-0.9,0.9,19)
  res1@gsnLeftStringFontHeightF = 0.015

  res2                 = True
  res2@gsnDraw         = False
  res2@gsnFrame        = False
  res2@cnLinesOn       = False
  res2@cnLineLabelsOn  = False
  res2@cnFillOn        = False
  res2@cnInfoLabelOn   = False
  res2@cnFillColor     = "gray70"
  res2@cnFillScaleF    = 0.6
  res2@cnLevelSelectionMode = "ExplicitLevels"

  res1@gsnLeftString = "(a) ~F33~D~F~SST~N~-induced uncertainty in the d4PDF 2K/4K "
  plot1(0) = gsn_csm_contour_map_ce(wks,cor_2K4K,res1)
  res2@cnLevels = (/-1.01,-1.,-0.99,0.99,1.,1.01/)*r1
  plot2(0) = gsn_csm_contour(wks,cor_2K4K,res2)
  plot2(0) = ShadeGeLeContour(plot2(0),-r1,r1,6)
  overlay(plot1(0),plot2(0))

  res1@gsnLeftString = "(b) Ensemble-mean in the d4PDF HIST"
  plot1(1) = gsn_csm_contour_map_ce(wks,cor_hist,res1)
  res2@cnLevels = (/-1.01,-1.,-0.99,0.99,1.,1.01/)*r2
  plot2(1) = gsn_csm_contour(wks,cor_hist,res2)
  plot2(1) = ShadeGeLeContour(plot2(1),-r2,r2,6)
  overlay(plot1(1),plot2(1))

  res1@gsnLeftString = "(c) Inter-model uncertainty in CMIP6 (19 models)"
  plot1(2) = gsn_csm_contour_map_ce(wks,cor_c6,res1)
  res2@cnLevels = (/-1.01,-1.,-0.99,0.99,1.,1.01/)*r3
  plot2(2) = gsn_csm_contour(wks,cor_c6,res2)
  plot2(2) = ShadeGeLeContour(plot2(2),-r3,r3,6)
  overlay(plot1(2),plot2(2))

  resbx = True
  resbx@gsLineColor = "magenta"
  resbx@gsLineThicknessF = 2
  resbx@tfPolyDrawOrder = "PostDraw"
  s1 = -20
  s2 =  -6
  s3 = 105
  s4 = 150
  do i = 0, 2
    box(i) = gsn_add_polyline(wks,plot1(i),(/s3,s4,s4,s3,s3/),(/s1,s1,s2,s2,s1/),resbx)
  end do

  resp = True
  resp@gsnFrame = False
  resp@gsnPanelYWhiteSpacePercent = 5
  resp@gsnPanelLabelBar    = True
  resp@lbBoxEndCapStyle    = "TriangleBothEnds"
  resp@lbLabelFontHeightF  = 0.011
  resp@lbTitleOn           = False
  resp@pmLabelBarWidthF    = 0.85
  resp@pmLabelBarHeightF   = 0.05
  resp@pmLabelBarOrthogonalPosF = -0.012
  gsn_panel(wks,plot1,(/3,1/),resp)

  restx = True
  restx@txFontHeightF = 0.02
  restx@txJust = "BottomLeft"
  gsn_text_ndc(wks,"Correlation b/w background SST & monsoon westerlies",0.04,0.89,restx)

;  drawNDCGrid(wks)
  frame(wks)
  
end