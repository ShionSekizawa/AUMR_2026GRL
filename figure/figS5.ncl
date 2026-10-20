;;; File: figS5.ncl
;;; Description: create Figure S5
;;;
;;; Input: sst_pc_historical.nc created by CMIP6/sst_pc.ncl
;;;        ausm_historical.nc created by CMIP6/ausm.ncl
;;;        uv10_historical.nc created by CMIP6/uv10.ncl
;;;        hfls_historical.nc created by CMIP6/hfls.ncl
;;;        CMIP6/u10_cov.txt created by CMIP6/cov.ncl
;;;        CMIP6/ModelMemList.txt
;;;
;;; Output: figS5.eps
;;;
;;; Running:
;;; ncl figS5.ncl
;;;

begin

  in = addfile("../CMIP6/sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  M = dimsizes(pc1(:,0))
  in = addfile("../CMIP6/ausm_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  ausm = dim_rmvmean(in->ausm(:,id))
  cov = asciiread("../CMIP6/u10_cov.txt",26,"float")
  M1 = ind(cov.le.0.12)
  M2 = ind(cov.gt.0.12)
  do i = 0, M-1
    ausm(i,:) = dtrend_msg(pc1(i,:),ausm(i,:),False,False)
    ausm(i,:) = dtrend_msg(yyyy(id),ausm(i,:),False,False)
    ausm(i,:) = ausm(i,:)/stddev(ausm(i,:))
  end do
  print(dimsizes(M1)+", "+dimsizes(M2))
  in = addfile("../CMIP6/uv10_historical.nc","r")
  lat = in->lat
  lon = in->lon
  NY = dimsizes(lat)
  NX = dimsizes(lon)
  u10 = new((/M,NY,NX/),float)
  v10 = new((/M,NY,NX/),float)
  do i = 0, M-1
    u10(i,:,:) = regCoef_n(ausm(i,:),dtrend_msg_n(yyyy(id),dtrend_msg_n(pc1(i,:),in->u10(i,id,:,:),False,False,0),False,False,0),0,0)
    v10(i,:,:) = regCoef_n(ausm(i,:),dtrend_msg_n(yyyy(id),dtrend_msg_n(pc1(i,:),in->v10(i,id,:,:),False,False,0),False,False,0),0,0)
  end do
  in = addfile("../CMIP6/hfls_historical.nc","r")
  lh  = new((/M,NY,NX/),float)
  lhx = new((/M,NY,NX/),float)
  do i = 0, M-1
    dum = regCoef_n(ausm(i,:),dtrend_msg_n(yyyy(id),dtrend_msg_n(pc1(i,:),in->hfls(i,id,:,:),False,False,0),False,False,0),0,0)
    lh(i,:,:)  = dum
    lhx(i,:,:) = cdft_p(reshape(dum@tval,dimsizes(dum)), \
                        reshape(dum@nptxy,dimsizes(dum))-2)
    delete(dum)
  end do

  m_tmp = asciiread("../CMIP6/ModelMemList.txt",M,"string")
  m_name = new(M,string)
  do i = 0, M-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do
  
  u10!0 = "model"
  u10!1 = "lat"
  u10!2 = "lon"
  u10&lat = lat
  u10&lon = lon
  copy_VarCoords(u10,v10)
  copy_VarCoords(u10,lh)
  copy_VarCoords(u10,lhx)

  wks = gsn_open_wks("eps","figS5")
  plot1 = new(30,graphic)
  plot2 = new(30,graphic)
  plot3 = new(30,graphic)

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
  res1@mpMaxLonF       = 180.
  res1@mpMinLonF       =  90.
  res1@mpCenterLonF    = 135.
  res1@gsnAddCyclic    = False
  res1@tmBorderThicknessF  = 1.25
  res1@tmXTOn = False
  res1@tmYROn = False
  res1@tmXBOn = False
  res1@tmYLOn = False
  res1@gsnLeftStringFontHeightF = 0.035
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMinLevelValF  = -24.
  res1@cnMaxLevelValF  =  24.
  res1@cnLevelSpacingF =  3.

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
  u10 = where(sqrt(u10*u10+v10*v10).lt.res3@vcRefMagnitudeF/5.,u10@_FillValue,u10)
  res3@vcRefAnnoOn             = False
  res3@vcRefAnnoFontHeightF    = 0.024
  res3@vcRefAnnoString1        = res3@vcRefMagnitudeF+" [m/s]"
  res3@vcRefAnnoString2On      = False
  res3@vcRefAnnoPerimOn        = False
  res3@vcRefAnnoOrthogonalPosF = 0.18
  res3@vcMinDistanceF          = 0.030
  res3@vcVectorDrawOrder       = "PostDraw"

  do i = 0, M-1
    res1@gsnLeftString = m_name(i)
    if (any(M1.eq.i)) then
      res1@gsnLeftString = res1@gsnLeftString+" *"
    end if
    plot1(i) = gsn_csm_contour_map_ce(wks,lh(i,:,:),res1)
    plot2(i) = gsn_csm_contour(wks,lhx(i,:,:),res2)
    plot2(i) = ShadeGeLeContour(plot2(i),0.025,0.975,6)
    overlay(plot1(i),plot2(i))
    if (i.eq.23) then
      res3@vcRefAnnoOn = True
    else
      res3@vcRefAnnoOn = False
    end if
    plot3(i) = gsn_csm_vector(wks,u10(i,:,:),v10(i,:,:),res3)
    overlay(plot1(i),plot3(i))
  end do

  resp = True
  resp@gsnFrame = False
  resp@gsnPanelLabelBar  = True
  resp@gsnPanelTop = 0.95
  resp@lbBoxEndCapStyle  = "TriangleBothEnds"
  resp@lbLabelFontHeightF   = 0.009
  resp@lbTitleOn            = True
  resp@lbTitleOffsetF       = -0.01
  resp@lbTitleString        = "[W/m~S~2~N~]"
  resp@lbTitlePosition      = "Right"
  resp@lbTitleDirection     = "Across" 
  resp@lbTitleFontHeightF   = 0.009
  resp@pmLabelBarWidthF  = 0.85
  resp@pmLabelBarHeightF = 0.03
  resp@pmLabelBarOrthogonalPosF = -0.01
  resp@gsnPanelXWhiteSpacePercent = 4
  resp@gsnPanelYWhiteSpacePercent = 4
  gsn_panel(wks,plot1(0:27),(/7,4/),resp)
  restx = True
  restx@txFontHeightF = 0.014
  restx@txJust = "CenterCenter"
  gsn_text_ndc(wks,"JF-mean 10-m wind and surface latent heat flux anomalies regressed onto AUMR in CMIP6 historical simulations",0.50,0.97,restx)
  frame(wks)
end