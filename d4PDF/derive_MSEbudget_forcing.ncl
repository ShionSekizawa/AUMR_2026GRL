;;; File: derive_MSEbudget_forcing.ncl
;;; Description: Calculate JF-mean surface heat flux and radiative heating terms of the MSE budget equation (Eq.(3)).
;;;
;;; Input: d4PDF data
;;; Output: MSEbudget_${experiment_type_name}_forcing.nc
;;;
;;; Running:
;;; ncl derive_MSEbudget_forcing.ncl I=${experiment_type_index}
;;;
;;; ${experiment_type_index} is given as an integer "I"
;;;   I= 0  HPB
;;;   I= 1  HPB_NAT
;;;   I= 2  HFB_2K_CC
;;;   I= 3  HFB_2K_GF
;;;   I= 4  HFB_2K_HA
;;;   I= 5  HFB_2K_MI
;;;   I= 6  HFB_2K_MP
;;;   I= 7  HFB_2K_MR
;;;   I= 8  HFB_4K_CC
;;;   I= 9  HFB_4K_GF
;;;   I=10  HFB_4K_HA
;;;   I=11  HFB_4K_MI
;;;   I=12  HFB_4K_MP
;;;   I=13  HFB_4K_MR
;;;

begin
  ;;; directory path referring to the d4PDF
  dir = "${d4PDF_dir_path}/Global/"
  
  ifn = "sfc_avr_mon_"

  ;;; select experiment type    
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(I,0)
  m1 = toint(data_list(I,1))
  m2 = toint(data_list(I,2))
  y1 = toint(data_list(I,3))
  y2 = toint(data_list(I,4))
  print("get dataset info from data_list: expr="+expr+" run=["+m1+".."+m2+"] year=["+y1+".."+y2+"]")

  ;;; constants
  cp = 1004.6  ;; J K^-1 kg^-1
  L  = 2.507e6 ;; J kg^-1

  ;;; read d4PDF data
  M = m2-m1+1
  Y = y2-y1+1
  do m = m1, m2
    mmm = "m"+sprinti("%0.3i",m)
    do yyyy = y1, y2
      d1 = days_in_month(yyyy,1)
      d2 = days_in_month(yyyy,2)
      do mn = 1, 2
        yyyymm = yyyy+sprinti("%0.2i",mn)
        print("  "+mmm+" "+yyyymm)
        f = dir+expr+"_"+mmm+"/"+yyyymm+"/"+ifn+expr+"_"+mmm+"_"+yyyymm+".grib"
        if (m.eq.m1.and.yyyy.eq.y1.and.mn.eq.1) then
          lat        = -90.+1.25*ispan(0,145-1,1)
          lat@units  = "degree_N"
          lon        = 0.+1.25*ispan(0,288-1,1)
          lon@units  = "degree_E"
          lat!0      = "lat"
          lat&lat    = lat
          lon!0      = "lon"
          lon&lon    = lon
          LAT = ind(lat.ge.-50.and.lat.le.10)
          LON = ind(lon.ge.75.and.lon.le.195)
          NY  = dimsizes(LAT)
          NX  = dimsizes(LON)
          cwv = new((/M,Y,NY,NX/),float) ;; Column-intagrated water vapor
          lwr = new((/M,Y,NY,NX/),float) ;; long-wave radiation heating
          swr = new((/M,Y,NY,NX/),float) ;; short-wave radiation heating
          lhf = new((/M,Y,NY,NX/),float) ;; surface latent heat flux
          shf = new((/M,Y,NY,NX/),float) ;; surface sensible heat flux
        end if
        in = addfile(f,"r")
        dum1 = g2fsh(in->PWATER_GDS4_ISBL,(/145,288/))
        ut = g2fsh(in->ULWT_GDS4_ISBL,(/145,288/)) ;; upward long-wave @TOA
        us = g2fsh(in->ULWB_GDS4_ISBL,(/145,288/)) ;; upward long-wave @Surface
        ds = g2fsh(in->DLWB_GDS4_ISBL,(/145,288/)) ;; downward long-wave @Surface
        dum2 = us - ds - ut  ;; net long-wave heating within the atmosphere
        delete([/ut,us,ds/])
        ut = g2fsh(in->USWT_GDS4_ISBL,(/145,288/)) ;; upward short-wave @TOA
        dt = g2fsh(in->DSWT_GDS4_ISBL,(/145,288/)) ;; downward short-wave @TOA
        us = g2fsh(in->USWB_GDS4_ISBL,(/145,288/)) ;; upward short-wave @Surface
        ds = g2fsh(in->DSWB_GDS4_ISBL,(/145,288/)) ;; downward short-wave @Surface
        dum3 = dt + us - ds - ut ;; net short-wave heating within the atmosphere
        delete([/ut,dt,us,ds/])
        dum4 = g2fsh(in->FLLH_GDS4_ISBL,(/145,288/)) ;; surface latent heat flux
        dum5 = g2fsh(in->FLSH_GDS4_ISBL,(/145,288/)) ;; surface sensible heat flux
        if (mn.eq.1) then
          cwv(m-m1,yyyy-y1,:,:) = dum1(LAT,LON)*d1
          lwr(m-m1,yyyy-y1,:,:) = dum2(LAT,LON)*d1
          swr(m-m1,yyyy-y1,:,:) = dum3(LAT,LON)*d1
          lhf(m-m1,yyyy-y1,:,:) = dum4(LAT,LON)*d1
          shf(m-m1,yyyy-y1,:,:) = dum5(LAT,LON)*d1
        else
          cwv(m-m1,yyyy-y1,:,:) = (cwv(m-m1,yyyy-y1,:,:) + dum1(LAT,LON)*d2)/(d1+d2)
          lwr(m-m1,yyyy-y1,:,:) = (lwr(m-m1,yyyy-y1,:,:) + dum2(LAT,LON)*d2)/(d1+d2)
          swr(m-m1,yyyy-y1,:,:) = (swr(m-m1,yyyy-y1,:,:) + dum3(LAT,LON)*d2)/(d1+d2)
          lhf(m-m1,yyyy-y1,:,:) = (lhf(m-m1,yyyy-y1,:,:) + dum4(LAT,LON)*d2)/(d1+d2)
          shf(m-m1,yyyy-y1,:,:) = (shf(m-m1,yyyy-y1,:,:) + dum5(LAT,LON)*d2)/(d1+d2)
        end if
        delete([/dum1,dum2,dum3,dum4,dum5/])
      end do
    end do
  end do
  
  cwv!0 = "run"
  cwv!1 = "year"
  cwv!2 = "lat"
  cwv!3 = "lon"
  cwv&run  = m1+ispan(0,M-1,1)
  cwv&year = y1+ispan(0,Y-1,1)
  cwv&lat  = lat(LAT)
  cwv&lon  = lon(LON)
  copy_VarCoords(cwv,lwr)
  copy_VarCoords(cwv,swr)
  copy_VarCoords(cwv,lhf)
  copy_VarCoords(cwv,shf)
  
  ;;; output
  o = "./MSEbudget_"+expr+"_forcing.nc" ;; output file path
  system("rm -f "+o)
  out = addfile(o,"c")
  out->cwv = cwv
  out->lwr = lwr
  out->swr = swr
  out->lhf = lhf
  out->shf = shf

end