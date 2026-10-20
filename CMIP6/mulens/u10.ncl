;;; File: u10.ncl
;;; Description: Calculate JF-mean surface zonal wind over the Indian Ocean northwest off Australia (monsoon westerlies) in CMIP6 simulations
;;;
;;; Input: CMIP6 data
;;; Output: u10_${experiment_index}_mulens.nc
;;;
;;; Running:
;;; ncl u10.ncl ex=${experiment_index}
;;;
;;; ${experiment_index} is given as an integer "ex"
;;;   ex=0  Historical
;;;   ex=4  SSP5-8.5
;;;

begin

  MM = 5
  list0 = "ModelMemList_mulens.txt"
  M = toint(systemfunc("wc -l "+list0))
  models = asciiread(list0,M,"string")

  list9 = "FileList_sftlf_mulens.txt"
  file_sftlf = asciiread(list9,M,"string")

  if (ex.eq.0) then
    expr = "historical"
    y1 = 1951
    y2 = 2014
  end if
  if (ex.eq.2) then
    expr = "ssp245"
    y1 = 2015
    y2 = 2100
  end if
  if (ex.eq.4) then
    expr = "ssp585"
    y1 = 2015
    y2 = 2100
  end if
  
  list1 = "FileList_uas_"+expr+"_mulens.txt"
  ofile = "u10_"+expr+"_mulens.nc"
  Y = y2-y1+1

  u10 = new((/M/MM,MM,Y/),float)
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1)+" "+model(2))
    if (i%MM.eq.0) then
      in = addfile(file_sftlf(i),"r")
      wgt = in->sftlf
      wgt = 1.-wgt/max(wgt)
      print("  sftlf: ["+min(wgt)+"..."+max(wgt)+"]")
      lat = in->lat
      lon = in->lon
      LAT = ind(lat.lt.-20.or.lat.gt.-6)
      LON = ind(lon.lt.105.or.lon.gt.145)
      wgt(LAT,:) = 0.0
      wgt(:,LON) = 0.0
      if (i.eq.0) then
        system("rm -f u10_wgt.nc")
        ox = addfile("u10_wgt.nc","c")
      end if
      wgt!0 = "lat_"+sprinti("%0.2i",i)
      wgt!1 = "lon_"+sprinti("%0.2i",i)
      wgt&$wgt!0$ = lat
      wgt&$wgt!1$ = lon
      ox->$("wgt_"+sprinti("%0.2i",i))$ = wgt
      delete([/in/])
    end if
    ifile = systemfunc("grep \/"+model(1)+"\/ "+list1+" | grep \/"+model(2)+"\/")
    if (dimsizes(ifile).eq.1) then
      in = addfile(ifile,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in->uas(t1:t2:12,:,:)*31. + \
             in->uas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    else
      in = addfiles(ifile,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in[:]->uas(t1:t2:12,:,:)*31. + \
             in[:]->uas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    if (sum(abs(dimsizes(dum(0,:,:))-dimsizes(wgt))).ne.0) then
      print("  ! invalid coordinate")
      exit
    end if
    u10(i/MM,i%MM,:) = wgt_areaave2(dum,wgt,0)
    if (i.eq.0) then
      u10!0 = "model"
      u10!1 = "ens"
      u10!2 = "time"
      u10&model = ispan(1,M/MM,1)
      u10&ens   = ispan(1,MM,1)
      u10&time  = time(t1:t2:12) 
    end if
    print("  u10: ["+min(u10(i/MM,i%MM,:))+"..."+max(u10(i/MM,i%MM,:))+"]")
    delete([/in,time,yyyymm,t1,t2,dum,ifile/])
    if (i%MM.eq.MM-1) then
      delete([/lat,lon,LAT,LON,wgt/])
    end if
  end do

  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->u10 = u10
  
end