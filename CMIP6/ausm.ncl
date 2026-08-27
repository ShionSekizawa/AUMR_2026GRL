;;; File: ausm.ncl
;;; Description: Calculate JF-mean rainfall over northern Australia in CMIP6 simulations
;;;
;;; Input: CMIP6 data
;;; Output: ausm_${experiment_index}.nc
;;;
;;; Running:
;;; ncl ausm.ncl ex=${experiment_index}
;;;
;;; ${experiment_index} is given as an integer "ex"
;;;   ex=0  Historical
;;;   ex=2  SSP2-4.5
;;;   ex=4  SSP5-8.5
;;;

begin

  ;;; read list of models
  list0 = "ModelMemList.txt"
  M = toint(systemfunc("wc -l "+list0))
  models = asciiread(list0,M,"string")

  ;;; read list of files for land-area fraction
  list9 = "FileList_sftlf.txt"
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

  ;;; read list of files for precipitation
  list1 = "FileList_pr_"+expr+".txt"
  ofile = "ausm_"+expr+".nc"
  Y = y2-y1+1

  ausm = new((/M,Y/),float)  ;;; JF-mean rainfall over northern Australia
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1))
    ;;; read land-area fraction and create weight array for spatial average
    in = addfile(file_sftlf(i),"r")
    wgt = in->sftlf
    print("  sftlf: [0..."+max(wgt)+"]")
    lat = in->lat
    lon = in->lon
    LAT = ind(lat.lt.-20.or.lat.gt.-11)
    LON = ind(lon.lt.115.or.lon.gt.150)
    wgt(LAT,:) = 0.0
    wgt(:,LON) = 0.0
    wgt = wgt/max(wgt)
    ; if (i.eq.0) then
    ;   system("rm -f ausm_wgt.nc")
    ;   ox = addfile("ausm_wgt.nc","c")
    ; end if
    ; wgt!0 = "lat_"+sprinti("%0.2i",i)
    ; wgt!1 = "lon_"+sprinti("%0.2i",i)
    ; wgt&$wgt!0$ = lat
    ; wgt&$wgt!1$ = lon
    ; ox->$("wgt_"+sprinti("%0.2i",i))$ = wgt
    delete([/in,lat,lon,LAT,LON/])
    ;;; read precipitation data and perform JF-mean
    ifile = systemfunc("grep \/"+model(1)+"\/ "+list1)
    if (dimsizes(ifile).eq.1) then
      in = addfile(ifile,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in->pr(t1:t2:12,:,:)*31. + \
             in->pr(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    else
      in = addfiles(ifile,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in[:]->pr(t1:t2:12,:,:)*31. + \
             in[:]->pr(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    if (sum(abs(dimsizes(dum(0,:,:))-dimsizes(wgt))).ne.0) then
      print("  ! invalid coordinate")
      exit
    end if
    ;;; spatial average
    ausm(i,:) = wgt_areaave2(dum,wgt,0)*60*60*24.
    if (i.eq.0) then
      ausm!0 = "model"
      ausm!1 = "time"
      ausm&model = ispan(1,M,1)
      ausm&time = time(t1:t2:12) 
    end if
    print("  ausm: ["+min(ausm(i,:))+"..."+max(ausm(i,:))+"]")
    delete([/wgt,in,time,yyyymm,t1,t2,dum,ifile/])
  end do

  ;;; output
  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->ausm = ausm
  
end