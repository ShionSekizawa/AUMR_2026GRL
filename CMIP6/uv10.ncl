;;; File: uv10.ncl
;;; Description: Calculate JF-mean 10-m wind in CMIP6 simulations interpolated into 1x1-deg grid
;;;
;;; Input: CMIP6 data
;;; Output: uv10_${experiment_index}.nc
;;;
;;; Running:
;;; ncl uv10.ncl ex=${experiment_index}
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

  ;;; read list of files for surface wind
  list1 = "FileList_uas_"+expr+".txt"
  list2 = "FileList_vas_"+expr+".txt"
  ofile = "uv10_"+expr+".nc"
  Y = y2-y1+1

  NY = 120
  NX = 360
  lat = fspan(-59.5,59.5,NY)
  lon = fspan(0.5,359.5,NX)
  lat@units = "degrees_N"
  lon@units = "degrees_E"
  
  u10 = new((/M,Y,NY,NX/),float)  ;;; JF-mean u
  v10 = new((/M,Y,NY,NX/),float)  ;;; JF-mean v
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1))
    ;;; read surface wind data and perform JF-mean
    ifile1 = systemfunc("grep \/"+model(1)+"\/ "+list1)
    ifile2 = systemfunc("grep \/"+model(1)+"\/ "+list2)
    if (dimsizes(ifile1).eq.1) then
      in = addfile(ifile1,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat1i = in->lat
      lon1i = in->lon
      dum1 = (in->uas(t1:t2:12,:,:)*31. + \
              in->uas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
      in = addfile(ifile2,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat2i = in->lat
      lon2i = in->lon        
      dum2 = (in->vas(t1:t2:12,:,:)*31. + \
              in->vas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)        
    else
      in = addfiles(ifile1,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat1i = in[0]->lat
      lon1i = in[0]->lon
      dum1 = (in[:]->uas(t1:t2:12,:,:)*31. + \
              in[:]->uas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
      in = addfiles(ifile2,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat2i = in[0]->lat
      lon2i = in[0]->lon
      dum2 = (in[:]->vas(t1:t2:12,:,:)*31. + \
              in[:]->vas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    ;;; interpolation
    u10(i,:,:,:) = area_conserve_remap(lon1i,lat1i,dum1,lon,lat,False)
    v10(i,:,:,:) = area_conserve_remap(lon2i,lat2i,dum2,lon,lat,False)
    if (i.eq.0) then
      u10!0 = "model"
      u10!1 = "time"
      u10!2 = "lat"
      u10!3 = "lon"
      u10&model = ispan(1,M,1)
      u10&time = time(t1:t2:12)
      u10&lat = lat
      u10&lon = lon
      copy_VarCoords(u10,v10)
    end if
    print("  u10: ["+min(u10(i,:,:,:))+"..."+max(u10(i,:,:,:))+"]")
    print("  v10: ["+min(v10(i,:,:,:))+"..."+max(v10(i,:,:,:))+"]")
    delete([/in,time,yyyymm,t1,t2,dum1,dum2,ifile1,ifile2,lat1i,lon1i,lat2i,lon2i/])
  end do

  ;;; output
  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->u10 = u10
  o->v10 = v10
  
end