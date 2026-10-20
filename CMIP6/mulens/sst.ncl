;;; File: sst.ncl
;;; Description: Calculate JF-mean SST in CMIP6 simulations interpolated into 1x1-deg grid
;;;
;;; Input: CMIP6 data
;;; Output: sst_${experiment_index}_mulens.nc
;;;
;;; Running:
;;; ncl sst.ncl ex=${experiment_index}
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
  
  list1 = "FileList_ts_"+expr+"_mulens.txt"
  ofile = "sst_"+expr+"_mulens.nc"
  Y = y2-y1+1

  NY = 120
  NX = 360
  lat = fspan(-59.5,59.5,NY)
  lon = fspan(0.5,359.5,NX)
  lat@units = "degrees_N"
  lon@units = "degrees_E"
  
  sst = new((/M/MM,MM,Y,NY,NX/),float)
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1)+" "+model(2))
    if (i%MM.eq.0) then
      in = addfile(file_sftlf(i),"r")
      sftlf = in->sftlf
      lati = in->lat
      loni = in->lon
      delete(in)
    end if
    ifile = systemfunc("grep \/"+model(1)+"\/ "+list1+" | grep \/"+model(2)+"\/")
    if (dimsizes(ifile).eq.1) then
      in = addfile(ifile,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in->ts(t1:t2:12,:,:)*31. + \
             in->ts(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    else
      in = addfiles(ifile,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in[:]->ts(t1:t2:12,:,:)*31. + \
             in[:]->ts(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    if (sum(abs(dimsizes(dum(0,:,:))-dimsizes(sftlf))).ne.0) then
      print("  ! invalid coordinate")
      exit
    end if
    dumx = area_conserve_remap(loni,lati,sftlf,lon,lat,False)
    sst(i/MM,i%MM,:,:,:) = area_conserve_remap(loni,lati,dum,lon,lat,False)
    sst(i/MM,i%MM,:,:,:) = mask(sst(i/MM,i%MM,:,:,:),dumx.lt.0.5,True)
    if (i.eq.0) then
      sst!0 = "model"
      sst!1 = "ens"
      sst!2 = "time"
      sst!3 = "lat"
      sst!4 = "lon"
      sst&model = ispan(1,M/MM,1)
      sst&ens  = ispan(1,MM,1)
      sst&time = time(t1:t2:12)
      sst&lat = lat
      sst&lon = lon
    end if
    print("  sst: ["+min(sst(i/MM,i%MM,:,:,:))+"..."+max(sst(i/MM,i%MM,:,:,:))+"]")
    delete([/in,time,yyyymm,t1,t2,dum,ifile,dumx/])
    if (i%MM.eq.MM-1) then
      delete([/lati,loni,sftlf/])
    end if
  end do

  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->sst = sst
  
end