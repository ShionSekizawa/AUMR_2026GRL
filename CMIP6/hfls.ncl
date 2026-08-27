;;; File: hfls.ncl
;;; Description: Calculate JF-mean surface latent heat flux in CMIP6 simulations interpolated into 1x1-deg grid
;;;
;;; Input: CMIP6 data
;;; Output: hfls_${experiment_index}.nc
;;;
;;; Running:
;;; ncl hfls.ncl ex=${experiment_index}
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

  ;;; read list of files for surface latent heat flux
  list1 = "FileList_hfls_"+expr+".txt"
  ofile = "hfls_"+expr+".nc"
  Y = y2-y1+1

  NY = 120
  NX = 360
  lat = fspan(-59.5,59.5,NY)
  lon = fspan(0.5,359.5,NX)
  lat@units = "degrees_N"
  lon@units = "degrees_E"
  
  hfls = new((/M,Y,NY,NX/),float)  ;;; JF-mean surface latent heat flux
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1))
    in = addfile(file_sftlf(i),"r")
    sftlf = in->sftlf
    lati = in->lat
    loni = in->lon
    delete(in)
    ;;; read surface latent heat flux data and perform JF-mean
    ifile = systemfunc("grep \/"+model(1)+"\/ "+list1)
    if (dimsizes(ifile).eq.1) then
      in = addfile(ifile,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in->hfls(t1:t2:12,:,:)*31. + \
             in->hfls(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    else
      in = addfiles(ifile,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      dum = (in[:]->hfls(t1:t2:12,:,:)*31. + \
             in[:]->hfls(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    if (sum(abs(dimsizes(dum(0,:,:))-dimsizes(sftlf))).ne.0) then
      print("  ! invalid coordinate")
      exit
    end if
    ;;; interpolation
    dumx = area_conserve_remap(loni,lati,sftlf,lon,lat,False)
    hfls(i,:,:,:) = area_conserve_remap(loni,lati,dum,lon,lat,False)
    if (i.eq.0) then
      hfls!0 = "model"
      hfls!1 = "time"
      hfls!2 = "lat"
      hfls!3 = "lon"
      hfls&model = ispan(1,M,1)
      hfls&time = time(t1:t2:12)
      hfls&lat = lat
      hfls&lon = lon
    end if
    print("  hfls: ["+min(hfls(i,:,:,:))+"..."+max(hfls(i,:,:,:))+"]")
    delete([/in,time,yyyymm,t1,t2,dum,ifile,lati,loni,sftlf,dumx/])
  end do

  ;;; output
  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->hfls = hfls
  
end