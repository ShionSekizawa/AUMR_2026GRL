;;; File: gmst.ncl
;;; Description: Calculate JF-mean global-mean surface air temperature (GMST) in CMIP6 simulations
;;;
;;; Input: CMIP6 data
;;; Output: gmst_${experiment_index}_mulens.nc
;;;
;;; Running:
;;; ncl gmst.ncl ex=${experiment_index}
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
  
  list1 = "FileList_tas_"+expr+"_mulens.txt"
  ofile = "gmst_"+expr+"_mulens.nc"
  Y = y2-y1+1

  gmst = new((/M/MM,MM,Y/),float)
  do i = 0, M-1
    model = str_split(models(i)," ")
    print(model(0)+": "+model(1)+" "+model(2))
    ifile = systemfunc("grep \/"+model(1)+"\/ "+list1+" | grep  \/"+model(2)+"\/")
    if (dimsizes(ifile).eq.1) then
      in = addfile(ifile,"r")
      time = in->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat = in->lat
      dum = (in->tas(t1:t2:12,:,:)*31. + \
             in->tas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    else
      in = addfiles(ifile,"r")
      ListSetType(in,"cat")
      time = in[:]->time
      yyyymm = cd_calendar(time,-1)
      t1 = ind(yyyymm.eq.y1*100+1)
      t2 = ind(yyyymm.eq.y2*100+1)
      lat = in[0]->lat
      dum = (in[:]->tas(t1:t2:12,:,:)*31. + \
             in[:]->tas(t1+1:t2+1:12,:,:)*28.25)/(31+28.25)
    end if
    dlat = lat(1:dimsizes(lat)-1)-lat(0:dimsizes(lat)-2)
    print("  dy:"+sprintf("%3.2f",min(dlat))+ \
             ", "+sprintf("%3.2f",max(dlat)))
    wgt = cos(lat*atan(1.)/45.)
    print("  ["+min(yyyymm)+"..."+max(yyyymm)+"]")
    gmst(i/MM,i%MM,:) = wgt_areaave(dum,wgt,1,0)
    if (i.eq.0) then
      gmst!0 = "model"
      gmst!1 = "ens"
      gmst!2 = "time"
      gmst&model = ispan(1,M/MM,1)
      gmst&ens   = ispan(1,MM,1)
      gmst&time  = time(t1:t2:12) 
    end if
    delete([/lat,dlat,wgt,in,time,yyyymm,t1,t2,dum,ifile/])
  end do

  system("rm -f "+ofile)
  o = addfile(ofile,"c")
  o->gmst = gmst
  
end