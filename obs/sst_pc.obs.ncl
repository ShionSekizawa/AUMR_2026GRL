;;; File: sst_pc.obs.ncl
;;; Description: Calculate correlation b/w AUMR and SST PC1 in the obervation
;;;
;;; Input: COBE-SST2 data and ausm_rain_obs.nc created by obs/make_ausm_rain.obs.ncl
;;; Output: no output file
;;;
;;; Running:
;;; ncl sst_pc.obs.ncl
;;;

begin

  y1 = 1951
  y2 = 2010

  ;; COBE-SST2
  sst = new((/y2-y1+1,180,360/),float)
  time = new(y2-y1+1,integer)
  do yyyy = y1, y2
    in = addfile("${COBE_SST2_dir_path}/mon/sst-glb."+yyyy+".dat.grb","r")
    a1 = in->WTMP_GDS0_DBSL(0,:,:)
    a2 = in->WTMP_GDS0_DBSL(1,:,:)
    d1 = days_in_month(yyyy,1)
    d2 = days_in_month(yyyy,2)
    sst(yyyy-y1,:,:) = (a1*d1+a2*d2)/(d1+d2)
    if (yyyy.eq.y1) then
      time(yyyy-y1) = d1
    else
      time(yyyy-y1) = time(yyyy-y1-1)+day_of_year(yyyy-1,12,31)
    end if
  end do
  lat = (/ in->g0_lat_1 /)
  lon = (/ in->g0_lon_2 /)
  lat@units = "degrees_N"
  lon@units = "degrees_E"
  ;; EOF (PCA) domain
  LAT = ind(lat.ge.-15.and.lat.le.15)
  LON = ind(lon.ge.0.and.lon.le.360)
  lat := lat(LAT)
  lon := lon(LON)
  sst := sst(:,LAT,LON)

  neof = 5
  wgt = sqrt(cos(lat*atan(1.)/45.))
  LON := ind(lon.ge.190.and.lon.le.270)
  sst = dim_rmvmean_n(sst,0)
  do j = 0, dimsizes(wgt)-1
    sst(:,j,:) = sst(:,j,:) * wgt(j)
  end do
  opt   = True
  opt@jopt = 0
  eof = eofunc_n(sst,neof,opt,0)
  pc  = eofunc_ts_n(sst,eof,False,0)
  pcvar = eof@pcvar
  do i = 0, neof-1
    eof(i,:,:) = eof(i,:,:)*sqrt(eof@eval(i))
    pc(i,:)  = pc(i,:)/sqrt(eof@eval(i))
    if (avg(eof(i,:,LON)).lt.0) then
      eof(i,:,:) = -eof(i,:,:)
      pc(i,:) = -pc(i,:)
    end if
    do j = 0, dimsizes(wgt)-1
      eof(:,j,:) = eof(:,j,:) / wgt(j)
    end do
  end do
  ;; set coordinates
  eof!0 = "eof_mode"
  eof!1 = "lat"
  eof!2 = "lon"
  eof&eof_mode = ispan(1,neof,1)
  eof&lat = lat
  eof&lon = lon
  pc!0 = "eof_mode"
  pc!1 = "time"
  pc&eof_mode = ispan(1,neof,1)
  pc&time = time
  pcvar!0 = "eof_mode"
  pcvar&eof_mode = ispan(1,neof,1)
  ; system("rm -f sst_pc.obs.nc")
  ; o = addfile("sst_pc.obs.nc","c")
  ; o->eof = eof
  ; o->pc  = pc
  ; o->pcvar = pcvar

  ;; AGCD-based AUMR
  in = addfile("ausm_rain_obs.nc","r")
  id = ind((in->year).ge.y1.and.(in->year).le.y2)
  ausm = in->ausm(id)

  ;; calculate correlation with its confidence interval
  r = escorc(ausm,pc(0,:))
  z = 0.5*log((1+r)/(1-r))
  zu = z + cdfnor_x(0.975,0.,1.)/sqrt(y2-y1+1-3)
  zl = z - cdfnor_x(0.975,0.,1.)/sqrt(y2-y1+1-3)
  ru = (exp(2*zu)-1)/(exp(2*zu)+1)
  rl = (exp(2*zl)-1)/(exp(2*zl)+1)
  print("r = "+r+" ( "+rl+" ... "+ru+" )")

end