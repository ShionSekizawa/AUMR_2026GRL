;;; File: derive_MSEbudget_advect.ncl
;;; Description: Calculate JF-mean horizontal and vertical advection terms of the MSE budget equation (Eq.(3)).
;;;
;;; Input: d4PDF data
;;; Output: MSEbudget_${experiment_type_name}_advect.nc
;;;
;;; Running:
;;; ncl derive_MSEbudget_advect.ncl I=${experiment_type_index}
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

;;; function for reading 4-byte binary data files of the d4PDF
function read_file(f:string,vname:string)
local xdef, ydef, zdef, msgval, lat, lon, lev, i, v, rec, var
begin
  if (vname.eq."U") then
    v = 0
  end if
  if (vname.eq."V") then
    v = 1
  end if
  if (vname.eq."W") then
    v = 2
  end if
  if (vname.eq."Z") then
    v = 3
  end if
  if (vname.eq."T") then
    v = 4
  end if
  if (vname.eq."Q") then
    v = 5
  end if
  
  xdef   = 288
  ydef   = 145
  zdef   = 12
  msgval = -9.99e33
  lat        = -90.+1.25*ispan(0,ydef-1,1)
  lat@units  = "degree_N"
  lon        = 0.+1.25*ispan(0,xdef-1,1)
  lon@units  = "degree_E"
  lev        = (/1000.,925.,850.,700.,600.,500.,400., \
                 300.,250.,200.,150.,100./)
  lev       := lev(::-1)
  lev@units  = "hPa"
  lat!0      = "lat"
  lat&lat    = lat
  lon!0      = "lon"
  lon&lon    = lon
  lev!0      = "lev"
  lev&lev    = lev
  
  setfileoption("bin", "ReadByteOrder", "BigEndian")
  var = new((/zdef,ydef,xdef/),float,msgval)
  do i = 0, zdef-1
    rec = v*24+(11-i)
    var(i,:,:) = fbindirread(f,rec,(/ydef,xdef/),"float")
  end do
  var!0 = "lev"
  var!1 = "lat"
  var!2 = "lon"
  var&lev = lev
  var&lat = lat
  var&lon = lon
  var@_FillValue = msgval
  return(var)
end

;;; function for performing vertical integration
function vint(z[*][*][*]:numeric,ps[*][*]:numeric)
local dim, lev, i, j, z_vint
begin
  dim  = dimsizes(z)
  lev  = (/1000.,925.,850.,700.,600.,500.,400., \
           300.,250.,200.,150.,100./)
  lev := lev(::-1)*100.
  z_vint = new(dim(1:2),float)
  do i = 0, dim(1)-1
    do j = 0, dim(2)-1
      z_vint(i,j) = ftcurvi(100*100.,ps(i,j),lev,z(:,i,j))
    end do
  end do
  z_vint = z_vint/9.80665
  return(z_vint)
end

;;; function for calculating vertical (pressure) derivative
function dp(h[*][*][*]:numeric)
local dim, lev, z, j, k
begin
  dim = dimsizes(h)
  lev  = (/1000.,925.,850.,700.,600.,500.,400., \
           300.,250.,200.,150.,100./)
  lev := lev(::-1)*100.
  dhdp = new(dim,float)
  do j = 0, dim(1)-1
    do k = 0, dim(2)-1
      dhdp(:,j,k) = ftcurvd(lev,h(:,j,k),lev)
    end do
  end do
  return(dhdp)
end

;;; function for calculating horizontal advection
function adv(z:numeric,u:numeric,v:numeric,LAT:integer,LON:integer)
local dx, dy
begin
  dx = z
  dy = z
  gradsf(z,dx,dy)
  return(-u(:,LAT,LON)*dx(:,LAT,LON)-v(:,LAT,LON)*dy(:,LAT,LON))
end

begin
  ;;; directory path referring to the d4PDF
  dir = "${d4PDF_dir_path}/Global/"
  
  ifn = "atm_avr_mon_1.25deg_"
  pfn = "sfc_avr_mon_"

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
        f = dir+expr+"_"+mmm+"/"+yyyymm+"/"+ifn+expr+"_"+mmm+"_"+yyyymm+".dr"
        p = dir+expr+"_"+mmm+"/"+yyyymm+"/"+pfn+expr+"_"+mmm+"_"+yyyymm+".grib"
        T = read_file(f,"T")
        if (yyyy.eq.y1.and.mn.eq.1) then
          lat = T&lat
          lon = T&lon
          LAT = ind(lat.ge.-50.and.lat.le.10)
          LON = ind(lon.ge.75.and.lon.le.195)
          NY  = dimsizes(LAT)
          NX  = dimsizes(LON) 
          shadv = new((/M,Y,NY,NX/),float) ;; horizontal advection of dry static energy (DSE)
          qhadv = new((/M,Y,NY,NX/),float) ;; horizontal advection of moisture
          svadv = new((/M,Y,NY,NX/),float) ;; vertical advection of DSE
          qvadv = new((/M,Y,NY,NX/),float) ;; vertical advection of moisture
        end if
        q = read_file(f,"Q")  ;; specific humidity
        z = read_file(f,"Z")  ;; geopotential
        s = cp*T + z  ;; DSE
        delete([/T,z/])
        u = read_file(f,"U")  ;; zonal wind
        v = read_file(f,"V")  ;; meridional wind
        poisson_grid_fill(s,True,1,1500,1.e-2,0.6,0)  ;; fill missing values with virtual values
        poisson_grid_fill(q,True,1,1500,1.e-2,0.6,0)  ;;      by solving the poisson equation
        poisson_grid_fill(u,True,1,1500,1.e-2,0.6,0)
        poisson_grid_fill(v,True,1,1500,1.e-2,0.6,0)
        in = addfile(p,"r")   ;; surface pressure
        ps = g2fsh(in->PS_GDS4_ISBL,(/145,288/))
        dum1 = vint(adv(s,u,v,LAT,LON),ps(LAT,LON))
        dum2 = vint(adv(q,u,v,LAT,LON),ps(LAT,LON))
        if (mn.eq.1) then
          shadv(m1-m,yyyy-y1,:,:) = dum1*d1
          qhadv(m1-m,yyyy-y1,:,:) = dum2*d1
        else
          shadv(m1-m,yyyy-y1,:,:) = (shadv(m1-m,yyyy-y1,:,:) + dum1*d2)/(d1+d2)
          qhadv(m1-m,yyyy-y1,:,:) = (qhadv(m1-m,yyyy-y1,:,:) + dum2*d2)/(d1+d2)
        end if
        delete([/u,v,dum1,dum2/])
        w = read_file(f,"W")  ;; vertical pressure velocity
        poisson_grid_fill(w,True,1,1500,1.e-2,0.6,0)
        dum1 = vint(-w(:,LAT,LON)*dp(s(:,LAT,LON)),ps(LAT,LON))
        dum2 = vint(-w(:,LAT,LON)*dp(q(:,LAT,LON)),ps(LAT,LON))
        if (mn.eq.1) then
          svadv(m1-m,yyyy-y1,:,:) = dum1*d1
          qvadv(m1-m,yyyy-y1,:,:) = dum2*d1
        else
          svadv(m1-m,yyyy-y1,:,:) = (svadv(m1-m,yyyy-y1,:,:) + dum1*d2)/(d1+d2)
          qvadv(m1-m,yyyy-y1,:,:) = (qvadv(m1-m,yyyy-y1,:,:) + dum2*d2)/(d1+d2)
        end if
        delete([/s,q,w,ps,dum1,dum2/])
      end do
    end do
  end do
  
  shadv!0 = "run"
  shadv!0 = "year"
  shadv!1 = "lat"
  shadv!2 = "lon"
  shadv&run  = m1+ispan(0,M-1,1)
  shadv&year = y1+ispan(0,Y-1,1)
  shadv&lat  = lat(LAT)
  shadv&lon  = lon(LON)
  copy_VarCoords(shadv,svadv)
  copy_VarCoords(shadv,qhadv)
  copy_VarCoords(shadv,qvadv)

  ;;; output
  o = "./MSEbudget_"+expr+"_advect.nc"  ;; output file path
  system("rm -f "+o)
  out = addfile(o,"c")
  out->shadv = shadv
  out->svadv = svadv
  out->qhadv = qhadv
  out->qvadv = qvadv

end