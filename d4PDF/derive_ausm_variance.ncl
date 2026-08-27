;;; File: derive_ausm_variance.ncl
;;;
;;; Input: ausm_rain_${experiment_type_name}.nc created by d4PDF/make_ausm_rain.ncl
;;; Output: ausm_variance.nc
;;;
;;; Running:
;;; ncl derive_ausm_variance.ncl
;;;

;;; function for calculating AUMR mean and variance based on Eqs. (1) and (2)
function calc_variance(x[*][*]:numeric,clim_opt[1]:integer)
local dim, n, N, x1, x2, x0, mean, vint, vem, vsst, vtot, a, l_mean, l_vint, l_vsst, l_vtot, u_mean, u_vint, u_vsst, u_vtot,ret
begin
  dim = dimsizes(x)
  n   = dim(0)
  N   = dim(1)
  mean = avg(x)
  if (clim_opt.eq.6) then
    ;; Climatology is derived by taking the average for members with each SST warming pattern
    do i = 0, clim_opt-1
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:) = \
      x(i*n/clim_opt:(i+1)*n/clim_opt-1,:) - avg(x(i*n/clim_opt:(i+1)*n/clim_opt-1,:))
    end do
  else
    ;; Climatology is derived by taking the average for all members
    x(:,:) = x(:,:) - avg(x(:,:))
  end if
  x1  = dim_rmvmean_n(x,0)
  vint = sum(x1*x1)/(N*(n-1.)) ;; Internal component
  x2  = dim_avg_n(x,0)
  vem  = sum(x2*x2)/(N-1.)
  vsst = vem - vint/n          ;; SST-forced component
  vtot = vsst + vint           ;; Total = SST-forced + Internal
  ;; lower limit of 95% confidence level
  a = 0.05
  l_mean = mean - avg(x*x)*cdft_t(1-a/2,N*n-1)/sqrt(N*n)
  l_vint = (N*(n-1))*vint/cdfchi_x(a/2,N*(n-1))
  l_vsst = (N-1)*vsst/cdfchi_x(a/2,N-1)
  l_vtot = (N*n-1)*vtot/cdfchi_x(a/2,N*n-1)
  ;; upper limit of 95% confidence level
  u_mean = mean + avg(x*x)*cdft_t(1-a/2,N*n-1)/sqrt(N*n)
  u_vint = (N*(n-1))*vint/cdfchi_x(1-a/2,N*(n-1))
  u_vsst = (N-1)*vsst/cdfchi_x(1-a/2,N-1)
  u_vtot = (N*n-1)*vtot/cdfchi_x(1-a/2,N*n-1)
  ret = (/ (/mean,vint,vsst,vtot/), \
           (/l_mean,l_vint,l_vsst,l_vtot/), \
           (/u_mean,u_vint,u_vsst,u_vtot/), \
           (/N*n,N*(n-1),N,N*n/) /)
  return(ret)
end

begin

  ;; get d4PDF experiment information
  data_list = str_split_csv(asciiread("./data_list",14,"string")," ",0)
  expr = data_list(:,0)
  m1 = toint(data_list(:,1))
  m2 = toint(data_list(:,2))
  y1 = toint(data_list(:,3))
  y2 = toint(data_list(:,4))

  ;; mean and variance 
  var = new((/4,4,15/),float)

  ;; coordinates with "var"
  var!0 = "stats_0"
  var!1 = "stats_1"
  var!2 = "exp_type"
  var&stats_0 = ispan(0,3,1)
    ;; = (/"estimate","lower_limit","upper_limit","number_of_samples"/)
  var&stats_1 = ispan(0,3,1)
    ;; = (/"mean","vint","vsst","vtot"/)
  var&exp_type = ispan(0,14,1)
    ;; = (/"HPB","HFB_2K","HFB_4K", \
    ;;     "HFB_2K_CC","HFB_2K_GF","HFB_2K_HA","HFB_2K_MI","HFB_2K_MP","HFB_2K_MR", \
    ;;     "HFB_4K_CC","HFB_4K_GF","HFB_4K_HA","HFB_4K_MI","HFB_4K_MP","HFB_4K_MR"/)

  ;;; read AUMR and calculate mean and variance
  ;;; HIST experiment
  in = addfile("./ausm_rain_"+expr(0)+".nc","r")
  var(:,:,0) = calc_variance(in->ausm,-1)
  delete(in)
  ;;; 2K experiments as a whole
  in = addfiles("./ausm_rain_"+expr(2:7)+".nc","r")
  ListSetType(in,"cat")
  var(:,:,1) = calc_variance(in[:]->ausm,6)
  ;;; 4K experiments as a whole
  in = addfiles("./ausm_rain_"+expr(8:13)+".nc","r")
  ListSetType(in,"cat")
  var(:,:,2) = calc_variance(in[:]->ausm,6)
  delete(in)
  do i = 2, 13
    ;;; 2K (i=2-7), 4K (i=8-13) experiments
    in = addfile("./ausm_rain_"+expr(i)+".nc","r")
    var(:,:,i+1) = calc_variance(in->ausm,-1)
  end do

  ;; output
  of = "ausm_variance.nc"
  system("rm -f "+of)
  o = addfile(of,"c")
  o->var = var
end