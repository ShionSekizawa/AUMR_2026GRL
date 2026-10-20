;;; File: cov.ncl
;;; Description: Calculate statistics of AUMR and U10 in CMIP6 historical simulations
;;;
;;; Input: sst_pc_historical.nc, ausm_historical.nc, u10_historical.nc
;;; Output: u10_cov.txt
;;;
;;; Running:
;;; ncl cov.ncl
;;;

begin

  m_tmp  = asciiread("ModelMemList.txt",26,"string")
  m_name = new(26,string)
  do i = 0, 26-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do

  in = addfile("sst_pc_historical.nc","r")
  pc1 = in->pc(:,0,:)
  in = addfile("ausm_historical.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  ausm_mean = dim_avg(in->ausm(:,id))
  dum = dim_rmvmean(in->ausm(:,id))
  dum = dtrend(dum,False)
  do i = 0, 26-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  ausm_std = dim_stddev(dum)
  in = addfile("u10_historical.nc","r")
  u10_mean = dim_avg(in->u10(:,id))
  dum2 = dim_rmvmean(in->u10(:,id))
  dum2 = dtrend(dum2,False)
  do i = 0, 26-1
    dum2(i,:) = dtrend_msg(pc1(i,:),dum2(i,:),False,False)
    dum(i,:) = dum(i,:)/stddev(dum(i,:))
  end do
  u10_reg = regCoef_n(dum,dum2,1,1)
  print(m_name+": aumr_mean="+sprintf("%3.2f",ausm_mean)+" aumr_stddev="+sprintf("%3.2f",ausm_std)+" u10_mean="+sprintf("%3.2f",u10_mean)+" u10_reg="+sprintf("%3.2f",u10_reg)+" ratio:"+sprintf("%4.3f",u10_reg/u10_mean))

  asciiwrite("u10_cov.txt",u10_reg/u10_mean)
end