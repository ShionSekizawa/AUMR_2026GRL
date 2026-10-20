;;; File: cov_mulens.ncl
;;; Description: Calculate statistics of AUMR and U10 in CMIP6 historical simulations
;;;
;;; Input: sst_pc_historical_mulens.nc, ausm_historical_mulens.nc, u10_historical_mulens.nc
;;; Output: u10_cov_mulens.txt
;;;
;;; Running:
;;; ncl cov_mulens.ncl
;;;

begin

  m_tmp  = asciiread("../ModelMemList.txt",26,"string")
  m_name = new(26,string)
  do i = 0, 26-1
    tmp = str_split(m_tmp(i)," ")
    m_name(i) = tmp(1)
    delete(tmp)
  end do
  mm_tmp = asciiread("ModelMemList_mulens.txt",13*5,"string")
  m_name_index = new(13,integer)
  do i = 0, 13-1
    tmp = str_split(mm_tmp(i*5)," ")
    m_name_index(i) = ind(m_name.eq.tmp(1))
    delete(tmp)
  end do
  m_name := m_name(m_name_index)

  in = addfile("sst_pc_historical_mulens.nc","r")
  pc1 = reshape(in->pc(:,:,0,:),(/13,5*60/))
  in = addfile("ausm_historical_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  ausm_mean = dim_avg_n(in->ausm(:,:,id),(/1,2/))
  dum = reshape(dim_rmvmean_n(in->ausm(:,:,id),(/1,2/)),(/13,5*60/))
  dum = dtrend_msg(ndtooned(conform_dims((/5,60/),yyyy(id),1)),dum,False,False)
  do i = 0, 13-1
    dum(i,:) = dtrend_msg(pc1(i,:),dum(i,:),False,False)
  end do
  ausm_std = dim_stddev(dum)
  in = addfile("u10_historical_mulens.nc","r")
  yyyy = cd_calendar(in->time,-1)/100
  id = ind(yyyy.ge.1951.and.yyyy.le.2010)
  u10_mean = dim_avg_n(in->u10(:,:,id),(/1,2/))
  dum2 = reshape(dim_rmvmean_n(in->u10(:,:,id),(/1,2/)),(/13,5*60/))
  dum2 = dtrend_msg(ndtooned(conform_dims((/5,60/),yyyy(id),1)),dum2,False,False)
  do i = 0, 13-1
    dum2(i,:) = dtrend_msg(pc1(i,:),dum2(i,:),False,False)
    dum(i,:) = dum(i,:)/stddev(dum(i,:))
  end do
  u10_reg = regCoef_n(dum,dum2,1,1)
  print(m_name+": aumr_mean="+sprintf("%3.2f",ausm_mean)+" aumr_stddev="+sprintf("%3.2f",ausm_std)+" u10_mean="+sprintf("%3.2f",u10_mean)+" u10_reg="+sprintf("%3.2f",u10_reg)+" ratio:"+sprintf("%4.3f",u10_reg/u10_mean))

  asciiwrite("u10_cov_mulens.txt",u10_reg/u10_mean)
end