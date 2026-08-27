This repository contains NCAR Command Language (NCL) scripts to reproduce the figures in the paper:
SST Warming Patterns Modulate the Amplitude of Interannual Variability of Australian Monsoon Rainfall under Global Warming (Sekizawa 2026, GRL).

## Data
The analyses uses:
  - d4PDF and JRA-3Q from the Data Integration and Analysis System (DIAS; https://search.diasjp.net/ja/dataset/d4PDF_GCM, https://search.diasjp.net/ja/dataset/JRA3Q)
  - CMIP6 from the Earth System Grid Federation (ESGF; https://metagrid.esgf-west.org/search/cmip6/)
  - AGCD from the Australian Bureau of Meteorology (https://www.bom.gov.au/climate/austmaps/about-agcd-maps.shtml)
  - COBE-SST2 from the Japan Meteorological Agency (https://ds.data.jma.go.jp/goos/data/rrtdb/jma-pro.html)

Due to data size, these data are not included here.


## Tools
NCAR Command Language (NCL; http://doi.org/10.5065/D6WD3XH5) version 6.5.0 has been used for the analyses and visualizations. 


##  File list
./d4PDF (includes files for making processed d4PDF data)
  (Text file including information on the d4PDF experiments)
  - data_list
  (NCL scripts)
  - clim_seasonality.ncl
  - derive_MSEbudget_advect.ncl
  - derive_MSEbudget_forcing.ncl
  - derive_ausm_variance.ncl
  - make_MSEbudget_area.ncl
  - make_ausm_rain.ncl
  - make_lhf_area.ncl
  - make_sst_clim.ncl
  - make_sst_hist.ncl
  - make_vadv_area.ncl
  - make_figS2_tmpfile_d4pdf.ncl

./CMIP6 (includes files for making processed CMIP6 data)
  (Text files of CMIP6 data file list)
  - FileList_hfls_historical.txt
  - FileList_hfls_ssp245.txt
  - FileList_hfls_ssp585.txt
  - FileList_pr_historical.txt
  - FileList_pr_ssp245.txt
  - FileList_pr_ssp585.txt
  - FileList_sftlf.txt
  - FileList_tas_historical.txt
  - FileList_tas_ssp245.txt
  - FileList_tas_ssp585.txt
  - FileList_ts_historical.txt
  - FileList_ts_ssp245.txt
  - FileList_ts_ssp585.txt
  - FileList_uas_historical.txt
  - FileList_uas_ssp245.txt
  - FileList_uas_ssp585.txt
  - FileList_vas_historical.txt
  - FileList_vas_ssp245.txt
  - FileList_vas_ssp585.txt
  - ModelMemList.txt
  (NCL scripts)
  - ausm.ncl
  - gmst.ncl
  - gmst_pi.ncl
  - sst.ncl
  - u10.ncl
  - gwlev.ncl
  - sst_pc.ncl
  - uv10.ncl
  - hfls.ncl

./obs (includes files for making processed observation data)
  (NCL scripts)
  - make_ausm_rain.obs.ncl
  - sst_pc.obs.ncl
  - make_figS2_tmpfile_obs.ncl

./figure: (includes files for creating figures)
  (NCL scripts)
  - fig1.ncl
  - fig2.ncl
  - fig3.ncl
  - fig4.ncl
  - figS1.ncl
  - figS2.ncl
  - figS3.ncl
  - figS4.ncl


## by Dr. Shion Sekizawa, Meteorological Research Institute, Japan Meteorological Agency (shion{at}mri-jma.go.jp)
## Aug 27, 2026