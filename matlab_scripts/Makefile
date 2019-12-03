######## preprocessing steps to generate grids, boundary conditions and forcing files for ROMS whole antarctic ocean applications #############

####### summaries ti run e.g. make 10km and get all the necescary files #############
10km: $(intdir)/waom10_grd.nc $(prodir)/waom10_ini.nc $(prodir)/waom10_bry.nc $(prodir)/waom10_sbc.nc $(prodir)/waom10_nudge.nc $(prodir)/waom10_tds.nc

5km: $(intdir)/waom5_grd.nc $(prodir)/waom5_ini.nc $(prodir)/waom5_bry.nc $(prodir)/waom5_sbc.nc $(prodir)/waom5_nudge.nc $(prodir)/waom5_tds.nc

grids: $(intdir)/waom10_grd.nc $(intdir)/waom5_grd.nc
ini: $(prodir)/waom10_ini.nc $(prodir)/waom5_ini.nc
bry: $(prodir)/waom10_bry.nc $(prodir)/waom5_bry.nc 
sbc: $(prodir)/waom10_sbc.nc  $(prodir)/waom10_nudge.nc $(prodir)/waom5_sbc.nc $(prodir)/waom5_nudge.nc
tds: $(prodir)/waom10_tds.nc $(prodir)/waom5_tds.nc

####### individual files ########

#### GRIDS ######
# see notebook reports/2.5-kuechenrole-grid_tuning.ipynb for additional tuning steps to increase time step


$(prodir)/waom1_grd.nc: $(extdir)/bedmap/bedmap2.mat $(extdir)/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc grid/make_grid.m
	sed -i 's/mr = ../mr =  1/g' grid/make_grid_wrapper.m
	matlab -nodesktop -nosplash -r 'run grid/make_grid_wrapper; quit()'
	mv $(intdir)/waom1_MinDepth20m_rx10.3_grd.nc $@

$(prodir)/waom5_grd.nc: $(intdir)/waom5_grd.nc grid/tune_waom5_grd.py
	python grid/tune_waom5_grd.py

$(intdir)/waom5_grd.nc: $(extdir)/bedmap/bedmap2.mat $(extdir)/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc grid/make_grid.m
	sed -i 's/mr = ../mr =  5/g' grid/make_grid_wrapper.m
	matlab -nodesktop -nosplash -r 'run grid/make_grid_wrapper; quit()'

$(prodir)/waom10_grd.nc: $(intdir)/waom10_grd.nc grid/tune_waom10_grd.py
	python grid/tune_waom10_grd.py

$(intdir)/waom10_grd.nc: $(extdir)/bedmap/bedmap2.mat $(extdir)/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc grid/make_grid.m
	sed -i 's/mr = ../mr = 10/g' grid/make_grid_wrapper.m
	matlab -nodesktop -nosplash -r 'run grid/make_grid_wrapper; quit()'


$(extdir)/bedmap/bedmap2.mat:
$(extdir)/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc:



##### Initial conditions ########

$(prodir)/waom10_ini.nc: $(extdir)/SALT.nc/SALT.1440x720x50.20070101.nc $(extdir)/THETA.nc/THETA.1440x720x50.20070101.nc ini/make_ini.py
	python ini/make_ini.py

$(extdir)/SALT.nc/SALT.1440x720x50.20070101.nc:
$(extdir)/THETA.nc/THETA.1440x720x50.20070101.nc:

$(prodir)/waom5_ini.nc: $(prodir)/waom5_grd.nc $(prodir)/waom10_grd.nc $(prodir)/waom10_his_0006.nc $(prodir)/waom10_ini.nc ini/make_hr_ini.py
	python ini/make_hr_ini.py

$(prodir)/waom10_his_0006.nc:



####### open boundary conditions #########

$(prodir)/waom10_bry.nc: $(prodir)/waom10_grd.nc $(extdir)/ecco2/SALT.nc/SALT.1440x720x50.200701.nc $(extdir)/ecco2/THETA.nc/THETA.1440x720x50.200701.nc $(extdir)/ecco2/SSH.nc/SSH.1440x720.2007_monmean.nc lbc/make_lbc.m
	sed -i 's/waom.*_grd.nc/waom10_grd.nc/g' lbc/make_lbc.m
	sed -i 's/waom.*_bry.nc/waom10_bry.nc/g' lbc/make_lbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom10'/g" lbc/make_lbc.m
	matlab -nodesktop -nosplash -r 'run lbc/make_lbc.m; quit()'

$(prodir)/waom5_bry.nc: $(prodir)/waom5_grd.nc $(extdir)/ecco2/SALT.nc/SALT.1440x720x50.200701.nc $(extdir)/ecco2/THETA.nc/THETA.1440x720x50.200701.nc $(extdir)/ecco2/SSH.nc/SSH.1440x720.2007_monmean.nc lbc/make_lbc.m
	sed -i 's/waom.*_grd.nc/waom5_grd.nc/g' lbc/make_lbc.m
	sed -i 's/waom.*_bry.nc/waom5_bry.nc/g' lbc/make_lbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom5'/g" lbc/make_lbc.m
	matlab -nodesktop -nosplash -r 'run lbc/make_lbc.m; quit()'

$(prodir)/waom2_bry.nc: $(prodir)/waom2_grd.nc $(extdir)/ecco2/SALT.nc/SALT.1440x720x50.200701.nc $(extdir)/ecco2/THETA.nc/THETA.1440x720x50.200701.nc $(extdir)/ecco2/SSH.nc/SSH.1440x720.2007_monmean.nc lbc/make_lbc.m
	sed -i 's/waom.*_grd.nc/waom2_grd.nc/g' lbc/make_lbc.m
	sed -i 's/waom.*_bry.nc/waom2_bry.nc/g' lbc/make_lbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom2'/g" lbc/make_lbc.m
	matlab -nodesktop -nosplash -r 'run lbc/make_lbc.m; quit()'

$(extdir)/ecco2/SALT.nc/SALT.1440x720x50.200701.nc:
$(extdir)/ecco2/THETA.nc/THETA.1440x720x50.200701.nc:
$(extdir)/ecco2/SSH.nc/SSH.1440x720.2007_monmean.nc:



######## surface boundary conditions ###########

### sbc files ###

$(prodir)/waom10_sbc.nc: $(prodir)/waom10_grd.nc $(extdir)/tamura/TSDM2hb_2007_apr.data $(extdir)/era_interim/ERA_Interim_1992_2011.2daily.uwinds.nc sbc/make_sbc.m sbc/do_ISOM_sbc_nc_daily_IAF.m sbc/ERA_interim_misom_grid_stress_annual.m
	sed -i 's/waom.*_grd.nc/waom10_grd.nc/g' sbc/make_sbc.m
	sed -i 's/waom.*_sbc.nc/waom10_sbc.nc/g' sbc/make_sbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom10'/g" sbc/make_sbc.m
	matlab -nodesktop -nosplash -r 'run sbc/make_sbc.m; quit()'

$(prodir)/waom5_sbc.nc: $(prodir)/waom5_grd.nc $(extdir)/tamura/TSDM2hb_2007_apr.data $(extdir)/era_interim/ERA_Interim_1992_2011.2daily.uwinds.nc sbc/make_sbc.m sbc/do_ISOM_sbc_nc_daily_IAF.m sbc/ERA_interim_misom_grid_stress_annual.m
	sed -i 's/waom.*_grd.nc/waom5_grd.nc/g' sbc/make_sbc.m
	sed -i 's/waom.*_sbc.nc/waom5_sbc.nc/g' sbc/make_sbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom5'/g" sbc/make_sbc.m
	matlab -nodesktop -nosplash -r 'run sbc/make_sbc.m; quit()'

$(prodir)/waom2_sbc.nc: $(prodir)/waom2_grd.nc $(extdir)/tamura/TSDM2hb_2007_apr.data $(extdir)/era_interim/ERA_Interim_1992_2011.2daily.uwinds.nc sbc/make_sbc.m sbc/do_ISOM_sbc_nc_daily_IAF.m sbc/ERA_interim_misom_grid_stress_annual.m
	sed -i 's/waom.*_grd.nc/waom2_grd.nc/g' sbc/make_sbc.m
	sed -i 's/waom.*_sbc.nc/waom2_sbc.nc/g' sbc/make_sbc.m
	sed -i "s/RunName = 'waom.*'/RunName = 'waom2'/g" sbc/make_sbc.m
	matlab -nodesktop -nosplash -r 'run sbc/make_sbc.m; quit()'

$(extdir)/tamura/TSDM2hb_2007_apr.data:
$(extdir)/era_interim/ERA_Interim_1992_2011.2daily.uwinds.nc:


#### surface nudging files ###

$(prodir)/waom10_nudge.nc: $(prodir)/waom10_grd.nc $(extdir)/sose/SALT_mnthlyBar.0000000100.data sbc/make_nudge.py
	sed -i 's/waom.*_grd.nc/waom10_grd.nc/g' sbc/make_nudge.py
	sed -i 's/waom.*_nudge.nc/waom10_nudge.nc/g' sbc/make_nudge.py
	python sbc/make_nudge.py

$(prodir)/waom5_nudge.nc: $(prodir)/waom5_grd.nc $(extdir)/sose/SALT_mnthlyBar.0000000100.data sbc/make_nudge.py
	sed -i 's/waom.*_grd.nc/waom5_grd.nc/g' sbc/make_nudge.py
	sed -i 's/waom.*_nudge.nc/waom5_nudge.nc/g' sbc/make_nudge.py
	python sbc/make_nudge.py

$(prodir)/waom2_nudge.nc: $(prodir)/waom2_grd.nc $(extdir)/sose/SALT_mnthlyBar.0000000100.data sbc/make_nudge.py
	sed -i 's/waom.*_grd.nc/waom2_grd.nc/g' sbc/make_nudge.py
	sed -i 's/waom.*_nudge.nc/waom2_nudge.nc/g' sbc/make_nudge.py
	python sbc/make_nudge.py

$(extdir)/sose/SALT_mnthlyBar.0000000100.data:




######## tidal forcing ###########

$(prodir)/waom10_tds.nc: $(intdir)/waom10_tds.nc tides/mask_onshelf.py
	sed -i 's/waom.*_grd.nc/waom10_grd.nc/g' tides/mask_onshelf.py
	sed -i 's/waom.*_tds.nc/waom10_tds.nc/g' tides/mask_onshelf.py
	python tides/mask_onshelf.py

$(intdir)/waom10_tds.nc: $(prodir)/waom10_grd.nc tides/DATA/Model_tpxo7.2 
	sed -i 's/waom.*_grd.nc/waom10_grd.nc/g' tides/otps2frc_wrapper.m
	sed -i 's/waom.*_tds.nc/waom10_tds.nc/g' tides/otps2frc_wrapper.m
	sed -i "s/'WAOM.*'/'WAOM10'/g" tides/otps2frc_wrapper.m
	matlab -nodesktop -nosplash -r 'run tides/otps2frc_wrapper;quit()'
	ncrename -O -v tide_Ephase,tide_Pphase $@ 
	ncrename -O -v tide_Eamp,tide_Pamp $@

$(prodir)/waom5_tds.nc: $(intdir)/waom5_tds.nc tides/mask_onshelf.py
	sed -i 's/waom.*_grd.nc/waom5_grd.nc/g' tides/mask_onshelf.py
	sed -i 's/waom.*_tds.nc/waom5_tds.nc/g' tides/mask_onshelf.py
	python tides/mask_onshelf.py

$(intdir)/waom5_tds.nc: $(prodir)/waom5_grd.nc tides/DATA/Model_tpxo7.2 
	sed -i 's/waom.*_grd.nc/waom5_grd.nc/g' tides/otps2frc_wrapper.m
	sed -i 's/waom.*_tds.nc/waom5_tds.nc/g' tides/otps2frc_wrapper.m
	sed -i "s/'WAOM.*'/'WAOM5'/g" tides/otps2frc_wrapper.m
	matlab -nodesktop -nosplash -r 'run tides/otps2frc_wrapper;quit()'
	ncrename -O -v tide_Ephase,tide_Pphase $@ 
	ncrename -O -v tide_Eamp,tide_Pamp $@

$(prodir)/waom2_tds.nc: $(intdir)/waom2_tds.nc tides/mask_onshelf.py
	sed -i 's/waom.*_grd.nc/waom2_grd.nc/g' tides/mask_onshelf.py
	sed -i 's/waom.*_tds.nc/waom2_tds.nc/g' tides/mask_onshelf.py
	python tides/mask_onshelf.py

$(intdir)/waom2_tds.nc: $(prodir)/waom2_grd.nc tides/DATA/Model_tpxo7.2 
	sed -i 's/waom.*_grd.nc/waom2_grd.nc/g' tides/otps2frc_wrapper.m
	sed -i 's/waom.*_tds.nc/waom2_tds.nc/g' tides/otps2frc_wrapper.m
	sed -i "s/'WAOM.*'/'WAOM2'/g" tides/otps2frc_wrapper.m
	matlab -nodesktop -nosplash -r 'run tides/otps2frc_wrapper;quit()'
	ncrename -O -v tide_Ephase,tide_Pphase $@ 
	ncrename -O -v tide_Eamp,tide_Pamp $@

tides/DATA/Model_tpxo7.2:

