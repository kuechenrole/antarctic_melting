run = 'waom4'
addpath ./TMD
addpath ./TMD/FUNCTIONS
addpath ./t_tide_v1.3beta

proj_dir = fullfile('..','..');

data_dir = fullfile(proj_dir,'data','preprocessing');
gfile=fullfile(data_dir,'processed',[run,'_grd.nc']);
base_date=datenum(2007,1,1);
pred_date=datenum(2015,10,2);
ofile=fullfile(data_dir,'processed',[run,'_tds_nodal_2015.nc']);
model_file=fullfile(data_dir,'external','tpxo','Model_tpxo7.2');
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,run)
