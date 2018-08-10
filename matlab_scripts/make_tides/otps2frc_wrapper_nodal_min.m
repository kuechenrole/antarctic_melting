run = 'waom4'
addpath ./TMD
addpath ./TMD/FUNCTIONS
addpath ./t_tide_v1.3beta
addpath [getenv('extdir'),'/tpxo']

gfile=[getenv('prodir'),'/',run,'_grd.nc']
base_date=datenum(2007,1,1);
pred_date=datenum(2015,10,2);
ofile=[getenv('prodir'),'/',run,'nodal_min_tds.nc'];
model_file=[getenv('extdir'),'/tpxo/Model_tpxo7.2'];
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,run)
