% make_lbc.m
% David Gwyther | 2014-5
% changelog:
% 2014:		 written
% 2015-Sep-08: 	 updated to remove 'interactive' mode. I don't EVER use it.
% 2015-Dec-01:   Updated to remove a bunch of redundant features
% 2017-Jun-21:   adapted to work on 360° whole antarctic domain

run = 'waom10';
%addpath(genpath('/ds/projects/iomp/matlab_scripts'))
proj_dir = fullfile('..','..');
data_dir = fullfile(proj_dir,'data_xmd','preprocessing');
interim_dir = fullfile(data_dir,'interim_xmd');
processed_dir = fullfile(data_dir,'processed_xmd');
external_dir = fullfile(data_dir,'external_xmd');

grdname = fullfile(processed_dir,[run,'_grd.nc']);
bryname = fullfile(processed_dir,[run,'_bry_1996_2016.nc']);
MinYear = 1996;
MaxYear = 2016;
ECCObounds = [1 1438 28 200];%as [xmin xmax ymin ymax]; 
RunName = run;

Vtransform = 2;
Vstretching = 4;
theta_s = 7;
theta_b = 8;
Tcline = 20;
N = 31;

%%%%%%%%%%%%%%%%%%%
% Totten is [410 525 80 125]; %as [xmin xmax ymin ymax];
% Amery (new) is [210 400 50 125];

%DataProduct = 4; 
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube84 
% (4) cube92 monthly
% (5) cube92 3-daily
% (6) Dinniman ACCIMA 5-km
% (7) O'Kane monthly ocean model (100 normal years) - COREv1 forced
% (8) O'Kane monthly ocean model (1948-2006) - COREv2 forced


%ForcingType = 3;
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube92 interannual w/ 365 day year
% (4) cube92 climatology
% (5) cube92 constant forcing
% (6) cube92 climatology w/ 364 day year
% (7) cube92 constant summer forcing
% (8) cube84 climatology
% (9) cube92-3day climatology w/ 364 day year
%(10) cube92-3day climatology
%(11) ACCIMA 5-km climatology
%(12) ACCIMA 5-km IAF
%(13) O'Kane-COREv1 monthly interannual (model yrs 1900-2000)
%(14) O'Kane-COREv1 monthly interannual climatology (model yrs 1900-2000)
%(15) O'Kane-COREv2 monthly interannual (1948-2006) - COREv2 forced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

force_id = {'ssh','salt','theta','uvel','vvel'}; %load and interp ecco2 data
xmin = ECCObounds(1); 
xmax = ECCObounds(2);
ymin = ECCObounds(3); 
ymax = ECCObounds(4);

%for ii = 5
%    disp(['loading ' force_id{ii} ' data'])
%    eval(['do_load_ecco2_',force_id{ii},'_cube92'])
%end


%% REGRID ECCO2 TO MODEL & run do_*isom_lbc.m



do_ISOM_lbc_cube92

%do_ISOM_lbc_nc_cube92(interim_dir,bryname,grdname,'Lateral Boundaries Salt and Temp',[1 1 1 1], Vtransform, Vstretching, Tcline, theta_s,theta_b,hc,N,[(DayPerYear/12/2):(DayPerYear/12):(DayPerYear*(MaxYear-MinYear+1))-(DayPerYear/12/2)],[DayPerYear*(MaxYear-MinYear+1)],'clobber')




disp('done!')
    
