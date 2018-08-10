%function otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,dname)
%otps2frc_v5 Generete OTPS tidal forcing file for ROMS
%
%otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,dname) generates a ROMS tidal
%forcing  ofile using the ROMS gridfile gfile and the tidal
%reference time base_date in matlab time. pred_date is the
%prediction timde for nodal corrections, typically the center time
%of a two year prediction period, also in MATLAB time.model_file is
%the path to the directory contained to appropriate OTPS model output
%dname is a domain name. 
%
% %Example:
% gfile='espresso_grid_c05.nc'
% base_date=datenum(2006,1,1);
% pred_date=datenum(2006,1,1);
% ofile='tidetest.nc';
% model_file='DATA/Model_EC';
% otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'ESPRESSO')
%
%Requirements:
%T_TIDE tidal analysis package
%"tidal_ellipse" (Zhigang Xu) package, ap2ep.m
%Requires Oregon TMD_toolbox
%
%Originally written by John Evans
%Revised by Eli Hunter 3/7/07
%Revised by Eli Hunter 5/25/07 (corrected phase lag error)
% REVISED by Eli Hunter 8/2/2011 Added 'nodal', to t_vuf and changed Tide.period
% REVISED by Eli Hunter 6/10/2016 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collect the information necessary to run extract_HC and
%create setup files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars={'z','u','v'};


if ~(exist(model_file,'file'))
    error(['No Such model file: ' model_file])
end

lon=ncread(gfile,'lon_rho');
lat=ncread(gfile,'lat_rho');


disp(['Mode parameter file used: ' model_file])

disp(['Interpolating ' vars{1}])
[z_amp,z_phase,Depth,conList]=tmd_extract_HC(model_file,lat,lon,vars{1},[1 2]);
disp(['Interpolating ' vars{2}])
[u_amp,u_phase,Depth,conList]=tmd_extract_HC(model_file,lat,lon,vars{2},[1 2]);
disp(['Interpolating ' vars{3}])
[v_amp,v_phase,Depth,conList]=tmd_extract_HC(model_file,lat,lon,vars{3},[1 2]);

conList
mask_rho = ncread( gfile, 'mask_rho' );
land = find(mask_rho==0);
water = find(mask_rho==1);
%
%
%
 cnames=upper(char(conList));


% % Make sure that the OTPS mask agrees with the ROMS mask.
% % Fill in any points that ROMS thinks is water but OTPS thinks is land.
num_constituents = 2;
a=t_getconsts;
%
for j = 1:num_constituents
    
        iconst(j)=strmatch(cnames(j,:), a.name);
    	Tide.period(j)=1/a.freq(iconst(j));
    
end
%
%
 Tide.names = cnames;
 Ntide = length(Tide.period);
 [Lp,Mp] = size(lon);
%
%
%
% %***********************************************************************
% % This is the call to t_vuf that
% % will correct the phase to be at the user specified time.  Also, the amplitude
% % is corrected for nodal adjustment.
%
% % Reference latitude for 3rd order satellites (degrees) is
% % set to 55.  You don't need to adjust this to your local latitude
% % It could also be set to NaN as in Xtide, with very little effect.
% % See T_VUF for more info.
%
reflat=55;
datestr(base_date)

[V,U,F]=t_vuf('nodal',base_date,iconst,reflat);
[Vp,Up,Fp]=t_vuf('nodal',pred_date,iconst,reflat);%Only used for nodal correction.

%vv and uu are returned in cycles, so * by 360 to get degrees or * by 2 pi to get radians


V=V*360;  % convert vv to phase in degrees
U=U*360;  % convert uu to phase in degrees
Vp=Vp*360;  % convert vv to phase in degrees
Up=Up*360;  % convert uu to phase in degrees


for k=1:2;
    z_phase(k,:,:) = z_phase(k,:,:) - Up(k)  - V(k);   % degrees
    z_amp(k,:,:) =z_amp(k,:,:) .* Fp(k);

    u_phase(k,:,:) =u_phase(k,:,:) - Up(k) - V(k);   % degrees
    u_amp(k,:,:) = u_amp(k,:,:) .* Fp(k);

    v_phase(k,:,:) = v_phase(k,:,:) - Up(k)  - V(k);   % degrees
    v_amp(k,:,:) =v_amp(k,:,:) .* Fp(k);

end
%
%
 
 z_phase=mod(z_phase,360);
 u_phase=mod(u_phase,360);
 v_phase=mod(v_phase,360);
 
 z_amp = zero_out_land ( z_amp, land );
 z_phase = zero_out_land ( z_phase, land );
%
 Tide.Ephase    = z_phase(:,:,:);
 Tide.Eamp      = z_amp(:,:,:);

%
%
%
% %---------------------------------------------------------------------
% %  Convert tidal current amplitude and phase lag parameters to tidal
% %  current ellipse parameters: Major axis, ellipticity, inclination,
% %  and phase.  Use "tidal_ellipse" (Zhigang Xu) package.
% %---------------------------------------------------------------------
%%
major = zeros(size(u_amp));
eccentricity = zeros(size(u_amp));
inclination = zeros(size(u_amp));
phase = zeros(size(u_amp));

for j = 1:num_constituents;
    [major(j,:,:),eccentricity(j,:,:),inclination(j,:,:),phase(j,:,:)]=ap2ep(u_amp(j,:,:),u_phase(j,:,:),v_amp(j,:,:),v_phase(j,:,:));
end
%[major,eccentricity,inclination,phase]=ap2ep(u_amp,u_phase,v_amp,v_phase);
major = zero_out_land ( major, land );
eccentricity = zero_out_land ( eccentricity, land );
major = major/100;
Tide.Cmax=major;
Tide.Cmin=major.*eccentricity;
Tide.Cangle= zero_out_land ( inclination, land );
Tide.Cphase = zero_out_land ( phase, land );
%
%
 write_roms_otps_ncfile_v2( Tide, gfile, ofile,base_date,model_file,dname);
