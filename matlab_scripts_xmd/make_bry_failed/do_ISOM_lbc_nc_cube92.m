function do_ISOM_lbc_nc_cube92(interim_dir,bryname,grdname,title,obc,...
    Vtransform, Vstretching, Tcline, theta_s,theta_b,hc,N,...
    time,cycle,cmode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Built on several decades of scripting, from:
%  Pierrick Penven (Original ROMSTOOLS from ROMS_AGRIF)
%  Ben Galton-Fenzi (Major conversion for our application of ROMS)
%  Eva Cougnon (Many bug fixes)
%  David Gwyther (Many bug fixes)
%  Kaitlin Alexander (Major rewrite to use native matlab NetCDF lib)
%  and David Gwyther (Rewrite to remove caisom refs and generalise into fcn)
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...[(DayPerYear/12/2):(DayPerYear/12):(DayPerYear*(MaxYear-MinYear+1))-(DayPerYear/12/2)],[DayPerYear*(MaxYear-MinYear+1)]
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(interim_dir,'isfc_cube92.mat'))

InterpSurface = isfc;
NST = 10;
%obc=[1 1 1 1];
%time=[(DayPerYear/12/2):(DayPerYear/12):(DayPerYear*(MaxYear-MinYear+1))-(DayPerYear/12/2)];
%cycle=[DayPerYear*(MaxYear-MinYear+1)];
%cmode='clobber';
disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')

% Get vertical coords
h=ncread(grdname,'h')';
zice=ncread(grdname,'zice')';
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
lat_u=ncread(grdname,'lat_u')';
lon_u=ncread(grdname,'lon_u')';
lat_v=ncread(grdname,'lat_v')';
lon_v=ncread(grdname,'lon_v')';
mask_rho=ncread(grdname,'mask_rho')';
%mask_zice=ncread(grdname,'mask_zice')';
mask_zice=zeros(size(zice));
mask_zice(zice<0.0)=1;
angle=ncread(grdname,'angle')';

%addpath(genpath('../matlab_tools'))

h = h.*mask_rho;
zice = zice.*mask_zice;
x = lon_rho;
y = lat_rho;
kgrid = 0;
column = 0;
plt = 0;
index = 1;
[z,s_r,Cs_r]=scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, kgrid, column, index, plt);
%ang = repmat([angle(2:end,1)',angle(end,:),fliplr(angle(1:end-1,end)')],length(Cs_r),1);
%disp('STILL HAVENT ADDED ANGLE!')

% Get dims from grdfile
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;



%% Save data to more useful names

%%% rotate coordinates from roms xi,eta grid directions to geographic east,north
% (but only after averaging/interpolating u,v to the rho points) 
%uv = (u_south_latlon+1i*v_south_latlon).*exp(1i*-angle(1,:)); % 1i = sqrt(-1) in Matlab
%u_south = real(uveitheta);
%v_south = imag(uveitheta);

if obc(1)==1
sstart = 1;
send = Lp;
angle_south=repmat(angle(1,:)',[1,31,12]);
u_south_latlon =  permute(InterpSurface.u(:,:,sstart:send),[3,2,1]);
v_south_latlon =  permute(InterpSurface.v(:,:,sstart:send),[3,2,1]);
uv_south = (u_south_latlon+1i*v_south_latlon).*exp(1i*-angle_south);
u_south_rho = real(uv_south);
v_south = imag(uv_south);
u_south = 0.5.*(u_south_rho(1:end-1,:,:)+u_south_rho(2:end,:,:));

[z,s_r,Cs_r] = scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, hc, N, kgrid, column, 1, plt);
z_ext = cat(2,z,zeros(size(z(:,1))));
dr = z_ext(:,2:end)-z_ext(:,1:end-1);
dr_sum = sum(dr,2);
weights_v = dr./repmat(dr_sum,1,31);
weights_u = 0.5*(weights_v(2:end,:)+weights_v(1:end-1,:));
ubar_south = inpaint_nans(squeeze(sum(u_south.*weights_u,2)));
vbar_south = inpaint_nans(squeeze(sum(v_south.*weights_v,2)));

%ubar_south =  zeros(L,length(time));
%vbar_south =  zeros(Lp,length(time));
%zeta_south =  zeros(Lp,length(time));
temp_south =  permute(InterpSurface.tmp(:,:,sstart:send),[3,2,1]);
salt_south =  permute(InterpSurface.slt(:,:,sstart:send),[3,2,1]);
zeta_south =  permute(InterpSurface.ssh(:,sstart:send),[2,1]);
end

if obc(2)==1
estart=Lp+1;
eend=Lp+Mp;
angle_east=repmat(angle(:,end),[1,31,12]);
u_east_latlon =  permute(InterpSurface.u(:,:,estart:eend),[3,2,1]);
v_east_latlon =  permute(InterpSurface.v(:,:,estart:eend),[3,2,1]);
uv_east = (u_east_latlon+1i*v_east_latlon).*exp(1i*-angle_east);
u_east = real(uv_east);
v_east_rho = imag(uv_east);
v_east = 0.5.*(v_east_rho(1:end-1,:,:)+v_east_rho(2:end,:,:));

%ubar_east =  zeros(Mp,length(time));
%vbar_east =  zeros(M,length(time));
[z,s_r,Cs_r] = scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, hc, N, kgrid, 1, Lp, plt);
z_ext = cat(2,z,zeros(size(z(:,1))));
dr = z_ext(:,2:end)-z_ext(:,1:end-1);
dr_sum = sum(dr,2);
weights_u = dr./repmat(dr_sum,1,31);
weights_v = 0.5*(weights_u(2:end,:)+weights_u(1:end-1,:));
ubar_east = inpaint_nans(squeeze(sum(u_east.*weights_u,2)));
vbar_east = inpaint_nans(squeeze(sum(v_east.*weights_v,2)));
%zeta_east =  zeros(Mp,length(time));
temp_east =  permute(InterpSurface.tmp(:,:,estart:eend),[3,2,1]);
salt_east =  permute(InterpSurface.slt(:,:,estart:eend),[3,2,1]);
zeta_east =  permute(InterpSurface.ssh(:,estart:eend),[2,1]);
end
if obc(3)==1
nstart=Lp+Mp+1;
nend=Lp+Mp+Lp;
angle_north=repmat(angle(end,:)',[1,31,12]);
u_north_latlon =  permute(InterpSurface.u(:,:,nstart:nend),[3,2,1]);
v_north_latlon =  permute(InterpSurface.v(:,:,nstart:nend),[3,2,1]);
uv_north = (u_north_latlon+1i*v_north_latlon).*exp(1i*-angle_north);
u_north_rho = real(uv_north);
v_north = imag(uv_north);
u_north = 0.5.*(u_north_rho(1:end-1,:,:)+u_north_rho(2:end,:,:));

[z,s_r,Cs_r] = scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, hc, N, kgrid, column, Mp, plt);
z_ext = cat(2,z,zeros(size(z(:,1))));
dr = z_ext(:,2:end)-z_ext(:,1:end-1);
dr_sum = sum(dr,2);
weights_v = dr./repmat(dr_sum,1,31);
weights_u = 0.5*(weights_v(2:end,:)+weights_v(1:end-1,:));
ubar_north = inpaint_nans(squeeze(sum(u_north.*weights_u,2)));
vbar_north = inpaint_nans(squeeze(sum(v_north.*weights_v,2)));
%ubar_north =  zeros(L,length(time));
%vbar_north =  zeros(Lp,length(time));
%zeta_north =  zeros(Lp,length(time));
temp_north =  permute(InterpSurface.tmp(:,:,nstart:nend),[3,2,1]);
salt_north =  permute(InterpSurface.slt(:,:,nstart:nend),[3,2,1]);
zeta_north =  permute(InterpSurface.ssh(:,nstart:nend),[2,1]);
end
if obc(4)==1
wstart=Lp+Mp+Lp+1;
wend=Lp+Mp+Lp+Mp;
angle_west=repmat(angle(:,1),[1,31,12]);
u_west_latlon =  permute(InterpSurface.u(:,:,wstart:wend),[3,2,1]);
v_west_latlon =  permute(InterpSurface.v(:,:,wstart:wend),[3,2,1]);
uv_west = (u_west_latlon+1i*v_west_latlon).*exp(1i*-angle_west);
u_west = real(uv_west);
v_west_rho = imag(uv_west);
v_west = 0.5.*(v_west_rho(1:end-1,:,:)+v_west_rho(2:end,:,:));
%ubar_west =  zeros(Mp,length(time));
%vbar_west =  zeros(M,length(time));
[z,s_r,Cs_r] = scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, hc, N, kgrid, 1, 1, plt);
z_ext = cat(2,z,zeros(size(z(:,1))));
dr = z_ext(:,2:end)-z_ext(:,1:end-1);
dr_sum = sum(dr,2);
weights_u = dr./repmat(dr_sum,1,31);
weights_v = 0.5*(weights_u(2:end,:)+weights_u(1:end-1,:));
ubar_west = inpaint_nans(squeeze(sum(u_west.*weights_u,2)));
vbar_west = inpaint_nans(squeeze(sum(v_west.*weights_v,2)));
%zeta_west =  zeros(Mp,length(time));
temp_west =  permute(InterpSurface.tmp(:,:,wstart:wend),[3,2,1]);
salt_west =  permute(InterpSurface.slt(:,:,wstart:wend),[3,2,1]);
zeta_west =  permute(InterpSurface.ssh(:,wstart:wend),[2,1]);
end

clear InterpSurface
%%
%cmode = netcdf.getConstant('NETCDF4');
%cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
%cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));

% Save results in NetCDF file
id = netcdf.create(bryname, cmode);

% Set up dimensions
xi_u_dim = netcdf.defDim(id, 'xi_u', L);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', Lp);
eta_v_dim = netcdf.defDim(id, 'eta_v', M);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', Mp);
s_rho_dim = netcdf.defDim(id, 's_rho', N);
temp_time_dim = netcdf.defDim(id, 'temp_time', length(time));
salt_time_dim = netcdf.defDim(id, 'salt_time', length(time));
v3d_time_dim = netcdf.defDim(id, 'v3d_time', length(time));
v2d_time_dim = netcdf.defDim(id, 'v2d_time', length(time));
zeta_time_dim = netcdf.defDim(id, 'zeta_time', length(time));
one_dim = netcdf.defDim(id, 'one', 1);
num_lat_v=M;

% Set up variables
theta_s_id = netcdf.defVar(id, 'theta_s', 'double', one_dim);
netcdf.putAtt(id, theta_s_id, 'long_name', 'S-coordinate surface control parameter');
netcdf.putAtt(id, theta_s_id, 'units', 'nondimensional');
theta_b_id = netcdf.defVar(id, 'theta_b', 'double', one_dim);
netcdf.putAtt(id, theta_b_id, 'long_name', 'S-coordinate bottom control parameter');
netcdf.putAtt(id, theta_b_id, 'units', 'nondimensional');
tcline_id = netcdf.defVar(id, 'Tcline', 'double', one_dim);
netcdf.putAtt(id, tcline_id, 'long_name', 'S-coordinate surface/bottom layer width');
netcdf.putAtt(id, tcline_id, 'units', 'meter');
hc_id = netcdf.defVar(id, 'hc', 'double', one_dim);
netcdf.putAtt(id, hc_id, 'long_name', 'S-coordinate parameter, critical depth');
netcdf.putAtt(id, hc_id, 'units', 'meter');
sc_r_id = netcdf.defVar(id, 'sc_r', 'double', s_rho_dim);
netcdf.putAtt(id, sc_r_id, 'long_name', 'S-coordinate at RHO-points');
netcdf.putAtt(id, sc_r_id, 'units', 'nondimensional');
cs_r_id = netcdf.defVar(id, 'Cs_r', 'double', s_rho_dim);
netcdf.putAtt(id, cs_r_id, 'long_name', 'S-coordinate stretching curves at RHO-points');
netcdf.putAtt(id, cs_r_id, 'units', 'nondimensional');
netcdf.putAtt(id, cs_r_id, 'valid_min', -1);
netcdf.putAtt(id, cs_r_id, 'valid_max', 0);
temp_time_id = netcdf.defVar(id, 'temp_time', 'double', temp_time_dim);
netcdf.putAtt(id, temp_time_id, 'long_name', 'time for temperature climatology');
netcdf.putAtt(id, temp_time_id, 'units', 'day');
netcdf.putAtt(id, temp_time_id, 'cycle_length', cycle);
salt_time_id = netcdf.defVar(id, 'salt_time', 'double', salt_time_dim);
netcdf.putAtt(id, salt_time_id, 'long_name', 'time for temperature climatology');
netcdf.putAtt(id, salt_time_id, 'units', 'day');
netcdf.putAtt(id, salt_time_id, 'cycle_length',cycle);
v3d_time_id = netcdf.defVar(id, 'v3d_time', 'double', v3d_time_dim);
netcdf.putAtt(id, v3d_time_id, 'long_name', 'time for temperature climatology');
netcdf.putAtt(id, v3d_time_id, 'units', 'day');
netcdf.putAtt(id, v3d_time_id, 'cycle_length',cycle);
v2d_time_id = netcdf.defVar(id, 'v2d_time', 'double', v2d_time_dim);
netcdf.putAtt(id, v2d_time_id, 'long_name', 'time for temperature climatology');
netcdf.putAtt(id, v2d_time_id, 'units', 'day');
netcdf.putAtt(id, v2d_time_id, 'cycle_length', cycle);
zeta_time_id = netcdf.defVar(id, 'zeta_time', 'double', zeta_time_dim);
netcdf.putAtt(id, zeta_time_id, 'long_name', 'time for sea surface');
netcdf.putAtt(id, zeta_time_id, 'units', 'day');
netcdf.putAtt(id, zeta_time_id, 'cycle_length',cycle);
if obc(3)==1 %nth bnd
temp_north_id = netcdf.defVar(id, 'temp_north', 'double', [xi_rho_dim, s_rho_dim, temp_time_dim]);
netcdf.putAtt(id, temp_north_id, 'long_name', 'northern boundary potential temperature');
netcdf.putAtt(id, temp_north_id, 'units', 'Celsius');
netcdf.putAtt(id, temp_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, temp_north_id, '_FillValue', -1e34);
salt_north_id = netcdf.defVar(id, 'salt_north', 'double', [xi_rho_dim, s_rho_dim, salt_time_dim]);
netcdf.putAtt(id, salt_north_id, 'long_name', 'northern boundary salinity');
netcdf.putAtt(id, salt_north_id, 'units', 'PSU');
netcdf.putAtt(id, salt_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, salt_north_id, '_FillValue', -1e34);
u_north_id = netcdf.defVar(id, 'u_north', 'double', [xi_u_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, u_north_id, 'long_name', 'northern boundary u-momentum component');
netcdf.putAtt(id, u_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, u_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, u_north_id, '_FillValue', -1e34);
v_north_id = netcdf.defVar(id, 'v_north', 'double', [xi_rho_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, v_north_id, 'long_name', 'northern boundary v-momentum component');
netcdf.putAtt(id, v_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, v_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, v_north_id, '_FillValue', -1e34);
ubar_north_id = netcdf.defVar(id, 'ubar_north', 'double', [xi_u_dim, v2d_time_dim]);
netcdf.putAtt(id, ubar_north_id, 'long_name', 'northern boundary vertically integrated u-momentum component');
netcdf.putAtt(id, ubar_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, ubar_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, ubar_north_id, '_FillValue', -1e34);
vbar_north_id = netcdf.defVar(id, 'vbar_north', 'double', [xi_rho_dim, v2d_time_dim]);
netcdf.putAtt(id, vbar_north_id, 'long_name', 'northern boundary vertically integrated v-momentum component');
netcdf.putAtt(id, vbar_north_id, 'units', 'meter second-1');
netcdf.putAtt(id, vbar_north_id, 'missing_value', -1e34);
netcdf.putAtt(id, vbar_north_id, '_FillValue', -1e34);
zeta_north_id = netcdf.defVar(id, 'zeta_north', 'double', [xi_rho_dim, zeta_time_dim]);
netcdf.putAtt(id, zeta_north_id, 'long_name', 'northern boundary sea surface height');
netcdf.putAtt(id, zeta_north_id, 'units', 'meter');
end
if obc(1)==1 %sth bnd
temp_south_id = netcdf.defVar(id, 'temp_south', 'double', [xi_rho_dim, s_rho_dim, temp_time_dim]);
netcdf.putAtt(id, temp_south_id, 'long_name', 'southern boundary potential temperature');
netcdf.putAtt(id, temp_south_id, 'units', 'Celsius');
netcdf.putAtt(id, temp_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, temp_south_id, '_FillValue', -1e34);
salt_south_id = netcdf.defVar(id, 'salt_south', 'double', [xi_rho_dim, s_rho_dim, salt_time_dim]);
netcdf.putAtt(id, salt_south_id, 'long_name', 'southern boundary salinity');
netcdf.putAtt(id, salt_south_id, 'units', 'PSU');
netcdf.putAtt(id, salt_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, salt_south_id, '_FillValue', -1e34);
u_south_id = netcdf.defVar(id, 'u_south', 'double', [xi_u_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, u_south_id, 'long_name', 'southern boundary u-momentum component');
netcdf.putAtt(id, u_south_id, 'units', 'meter second-1');
netcdf.putAtt(id, u_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, u_south_id, '_FillValue', -1e34);
v_south_id = netcdf.defVar(id, 'v_south', 'double', [xi_rho_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, v_south_id, 'long_name', 'southern boundary v-momentum component');
netcdf.putAtt(id, v_south_id, 'units', 'meter second-1');
netcdf.putAtt(id, v_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, v_south_id, '_FillValue', -1e34);
ubar_south_id = netcdf.defVar(id, 'ubar_south', 'double', [xi_u_dim, v2d_time_dim]);
netcdf.putAtt(id, ubar_south_id, 'long_name', 'southern boundary vertically integrated u-momentum component');
netcdf.putAtt(id, ubar_south_id, 'units', 'meter second-1');
netcdf.putAtt(id, ubar_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, ubar_south_id, '_FillValue', -1e34);
vbar_south_id = netcdf.defVar(id, 'vbar_south', 'double', [xi_rho_dim, v2d_time_dim]);
netcdf.putAtt(id, vbar_south_id, 'long_name', 'southern boundary vertically integrated v-momentum component');
netcdf.putAtt(id, vbar_south_id, 'units', 'meter second-1');
netcdf.putAtt(id, vbar_south_id, 'missing_value', -1e34);
netcdf.putAtt(id, vbar_south_id, '_FillValue', -1e34);
zeta_south_id = netcdf.defVar(id, 'zeta_south', 'double', [xi_rho_dim, zeta_time_dim]);
netcdf.putAtt(id, zeta_south_id, 'long_name', 'southern boundary sea surface height');
netcdf.putAtt(id, zeta_south_id, 'units', 'meter');
end
if obc(2)==1 %est bnd
temp_east_id = netcdf.defVar(id, 'temp_east', 'double', [eta_rho_dim, s_rho_dim, temp_time_dim]);
netcdf.putAtt(id, temp_east_id, 'long_name', 'eastern boundary potential temperature');
netcdf.putAtt(id, temp_east_id, 'units', 'Celsius');
netcdf.putAtt(id, temp_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, temp_east_id, '_FillValue', -1e34);
salt_east_id = netcdf.defVar(id, 'salt_east', 'double', [eta_rho_dim, s_rho_dim, salt_time_dim]);
netcdf.putAtt(id, salt_east_id, 'long_name', 'eastern boundary salinity');
netcdf.putAtt(id, salt_east_id, 'units', 'PSU');
netcdf.putAtt(id, salt_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, salt_east_id, '_FillValue', -1e34);
u_east_id = netcdf.defVar(id, 'u_east', 'double', [eta_rho_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, u_east_id, 'long_name', 'eastern boundary u-momentum component');
netcdf.putAtt(id, u_east_id, 'units', 'meter second-1');
netcdf.putAtt(id, u_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, u_east_id, '_FillValue', -1e34);
v_east_id = netcdf.defVar(id, 'v_east', 'double', [eta_v_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, v_east_id, 'long_name', 'eastern boundary v-momentum component');
netcdf.putAtt(id, v_east_id, 'units', 'meter second-1');
netcdf.putAtt(id, v_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, v_east_id, '_FillValue', -1e34);
ubar_east_id = netcdf.defVar(id, 'ubar_east', 'double', [eta_rho_dim, v2d_time_dim]);
netcdf.putAtt(id, ubar_east_id, 'long_name', 'eastern boundary vertically integrated u-momentum component');
netcdf.putAtt(id, ubar_east_id, 'units', 'meter second-1');
netcdf.putAtt(id, ubar_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, ubar_east_id, '_FillValue', -1e34);
vbar_east_id = netcdf.defVar(id, 'vbar_east', 'double', [eta_v_dim, v2d_time_dim]);
netcdf.putAtt(id, vbar_east_id, 'long_name', 'eastern boundary vertically integrated v-momentum component');
netcdf.putAtt(id, vbar_east_id, 'units', 'meter second-1');
netcdf.putAtt(id, vbar_east_id, 'missing_value', -1e34);
netcdf.putAtt(id, vbar_east_id, '_FillValue', -1e34);
zeta_east_id = netcdf.defVar(id, 'zeta_east', 'double', [eta_rho_dim, zeta_time_dim]);
netcdf.putAtt(id, zeta_east_id, 'long_name', 'eastern boundary sea surface height');
netcdf.putAtt(id, zeta_east_id, 'units', 'meter');
end
if obc(4)==1 %wst bnd
temp_west_id = netcdf.defVar(id, 'temp_west', 'double', [eta_rho_dim, s_rho_dim, temp_time_dim]);
netcdf.putAtt(id, temp_west_id, 'long_name', 'western boundary potential temperature');
netcdf.putAtt(id, temp_west_id, 'units', 'Celsius');
netcdf.putAtt(id, temp_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, temp_west_id, '_FillValue', -1e34);
salt_west_id = netcdf.defVar(id, 'salt_west', 'double', [eta_rho_dim, s_rho_dim, salt_time_dim]);
netcdf.putAtt(id, salt_west_id, 'long_name', 'western boundary salinity');
netcdf.putAtt(id, salt_west_id, 'units', 'PSU');
netcdf.putAtt(id, salt_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, salt_west_id, '_FillValue', -1e34);
u_west_id = netcdf.defVar(id, 'u_west', 'double', [eta_rho_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, u_west_id, 'long_name', 'western boundary u-momentum component');
netcdf.putAtt(id, u_west_id, 'units', 'meter second-1');
netcdf.putAtt(id, u_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, u_west_id, '_FillValue', -1e34);
v_west_id = netcdf.defVar(id, 'v_west', 'double', [eta_v_dim, s_rho_dim, v3d_time_dim]);
netcdf.putAtt(id, v_west_id, 'long_name', 'western boundary v-momentum component');
netcdf.putAtt(id, v_west_id, 'units', 'meter second-1');
netcdf.putAtt(id, v_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, v_west_id, '_FillValue', -1e34);
ubar_west_id = netcdf.defVar(id, 'ubar_west', 'double', [eta_rho_dim, v2d_time_dim]);
netcdf.putAtt(id, ubar_west_id, 'long_name', 'western boundary vertically integrated u-momentum component');
netcdf.putAtt(id, ubar_west_id, 'units', 'meter second-1');
netcdf.putAtt(id, ubar_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, ubar_west_id, '_FillValue', -1e34);
vbar_west_id = netcdf.defVar(id, 'vbar_west', 'double', [eta_v_dim, v2d_time_dim]);
netcdf.putAtt(id, vbar_west_id, 'long_name', 'western boundary vertically integrated v-momentum component');
netcdf.putAtt(id, vbar_west_id, 'units', 'meter second-1');
netcdf.putAtt(id, vbar_west_id, 'missing_value', -1e34);
netcdf.putAtt(id, vbar_west_id, '_FillValue', -1e34);
zeta_west_id = netcdf.defVar(id, 'zeta_west', 'double', [eta_rho_dim, zeta_time_dim]);
netcdf.putAtt(id, zeta_west_id, 'long_name', 'western boundary sea surface height');
netcdf.putAtt(id, zeta_west_id, 'units', 'meter');

end

netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'Lateral Boundaries');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', bryname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'BOUNDARY file');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ROMS');
netcdf.endDef(id);



% Write variables
netcdf.putVar(id, theta_s_id, theta_s);
netcdf.putVar(id, theta_b_id, theta_b);
netcdf.putVar(id, tcline_id, Tcline);
netcdf.putVar(id, hc_id, hc);
netcdf.putVar(id, sc_r_id, s_r);
netcdf.putVar(id, cs_r_id, Cs_r);
netcdf.putVar(id, temp_time_id, time);
netcdf.putVar(id, salt_time_id, time);
netcdf.putVar(id, v3d_time_id, time);
netcdf.putVar(id, v2d_time_id, time);
netcdf.putVar(id, zeta_time_id, time);

if obc(1)==1
netcdf.putVar(id, temp_south_id, temp_south);
netcdf.putVar(id, salt_south_id, salt_south);
netcdf.putVar(id, u_south_id, u_south);
netcdf.putVar(id, v_south_id, v_south);
netcdf.putVar(id, ubar_south_id, ubar_south);
netcdf.putVar(id, vbar_south_id, vbar_south);
netcdf.putVar(id, zeta_south_id, zeta_south);
end

if obc(2)==1
netcdf.putVar(id, temp_east_id, temp_east);
netcdf.putVar(id, salt_east_id, salt_east);
netcdf.putVar(id, u_east_id, u_east);
netcdf.putVar(id, v_east_id, v_east);
netcdf.putVar(id, ubar_east_id, ubar_east);
netcdf.putVar(id, vbar_east_id, vbar_east);
netcdf.putVar(id, zeta_east_id, zeta_east);
end

if obc(3)==1
netcdf.putVar(id, temp_north_id, temp_north);
netcdf.putVar(id, salt_north_id, salt_north);
netcdf.putVar(id, u_north_id, u_north);
netcdf.putVar(id, v_north_id, v_north);
netcdf.putVar(id, ubar_north_id, ubar_north);
netcdf.putVar(id, vbar_north_id, vbar_north);
netcdf.putVar(id, zeta_north_id, zeta_north);
end
if obc(4)==1
netcdf.putVar(id, temp_west_id, temp_west);
netcdf.putVar(id, salt_west_id, salt_west);
netcdf.putVar(id, u_west_id, u_west);
netcdf.putVar(id, v_west_id, v_west);
netcdf.putVar(id, ubar_west_id, ubar_west);
netcdf.putVar(id, vbar_west_id, vbar_west);
netcdf.putVar(id, zeta_west_id, zeta_west);
end

netcdf.close(id);

disp('done!')
