%% Load ECCO data
disp('Load ECCO data')



m1 = matfile(fullfile(interim_dir,['cube92_iaf_salt_',RunName,'.mat']),'Writable',true)
m2 = matfile(fullfile(interim_dir,['cube92_iaf_theta_',RunName,'.mat']),'Writable',true)
m3 = matfile(fullfile(interim_dir,['cube92_iaf_uvel_',RunName,'.mat']),'Writable',true)
m4 = matfile(fullfile(interim_dir,['cube92_iaf_vvel_',RunName,'.mat']),'Writable',true)
m5 = matfile(fullfile(interim_dir,['cube92_iaf_ssh_',RunName,'.mat']),'Writable',true)

%salt = m1.salt(:,:,:,:);
%theta = m2.theta(:,:,:,:);
%ssh = m5.ssh(:,:,:);
%uvel = m3.uvel(:,:,:,:);
%vvel = m4.vvel(:,:,:,:);

disp('done')

%theta(theta < -10) = NaN;
%salt(salt < 0) = NaN;

%%
% Load model grid:
%ncload(grdname,'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v')

% Find nearest indexes for locations:
%disp('Find nearest indexes for locations')
%ncload /u/crcdata/ECCO2/cube84/THETA/THETA.1440x720x50.001.nc LATITUDE_T LONGITUDE_T DEPTH_T
ecco_grd =  [external_dir,'/ecco2/THETA.nc/THETA.1440x720x50.19960102.nc'];
LATITUDE_T = ncread(ecco_grd,'LATITUDE_T')';
LONGITUDE_T = ncread(ecco_grd,'LONGITUDE_T')';
DEPTH_T = ncread(ecco_grd,'DEPTH_T')';

disp('done')

%[lons lats] = meshgrid(LONGITUDE_T,LATITUDE_T);

%depth = DEPTH_T;
% % For TISOM:
% xmax = 525;
% xmin = 410; % Indices for this are made in make_lbc.m
% ymax = 125;
% ymin = 80;

Xloc = [xmin:xmax];
Yloc = [ymin:ymax];
%{
Var1 = theta;
Var2 = salt;
%}
%X = lons(Yloc,Xloc);
%Y = lats(Yloc,Xloc);

%% Calculate vertical levels:
%ncload(grdname,'h','zice','lat_rho','lon_rho','mask_rho','mask_zice')
h=ncread(grdname,'h')';
zice=ncread(grdname,'zice')';
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
%lat_u=ncread(grdname,'lat_u')';
%lon_u=ncread(grdname,'lon_u')';
%lat_v=ncread(grdname,'lat_v')';
%lon_v=ncread(grdname,'lon_v')';
mask_rho=ncread(grdname,'mask_rho')';
%mask_zice=ncread(grdname,'mask_zice')';
mask_zice=zeros(size(zice));
mask_zice(zice<0.0)=1;

addpath(genpath('../matlab_tools/'));

h = h.*mask_rho;
zice = zice.*mask_zice;
hc = Tcline;
x = lon_rho;
y = lat_rho;
% Vtransform = 1;
% Vstretching = 2;
% theta_s = 0.9;
% theta_b = 4;
%Vtransform = 2;
%Vstretching = 4;
%theta_s = 0.9;
%theta_b = 4;
%N = 31;
kgrid = 0;
column = 0;
plt = 0;
index = 1;
[z,s,Cs_r]=scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, kgrid, column, index, plt);
kgrid = 1;
[z,s,Cs_w]=scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, 1, column, index, plt);

%% For rho points:

InterpSurfaceIni    = nan(length(Cs_r),size(h,1)*2+size(h,2)-2);
TotalWCT            = [h(1,1:end),h(1:end,end)',h(end,1:end),h(1:end,1)'];%[h(2:end,1)',h(end,:),fliplr(h(1:end-1,end)')];
MaskVerSec          = repmat([mask_rho(2:end,1)',mask_rho(end,:),fliplr(mask_rho(1:end-1,end)')],length(Cs_r),1);
ind1 = find(MaskVerSec == 0); MaskVerSecNaN = MaskVerSec;
MaskVerSecNaN(ind1)= NaN;
DepthsInterpSurface = repmat(TotalWCT,length(Cs_r),1).*repmat(Cs_r',1,size(TotalWCT,2));

%LatInterpSurface = repmat([lat_rho(1,1:end-1),lat_rho(1:end-1,end)',fliplr(lat_rho(end,2:end)),fliplr(lat_rho(2:end,1))'],length(Cs_r),1);
%LonInterpSurface = repmat([lon_rho(1,1:end-1),lon_rho(1:end-1,end)',fliplr(lon_rho(end,2:end)),fliplr(lon_rho(2:end,1))'],length(Cs_r),1);

LatInterpSurface = repmat([lat_rho(1,1:end),lat_rho(1:end,end)',lat_rho(end,1:end),lat_rho(1:end,1)'],length(Cs_r),1);
LonInterpSurface = repmat([lon_rho(1,1:end),lon_rho(1:end,end)',lon_rho(end,1:end),lon_rho(1:end,1)'],length(Cs_r),1);


%LatInterpSurface    = repmat([lat_rho(2:end,1)',lat_rho(end,:),fliplr(lat_rho(1:end-1,end)')],length(Cs_r),1);
%LonInterpSurface    = repmat([lon_rho(2:end,1)',lon_rho(end,:),fliplr(lon_rho(1:end-1,end)')],length(Cs_r),1);

LonNegIdx = find(LonInterpSurface<0.0);
LonInterpSurface(LonNegIdx)=LonInterpSurface(LonNegIdx)+360.0;

%{
Bathy = [];
for i = 1:size(theta,3),
  for j = 1:size(theta,4),

ii = find(isnan(squeeze(theta(1,:,i,j))) == 1);

if length(ii) > 0
Bathy(j,i) = depth(ii(1));
else
Bathy(j,i) = depth(end);
end
end
end

deg2rad = pi/180.0;                                                                                                          
roms_x_bathy = -(LatInterpSurface+90).*cos(LonInterpSurface*deg2rad+pi/2);                                                    
roms_y_bathy = (LatInterpSurface+90).*sin(LonInterpSurface*deg2rad+pi/2);                                                     ecco_x_bathy = -(Y+90).*cos(X*deg2rad+pi/2);                                                                                  ecco_y_bathy = (Y+90).*sin(X*deg2rad+pi/2);                                                                                       

BathyInterp = interp2(ecco_x_bathy,ecco_y_bathy,Bathy,roms_x_bathy,roms_y_bathy,'linear');

L1 = length(mask_rho(2:end,1));
L2 = length(lat_rho(end,:))+L1;
L3 = length(lat_rho(1:end-1,end))+L2;


XX = [];
YY = [];
for i = 1:50;
XX(i,:,:) = X';
YY(i,:,:) = Y';
end

ZZ = [];
for i = 1:size(Bathy,1);
    for j = 1:size(Bathy,2);
ZZ(:,j,i) = -depth;
    end
end
%}
%not a falid meshgrid! check roms coords!

%theta(theta < -10) = NaN;
%salt(salt < 0) = NaN;


isfc.tmp = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.slt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.ssh = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,2));
isfc.u =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.v =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.dpt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.prs =  nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));

[Xm,Zm,Ym] = meshgrid(LONGITUDE_T(Xloc),-DEPTH_T,LATITUDE_T(Yloc));
[Xm2,Ym2] = meshgrid(LATITUDE_T(Yloc),LONGITUDE_T(Xloc));
%%
InterpFcn = 'linear';
%adpath('/ds/projects/iomp/matlab_scripts')
disp('Interpolate to roms boundary')
for i = 1:12*(MaxYear-MinYear+1); 

%salt = m1.salt(i,:,:,:);

theta = m2.theta(i,:,:,:);
theta(theta < -10) = NaN;
In1 = squeeze(theta);
clear theta

salt = m1.salt(i,:,:,:);
salt(salt < 0) = NaN;
In2 = squeeze(salt);
clear salt
%load(fullfile(interim_dir,['cube92_iaf_ssh_',RunName,'.mat']))
ssh = m5.ssh(i,:,:);
In3 = squeeze(ssh);
clear ssh
%load(fullfile(interim_dir,['cube92_iaf_uvel_',RunName,'.mat']))
uvel = m3.uvel(i,:,:,:);
In4 = squeeze(uvel);
clear uvel
%load(fullfile(interim_dir,['cube92_iaf_vvel_',RunName,'.mat']))
vvel = m4.vvel(i,:,:,:);
In5 = squeeze(vvel);
clear vvel

isfc.tmp(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In1,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.slt(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In2,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.ssh(i,:) = inpaint_nans(interp2(Xm2,Ym2,In3,squeeze(LatInterpSurface(31,:)),squeeze(LonInterpSurface(31,:)),InterpFcn));
isfc.u(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In4,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.v(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In5,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.dpt(i,:,:) = -DepthsInterpSurface(:,:);
isfc.prs(i,:,:) = -DepthsInterpSurface(:,:);


disp(['month ' num2str(i) ' of ' num2str(12*(MaxYear-MinYear+1)) ' done.'])
save isfcMonth.dat i -ascii


end
%%
%temp1 = reshape(isfc.tmp,12,[],size(isfc.tmp,2),size(isfc.tmp,3));
%temp2 = reshape(isfc.slt,12,[],size(isfc.slt,2),size(isfc.slt,3));
%temp3 = reshape(isfc.u,12,[],size(isfc.u,2),size(isfc.u,3));
%temp4 = reshape(isfc.v,12,[],size(isfc.v,2),size(isfc.v,3));
%temp5 = reshape(isfc.dpt,12,[],size(isfc.dpt,2),size(isfc.dpt,3));
%temp6 = reshape(isfc.prs,12,[],size(isfc.prs,2),size(isfc.prs,3));
%
%isfc.tmp = squeeze(nanmean(temp1,2));
%isfc.slt = squeeze(nanmean(temp2,2));
%isfc.u = squeeze(nanmean(temp3,2));
%isfc.v = squeeze(nanmean(temp4,2));
%isfc.dpt = squeeze(nanmean(temp5,2));
%isfc.prs = squeeze(nanmean(temp6,2));

%%
save(fullfile(interim_dir,'isfc_cube92.mat'),'isfc','-v7.3')

DayPerYear=365;


%load(fullfile(interim_dir,['cube92_iaf_theta_',RunName,'.mat']))

%do_ISOM_lbc_nc_cube92(interim_dir,bryname,grdname,'Lateral Boundaries Salt and Temp',[1 1 1 1], Vtransform, Vstretching, Tcline, theta_s,theta_b,hc,N,[(DayPerYear/12/2):(DayPerYear/12):(DayPerYear*(MaxYear-MinYear+1))-(DayPerYear/12/2)],[DayPerYear*(MaxYear-MinYear+1)],'clobber')

