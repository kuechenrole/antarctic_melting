file='tidetest.nc';
fileo='espresso_tide_c05_20060101.nc';

lon=ncread(file,'lon_rho');
lat=ncread(file,'lat_rho');


zamp=ncread(file,'tide_Eamp');
zph=ncread(file,'tide_Ephase');

zamp_o=ncread(fileo,'tide_Eamp');
zph_o=ncread(fileo,'tide_Ephase');

tperiod=ncread(file,'tide_period');
tperiodo=ncread(fileo,'tide_period');

I=1;
Io=4;

tperiod(I)
tperiodo(Io)

figure(1)
colormap jet
subplot(1,2,1)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zamp(:,:,I)))
shading interp
hold off
colorbar
caxis([0 0.75])
subplot(1,2,2)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zph(:,:,I)))
shading interp
hold off
colorbar
caxis([0 50])

figure(2)
colormap jet
subplot(1,2,1)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zamp_o(:,:,Io)))
shading interp
hold off
colorbar
caxis([0 0.75])
subplot(1,2,2)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zph_o(:,:,Io)))
shading interp
hold off
colorbar
caxis([0 50])

figure(3)
colormap jet
subplot(1,2,1)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zamp_o(:,:,Io))./squeeze(zamp(:,:,I)))
shading interp
hold off
colorbar
caxis([0.9 1.1])
subplot(1,2,2)
EC_map(0)
hold on
pcolor(lon,lat,squeeze(zph_o(:,:,Io))./squeeze(zph(:,:,I)))
shading interp
hold off
colorbar
caxis([-5 5])


