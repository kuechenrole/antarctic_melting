addpath antarctic_mapping_tools/code/AntarcticMappingTools/
addpath antbounds/

load iceshelves_2008_v2.mat


[xgrid,ygrid] = psgrid(-90,0,5510,5,'xy');
%plot(xgrid,ygrid,'.','color',0.8*[1 1 1])
%hold on
%antbounds('gl','k')
%antbounds('coast','k')
%antbounds('shelves','k')
shelf = isiceshelf(xgrid,ygrid);
%plot(xgrid(shelf),ygrid(shelf),'kx')


for i=1:length(name)
shelf=name{i};
shelves(i).name=shelf;
[wx,wy] = antbounds_data(shelf,'xy');
shelf_poly = inpolygon(xgrid,ygrid,wx,wy);
%plot(xgrid(shelf_poly),ygrid(shelf_poly),'ro')
[shelves(i).lat,shelves(i).lon] = ps2ll(xgrid(shelf_poly),ygrid(shelf_poly));
end

save('shelves5.mat','shelves')