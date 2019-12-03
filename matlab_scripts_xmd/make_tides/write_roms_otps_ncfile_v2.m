function write_roms_otps_ncfile_v2( Tide, grid_ncfile, output_tides_ncfile,base_date,model_file,dname)
% WRITE_ROMS_OTPS_NCFILE:
%
% PARAMETERS:
% Input:
%     Tide:
%        Structure with arrays to be written to file
%     grid_ncfile:
%        Existing grid file.  Will not be written to.
%     output_tides_ncfile:
%        To be created.



fprintf ( 1, '%s:  writing %s...\n', mfilename, output_tides_ncfile );

if exist(output_tides_ncfile)
    disp(['Deleting ' output_tides_ncfile])
    delete(output_tides_ncfile)
end

grid = ncinfo ( grid_ncfile );

dnames={grid.Dimensions.Name};
Ier=find(strcmp('eta_rho',dnames));
Ixr=find(strcmp('xi_rho',dnames));
lon_rho = ncread( grid_ncfile, 'lon_rho' );
lat_rho = ncread ( grid_ncfile, 'lat_rho' );
mask_rho = ncread ( grid_ncfile, 'mask_rho' );


disp('writing lat')
nccreate(output_tides_ncfile,'lat_rho', 'Dimensions',{'xi_rho',grid.Dimensions(Ixr).Length,'eta_rho',grid.Dimensions(Ier).Length})
ncwriteatt(output_tides_ncfile,'lat_rho','long_name','latitude of RHO-points')
ncwriteatt(output_tides_ncfile,'lat_rho','units','degree_north')
ncwriteatt(output_tides_ncfile,'lat_rho','field', 'lat_rho, scalar')
ncwrite(output_tides_ncfile,'lat_rho',lat_rho)

disp('writing lon')
nccreate(output_tides_ncfile,'lon_rho', 'Dimensions',{'xi_rho','eta_rho'})
ncwriteatt(output_tides_ncfile,'lon_rho','long_name','longitude of RHO-points')
ncwriteatt(output_tides_ncfile,'lon_rho','units','degree_east')
ncwriteatt(output_tides_ncfile,'lon_rho','field', 'lon_rho, scalar')
ncwrite(output_tides_ncfile,'lon_rho',lon_rho)


disp('writing mask')
nccreate(output_tides_ncfile,'mask_rho', 'Dimensions',{'xi_rho','eta_rho'})
ncwriteatt(output_tides_ncfile,'mask_rho','long_name','mask on RHO-points')
ncwriteatt(output_tides_ncfile,'mask_rho','option_0','land')
ncwriteatt(output_tides_ncfile,'mask_rho','option_1','water')
ncwrite(output_tides_ncfile,'mask_rho',mask_rho)


disp('TIDAL PERIOD')
nccreate(output_tides_ncfile,'tide_period', 'Dimensions',{'tide_period',inf})
ncwriteatt(output_tides_ncfile,'tide_period','long_name','tide angular period')
ncwriteatt(output_tides_ncfile,'tide_period','units','hours')
ncwriteatt(output_tides_ncfile,'tide_period','field','tide_period, scalar ')
ncwrite(output_tides_ncfile,'tide_period',Tide.period)


disp('writing Eamp')
nccreate(output_tides_ncfile,'tide_Eamp', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Eamp','long_name','tidal elevation amplitude')
ncwriteatt(output_tides_ncfile,'tide_Eamp','units','meter')
ncwriteatt(output_tides_ncfile,'tide_Eamp','field', 'tide_Eamp, scalar')
ncwrite(output_tides_ncfile,'tide_Eamp',shiftdim(Tide.Eamp,1))



disp('writing Ephase')
nccreate(output_tides_ncfile,'tide_Ephase', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Ephase','long_name','tidal elevation phase angle')
ncwriteatt(output_tides_ncfile,'tide_Ephase','units','degrees, time of maximum elevation with respect to chosen time origin')
ncwriteatt(output_tides_ncfile,'tide_Ephase','field', 'tide_Ephase, scalar ')
ncwrite(output_tides_ncfile,'tide_Ephase',shiftdim(Tide.Ephase,1))



disp('writing Cphase')
nccreate(output_tides_ncfile,'tide_Cphase', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Cphase','long_name','tidal current phase angle')
ncwriteatt(output_tides_ncfile,'tide_Cphase','units','degrees, time of maximum velocity with respect chosen time origin')
ncwriteatt(output_tides_ncfile,'tide_Cphase','field', 'tide_Cphase, scalar ')
ncwrite(output_tides_ncfile,'tide_Cphase',shiftdim(Tide.Cphase,1))


disp('writing tide_Cangle')
nccreate(output_tides_ncfile,'tide_Cangle', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Cangle','long_name','tidal current inclination angle')
ncwriteatt(output_tides_ncfile,'tide_Cangle','units','degrees between semi-major axis and East')
ncwriteatt(output_tides_ncfile,'tide_Cangle','field', 'tide_Cangle, scalar ')
ncwrite(output_tides_ncfile,'tide_Cangle',shiftdim(Tide.Cangle,1))


disp('writing tide_Cmin')
nccreate(output_tides_ncfile,'tide_Cmin', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Cmin','long_name','minimum tidal current, ellipse semi-minor axis')
ncwriteatt(output_tides_ncfile,'tide_Cmin','units','meter second-1')
ncwriteatt(output_tides_ncfile,'tide_Cmin','field', 'tide_Cmin, scalar ')
ncwrite(output_tides_ncfile,'tide_Cmin',shiftdim(Tide.Cmin,1))


disp('writing tide_Cmax')
nccreate(output_tides_ncfile,'tide_Cmax', 'Dimensions',{'xi_rho','eta_rho','tide_period'})
ncwriteatt(output_tides_ncfile,'tide_Cmax','long_name','maximum tidal current, ellipse semi-major axis')
ncwriteatt(output_tides_ncfile,'tide_Cmax','units','meter second-1')
ncwriteatt(output_tides_ncfile,'tide_Cmax','field', 'tide_Cmax, scalar ')
ncwrite(output_tides_ncfile,'tide_Cmax',shiftdim(Tide.Cmax,1))

nccreate(output_tides_ncfile,'tidal_constituents', 'Dimensions',{'tide_period','string',4},'Datatype','char')
ncwriteatt(output_tides_ncfile,'tidal_constituents','long_name','Tidal Constituent Names')
ncwrite(output_tides_ncfile,'tidal_constituents',Tide.names)

ncwriteatt(output_tides_ncfile,'/','Type','ROMS Tidal Forcing File')
ncwriteatt(output_tides_ncfile,'/','Title',['Forcing for ' dname '  domain'])
ncwriteatt(output_tides_ncfile,'/','grid_file',grid_ncfile)
ncwriteatt(output_tides_ncfile,'/','Source','OTPS-http://volkov.oce.orst.edu/tides/global.html')
ncwriteatt(output_tides_ncfile,'/','history',sprintf ( '%s:  Created by %s with %s.\n', datestr(now), getenv('USER'), mfilename ))



%
%
% clear v;
% v.Name = 'tidal_constituents';
% v.Nctype = 'char';
% v.Dimension = { 'tide_period', 'two' };
% v.Attribute(1).Name = 'long_name';
% v.Attribute(1).Value = 'Tidal Constituent Names';
%
% nc_addvar ( output_tides_ncfile, v );
%
%
% %
% % Now write the unlimited dimension data
% output_buffer.tide_period = Tide.period;
% output_buffer.tide_Ephase = Tide.Ephase;
% output_buffer.tide_Eamp = Tide.Eamp;
% output_buffer.tide_Cphase = Tide.Cphase;
% output_buffer.tide_Cangle = Tide.Cangle;
% output_buffer.tide_Cmin = Tide.Cmin;
% output_buffer.tide_Cmax = Tide.Cmax;
%
% %
% % Construct the names
% for j = 1:length(Tide.names)
% 	output_buffer.tidal_constituents(j,:) = Tide.names{j}
% end
%
%
% nc_addnewrecs ( output_tides_ncfile, output_buffer, 'tide_period' );
%
% dname=input('Enter domain namefor NetCDf attribute:');
%
% %
% % And finally add some comments
% Attribute(1).Name = 'type';
% Attribute(1).Value = 'ROMS Forcing File';
% Attribute(2).Name = 'title';
% Attribute(2).Value =['Forcing for ' dname '  domain'];
% Attribute(3).Name = 'base_date';
% Attribute(3).Value = ['days since ' datestr(base_date,31)];
% Attribute(4).Name = 'grid_file';
% Attribute(4).Value = grid_ncfile;
% Attribute(5).Name = 'source';
% Attribute(5).Value = 'OTPS';
% Attribute(6).Name = 'source_url';
% Attribute(6).Value = 'http://www.coas.oregonstate.edu/research/po/research/tide/region.html';
% Attribute(7).Name = 'history';
% Attribute(7).Value = sprintf ( '%s:  Created by %s with %s.\n', datestr(now), getenv('USER'), mfilename );
% Attribute(8).Name = 'comment';
% Attribute(8).Value = ['Inputs for OTPS executable "extract_HC" ' ...
%                     'created with m-file roms2ll.m using the grid ' ...
%                     'file as input.' model_file ' was used as the OTPS regional model.   '];
%
% %
% % Name the tidal constituents.
% const_string = Tide.names{1};
% for j = 2:length(Tide.names)
% 	const_string = sprintf ( '%s, %s', const_string, Tide.names{j} );
% end
% Attribute(9).Name = 'tidal_constituents';
% Attribute(9).Value = const_string;
%
% for j = 1:length(Attribute)
% 	nc_attput ( output_tides_ncfile, nc_global, Attribute(j).Name, Attribute(j).Value );
% end
%
%
%
% %
% % transfer the lon and lats
% lon = nc_varget ( grid_ncfile, 'lon_rho' );
% nc_varput ( output_tides_ncfile, 'lon_rho', lon );
% lat = nc_varget ( grid_ncfile, 'lat_rho' );
% nc_varput ( output_tides_ncfile, 'lat_rho', lat );
% mask_rho = nc_varget ( grid_ncfile, 'mask_rho' );
% nc_varput ( output_tides_ncfile, 'mask_rho', mask_rho );

