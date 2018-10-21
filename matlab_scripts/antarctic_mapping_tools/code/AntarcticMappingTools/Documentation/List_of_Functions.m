%% List of Functions 
% The Antarctic Mapping Tools (AMT) package for Matlab contains functions for working with Antarctic geospatial data and creating maps in Matlab. These
% functions are called by several data-specific plugins for AMT. 

%% Lookup functions 
% These functions search databases for location names: 
% 
% * <scarloc_documentation.html |scarloc|> returns geographic or polar stereographic coordinates of any location in the SCAR database.   
% * <scarloc_documentation.html |coreloc|> returns geographic coordinates of any ice core in the ITASE database. 

%% Coordinate Transformations 
% Currently, polar stereographic coordinate transformations are supported
% by AMT. Let me know if you'd like to see more.  
% 
% * <ll2ps_documentation.html |ll2ps|> transforms geo coordinates (lat,lon) to polar stereographic (x,y) meters.  
% * <ps2ll_documentation.html |ps2ll|> performs the inverse of |ll2ps|, converting polar stereographic (x,y) meters to georeferenced (lat,lon) coordinates. 
% * <psgrid_documentation.html |psgrid|> returns georeferenced coordinates of polar stereographic equally-spaced grids.
% * <uv2vxvy_documentation.html |uv2vxvy|> transforms vector components from zonal and meridional components  cartesian grid vx and vy. 
% * <vxvy2uv_documentation.html |vxvy2uv|> transforms vector components from cartesian grid vx and vy to zonal and meridional components. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/54104.html |wgs2gl04c|> converts WGS84 ellipsoid-referenced elevations to the GL04c geoid. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/54104.html |gl04c2wgs|> converts sea level-referenced elevations to WGS84 ellipsoid-referenced elevations using the GL04c geoid. 

%% Other useful calculations  
% These functions often come in handy when analyzing Antarctic geospatial
% data: 
% 
% * <inpsquad_documentation.html |inpsquad|> determines whether points are bounded by a polar stereographic quadrangle. 
% * <pathdist_documentation.html |pathdist|> calculates cumulative distance traveled along a flight line, ship track, snowmobile traverse, or satellite ground track. (Requires Mapping Toolbox) 
% * <pathdistps_documentation.html |pathdistps|> calculates cumulative distance traveled along a flight line, ship track, snowmobile traverse, or satellite ground track. (Does NOT require Mapping Toolbox)    
% * <pspath_documentation.html |pspath|> creates a path of equal spacing  in polar stereographic coordinates, such as for interpolation to a common spacing along any repeat-track analysis. 
% * <pathcrossingps71_documentation.html |pathcrossingps71|> finds intersection points of two paths.  
% * <crossovers_documentation.html |crossovers|> efficiently calculates self intersections of a spaghetti-like path. 
% * <freeboard2thickness_documentation.html |freeboard2thickness|> assumes hydrostatic equilibrium to convert ice surface elevation to ice thickness. 
% * <thickness2freeboard_documentation.html |thickness2freeboard|> assumes hydrostatic equilibrium to convert ice thickness to ice surface elevation. 

%% Mapping in polar stereographic coordinates
% These functions plot georeferenced (lat,lon) data in polar stereographic meters or kilometers. They 
% do _not_ require Matlab's Mapping Toolbox. 
% 
% * <mapzoomps_documentation.html |mapzoomps|> zooms a polar stereographic map to a location given 
% by name, geo coordinates, or polar stereographic coordinates, and can create inset maps.  
% * <plotps_documentation.html |plotps|> acts like Matlab's |plot| function, but transforms
% georeferenced data to polar stereographic cartesian coordinates before plotting. 
% * <plot3ps_documentation.html |plot3ps|> acts like Matlab's |plot3| function, but transforms
% georeferenced data to polar stereographic cartesian coordinates before plotting. 
% * <pcolorps_documentation.html |pcolorps|> acts like Matlab's |pcolor| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <surfps_documentation.html |surfps|> acts like Matlab's |surf| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <patchps_documentation.html |patchps|> acts like Matlab's |patch| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <scatterps_documentation.html |scatterps|> acts like Matlab's |scatter| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <contourps_documentation.html |contourps|> acts like Matlab's |contour| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <contourfps_documentation.html |contourfps|> acts like Matlab's |contourf| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <quiverps_documentation.html |quiverps|> is similar to Matlab's |quiver| function, but transforms zonal and meridional
% components of vector fields to polar stereographic components before plotting. 
% * <textps_documentation.html |textps|> acts like Matlab's |text| function, but transforms
% georeferenced data to polar stereographic coordinates before plotting. 
% * <graticuleps_documentation.html |graticuleps|> places a georeferenced grid or graticule on a polar stereographic 
% cartesian coordinate grid. 
% * <circleps_documentation.html |circleps|> places circles of given radii on a polar stereographic map. 
% * <geoquadps_documentation.html |geoquadps|> plots a box bound by geographic limits, on polar stereographic map. 
% * <scalebarps_documentation.html |scalebarps|> places a graphical reference scale on polar stereographic coordinates. 
% * <scarlabel_documentation.html |scarlabel|> labels features by looking up their locations in the SCAR database. 
% * <shadem_documentation.html |shadem|> applies topographic relief shading to DEM data. 
% * <coord_documentation.html |coord|> returns geo- or polar stereographic coordinates from mouse clicks.  
% * <clickz_documentation.html |clickz|> returns gridded z data such as elevation or velocity interpolated mouse click locations.    

%% Mapping in unprojected coordinates 
% Oceanographers seem to like making maps in unprojected coordinates, where longitudes 
% are plotted as x values, latitudes are y values, and south is always on the bottom. 
% Unprojected coordinates tend to have quite a bit of spatial distortion near the poles, 
% but if your data are gridded in equally-spaced lats and lons, it might make sense to 
% give each grid cell equal real estate in your plot. If you want to plot in unprojected
% coordinates, simply use Matlab's standard plotting functions, but replace x with lon and 
% y with lat. For example, 
% 
%  plot(lon,lat) 
%  pcolor(lon,lat,z) 
% 
% and so forth. AMT does not provide much support for plotting in unprojected 
% coordinates. In fact, as of now there's just one function for unprojected coordinates: 
% 
% * <inset_unproj_documentation.html |inset_unproj|> creates an inset map for plotting 
% in unprojected coordinates. 

%% Mapping with Matlab's Mapping Toolbox
% I've recently been moving away from dependence on Matlab's Mapping Toolbox because most folks
% don't have it, it's expensive, it's often computationally inefficient, and troubleshooting 
% can be a bear. Nonetheless, if you do use Matlab's Mapping Toolbox, these functions might 
% make your life a little more convenient: 
% 
% * <antmap_documentation.html |antmap|> initializes a polar stereographic
% map centered on the South Pole. The |antmap| function can also be called to place a
% grid or "graticule". 
% * <mapzoom_documentation.html |mapzoom|> calls |antmap| and zooms to any
% location given by geo coordinates or by location name. The |mapzoom|
% function can also create inset maps for zoomed regions. 
% * <scarlabel_documentation.html |scarlabel|> labels features by looking up their locations in the SCAR database. 
% * <scarlabel_documentation.html |corelabel|> labels ice core locations. 
% * <scarclick_documentation.html |scarclick|> lets you click on a map to get the names of features.  
% * <scalebar_documentation.html |scalebar|> places a graphical reference scale
% on a map. 
% * <shadem_documentation.html |shadem|> applies topographic relief shading to
% DEM data. 
% * <quivermc_documentation.html |quivermc|> shows formatted vector fields such as wind, ice, or water motion.  
% * <coord_documentation.html |coord|> returns geo- or polar stereographic coordinates from mouse clicks.  
% * <clickz_documentation.html |clickz|> returns gridded z data such as elevation or velocity interpolated mouse click locations.    

%% Helper functions 
% A couple of functions exist mostly in the background, but were developed
% at some point or another to assist an AMT plugin.  
%  
% * |islatlon| is used for input parsing to determine whether function inputs are likely georeferenced coordinates. 
% * |find2drange| returns matrix indices corresponding to a range of input coordinates. This is sometimes used to find relevant pixel indices before loading a large gridded dataset, however
% for most applications |find2drange| can now be replaced by |inpsquad|. 
% * |C2xyz| converts a Matlab-formatted contour matrix to x,y,z values. 
% 
%% Data-specific plugins 
% Here's a list of data-specific plugins available for AMT. Newer datasets may be available, so check the Matlab File Exchange
% site for the AMT tag if you want to make sure you're up to date: 
% 
% *DEMs and masks*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60246.html Antarctic boundaries and masks> from Mouginot et al. 2017 has helpful masking functions |isgrounded|, |isiceshelf|, |isopenocean|, etc. This plugin is very useful and *highly recommended*.
% * <http://www.mathworks.com/matlabcentral/fileexchange/42353.html Bedmap2 Toolbox> contains surface, bed, and ice thickness DEMs. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/56333.html CryoSat-2 Toolbox> contains functions for fast loading and interpolating CryoSat-2 surface elevation DEM data. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/52333.html IBCSO Toolbox> contains functions for plotting, interpolating, or accessing raw IBCSO DEM data. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60105.html RTopo-2> contains ice bottom, bed, surface, masks, and water column thickness from Schaffer et al. 2016. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/57205.html Bamber et al. DEM> contains Antarctic surface elevations from ERS-1 and ICESat. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/67260.html CPOM Antarctic surface elevation change> Surface elevation changes from CryoSat-2 (Shepherd et al.).
% 
% *Grounding lines, coast lines, and ice shelf outlines:*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60246.html Antarctic boundaries and masks> from Mouginot et al. 2017 has an uninterrupted grounding line and coast line. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/47640.html |asaid|> ASAID-derived break-in-slope and hydrostatic lines. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/42353.html Bedmap2 Toolbox> contains a grounding line and coast line inferred from Bedmap2 masks. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60105.html RTopo-2> contains grounding lines and coast lines from Schaffer et al. 2016. Also contains ice shelf outlines given by name. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/47329.html MEaSUREs Toolbox> contains time-dependent landward limit of flexure (ground line) from the 1990s to today. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/46611.html ICESat Grounding Zones> plots an ICESat-derived interpretation of different parts of Antarctic grounding zones. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/67286.html Grounding Line Migration> estimated from Konrad et al., 2018.
% 
% *Satellite image mosaics:*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/47282.html MODIS Mosaic of Antarctica> plots MODIS Mosaic of Antarctica images. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/49842.html |lima|> plots Landsat Image Mosaic of Antarctica images. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/52031.html |ramp|> plots Radarsat Image Mosaic of Antarctica images. 
%
% *Other gridded datasets:*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/47329.html MEaSUREs Toolbox> plots ice speed or ice velocity. Contains tools for interpolation, including interpolation along or across an arbitrary track, such as for calculating ice flux across a grounding line. Also has a function for annual velocities from 2005 to 2016. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/53152.html |flowline|> predicts ice flow lines from seed locations using MEaSUREs surface velocity data.
% * <http://www.mathworks.com/matlabcentral/fileexchange/50126.html |seaice|> plots daily sea ice concentration grids. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/52347.html SODB Toolbox> tools for plotting, interpolating, and creating profiles of Southern Ocean Database data. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/50430.html |slr_interp|> interpolates local trends in sea level rise.  
% * <http://www.mathworks.com/matlabcentral/fileexchange/47168.html |pgr_interp|> interpolates local trends of glacial isostatic adjustment or post-glacial rebound. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/56182.html |grainsize_interp|> returns snow grain sizes from MODIS MOA. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/54728.html |heatflux_interp|> returns geothermal heat flux at any Antarctic location(s). 
% * <http://www.mathworks.com/matlabcentral/fileexchange/54915.html Antarctic Gravity Data> offers raw gridded or interpolated gravity anomalies. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60565.html DTU Mean Dynamic Topography> a global 1 min dataset of dynamic topography from Andersen and Knudsen 2009.
% * <http://www.mathworks.com/matlabcentral/fileexchange/62767.html Antarctic accumulation> mean accumulation from Arthern et al. 2006. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/67497.html Antarctic Gravimetric Mass Balance> from Groh and Horwath's processing of GRACE data.
% * <http://www.mathworks.com/matlabcentral/fileexchange/67146.html Antarctic Ice Sheet basal properties> basal slipperiness and englacial rate factor from Gudmundsson. 
% 
% *Other polygons:*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/47639.html |basins|> outlines ice drainage basins from Zwally et al., 2012. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/60246.html Antarctic boundaries and masks> from Mouginot et al. 2017 contains IMBIE ice basins and IMBIE refined basins. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/52057.html |smithlakes|> plots locations of ICESat-detected active subglacial lakes.
% * <http://www.mathworks.com/matlabcentral/fileexchange/40627.html ACC Fronts> plots Antarctic Circumpolar Current Fronts. 
% 
% *Tools developed for specific types of analyses:*
% 
% * <http://www.mathworks.com/matlabcentral/fileexchange/55018.html |bedhead|> calculates subglacial pressure head or hydrostatic potential. 
% * <http://www.mathworks.com/matlabcentral/fileexchange/54264.html |iceflex_interp|> estimates tide-induced vertical flexure in the grounding zone of Antarctic ice shelves.   
% * <http://www.mathworks.com/matlabcentral/fileexchange/53766.html |ICPcampaign|> returns UTIG ICECAP field seasons based on timestamps. 
% 
% Please, feel free to develop and share your own AMT plugins.  

%% Getting started with AMT
% If you're new to AMT, new to glaciology, or new to Matlab, I recommend going through 
% the tutorial called <AMT_getting_started.html AMT Getting Started>. It will walk you 
% through some of the basics. Then go though the documentation files for each of the functions
% and make sure you understand how the work. If anything doesn't make sense or doesn't work, drop me a line
% and I'll be glad to help troubleshoot. 

%% More examples 
% Taking inspiration from musical <https://en.wikipedia.org/wiki/Fake_book
% _fake books_>, I have written a number of AMT-related example files and I'm
% calling them a <http://www.mathworks.com/matlabcentral/fileexchange/52377
% Glaciology Fake Book>. The Glaciology Fake Book is intended to 
% 
% # teach how scientific results were obtained, 
% # ensure scientific repeatability, and 
% # provide examples of how to use AMT functions. 
% 
% In addition to the Glaciology Fake Book, I've written a few other examples which rely 
% heavily on AMT. <http://www.mathworks.com/matlabcentral/fileexchange/55352-how-to-estimate-subglacial-water-routes/content/subglacial%20water%20routing/html/How_to_estimate_subglacial_water_routes.html
% How to estimate subglacial flow accumulation> uses TopoToolbox to generate subglacial water routing maps. 
% and <http://www.mathworks.com/matlabcentral/fileexchange/49693-how-to-drape-landsat-images-over-bedmap2-topography/content/html/LandsatElevationDrape.html 
% How to drape Landsat images over Bedmap2 topography> explains, well, how to drape Landsat images over Bedmap2 topography. 

%% Citing Antarctic Mapping Tools for Matlab
% If Antarctic Mapping Tools is useful for you, please cite the paper that describes AMT. If you use any 
% datasets that are available as plugins for AMT, cite those datasets too. Here's how you can cite AMT: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
% Or you may prefer BibTeX: 
% 
%  @article{amt,
%    title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%    author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%    journal={Computers \& Geosciences},
%    year={2017},
%    volume={104},
%    pages={151--157},
%    publisher={Elsevier}, 
%    doi={10.1016/j.cageo.2016.08.003}, 
%    url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
%  }
%   
%% Author Info
% Antarctic Mapping Tools, supporting plugins, and all documentation were
% written by <http://www.chadagreene.com Chad Greene> of the University of Texas at Austin's Institute
% for Geophysics with some help from David Gwyther of the University of Tasmania. 
% Some bits in AMT were adapted from codes written by Andrew Bliss, Kelly Kearney, and Andrew
% Roberts. 

