%This function will convert the given variable named varName from rotated
%latitude and longitude coordinates to normal latitude and longitude
%coordinates using Delaunay triangularization (linear) interpolation and
%nearest neighbor extrapolation.  The egLat and egLon vector inputs contain
%the desired unrotated lat and lon dimensions for the rotated variable to
%be interpolated to.

%This function assumes the data for the given variable is
%already limited to a single time step.  Output nc file has the same name
%as the original, but with [varName, 'Unrotated'] as the file name prefix
%instead of just varName.  Variables egLat and egLon are the latitude and
%longitude vectors the user wants interpolated variable values at
%(interpolated coordinates will be cross-product of egLon with egLat).

%This function assumes the given variable is 2D (3D including time), and
%only has data for 1 timestep in the given file

%This function should be called after split_nc_files and before
%combined_nc_files if necessary to convert sic from rotated to normal
%coordinates.

%NOTE: sometimes the longitude in the rotated lat/lon file is shifted to a
%different range than [0, 360].  In this case you should shift the lon
%variable before or after calling this function.

%NOTE: because interpolation and extrapolation is performed, there are
%values over land that should be set to 0.  This is performed in
%format_nc_files.

function rotated_to_normal_coords(fileName, varName, egLat, egLon)

%get all dimensions and variable data
lat = get_nc_variable(fileName, 'lat');
lon = get_nc_variable(fileName, 'lon'); %NOTE: sometimes lon must be shifted
rvar = get_nc_variable(fileName, varName);

%convert lat and lon to vector data
latVec = reshape(lat, numel(lat), 1);
lonVec = reshape(lon, numel(lon), 1);
coordVec = double([lonVec, latVec]);
varVec = reshape(rvar, numel(rvar), 1);

%get rid of nans
latNan = isnan(latVec);
lonNan = isnan(lonVec);
varNan = isnan(varVec);
inputNans = latNan | lonNan | varNan;
coordVec = coordVec(~inputNans, :);

%generate interpolant using Delaunay triangularization
F = scatteredInterpolant(coordVec, varVec(~inputNans), 'linear',  'nearest'); %or maybe nearest instead of none

%create vector of interpolation coordinates
[egLonGrid, egLatGrid] = ndgrid(egLon, egLat);

%interpolate at the interpolation coordinates
var = single(F(egLonGrid, egLatGrid));

%make sure variable values are within maximum and minimum of data
varMin = single(min(varVec(~varNan)));
varMax = single(max(varVec(~varNan)));
var(var > varMax) = varMax;
var(var < varMin) = varMin;

%convert variable data to 2D matrix with (lon, lat) coords, singleton
%time dimension
var = reshape(var, length(egLon), length(egLat));

%compute output file name
breaks = strfind(fileName, '_');
outFile = [varName, 'Unrotated', fileName(breaks(1):end)];

%create new nc file with correct variable data, delete associated old data
system(['cp ', fileName, ' ', outFile]);
newVarName = [varName, 'Unrotated'];
delete_nc_variable(outFile, varName);
overwrite_nc_variable(outFile, 'lat', egLat, 'lat', 1);
overwrite_nc_variable(outFile, 'lon', egLon, 'lon', 1);
overwrite_nc_variable(outFile, newVarName, var, newVarName, 3);

varsToDelete = 'i,j,i_vertices,j_vertices,i_bnds,j_bnds,rlat,rlon,lat_vertices,lon_vertices,rlat_vertices,rlon_vertices,lat_bnds,lon_bnds,rlat_bnds,rlon_bnds';
delete_nc_variable(outFile, varsToDelete);
end