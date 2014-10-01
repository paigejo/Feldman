%This function will convert the given variable named varName from rotated
%latitude and longitude coordinates to normal latitude and longitude
%coordinates using Delaunay triangularization (linear) interpolation.  This
%function assumes the data for the given variable is already limited to a
%single time step.  Output nc file has the same name as the original, but
%with [varName, 'Unrotated'] as the file name prefix instead of just
%varName.  Variables egLat and egLon are the latitude and longitude vectors
%the user wants interpolated variable values at (interpolated coordinates
%will be cross-product of egLon with egLat).

%This function assumes the given variable is 2D (3D including time), and
%only has data for 1 timestep in the given file

%This function should be called after split_nc_files and before
%combined_nc_files if necessary to convert sic from rotated to normal
%coordinates.

function rotated_to_normal_coords(fileName, varName)

%Here's the correct latitude and longitude coordinates to interpolate at
egLat = [-88.9277353522959, -87.5387052130272, -86.1414721015279, ...
    -84.7423855907142, -83.3425960440704, -81.9424662991732, ...
    -80.5421464346171, -79.1417096486217, -77.7411958655138, ...
    -76.3406287023715, -74.9400230196494, -73.5393886337675, ...
    -72.1387322891624, -70.7380587725176, -69.3373715749609, ...
    -67.9366733025785, -66.5359659401756, -65.1352510260352, ...
    -63.7345297708429, -62.3338031405324, -60.9330719152074, ...
    -59.5323367318266, -58.1315981156439, -56.7308565037137, ...
    -55.3301122627028, -53.9293657025561, -52.5286170870997, ...
    -51.1278666423533, -49.7271145631097, -48.3263610181882, ...
    -46.9256061546646, -45.5248501013023, -44.1240929713558, ...
    -42.723334864877, -41.3225758706231, -39.9218160676465, ...
    -38.5210555266244, -37.1202943109788, -35.719532477824, ...
    -34.3187700787707, -32.918007160614, -31.5172437659226, ...
    -30.1164799335463, -28.7157156990552, -27.3149510951204, ...
    -25.9141861518467, -24.5134208970629, -23.1126553565776, ...
    -21.7118895544042, -20.3111235129604, -18.9103572532454, ...
    -17.5095907949986, -16.1088241568413, -14.7080573564048, ...
    -13.3072904104462, -11.9065233349538, -10.5057561452436, ...
    -9.10498885604852, -7.70422148160049, -6.3034540357076, ...
    -4.90268653182654, -3.5019189831313, -2.10115140257898, ...
    -0.700383802973324, 0.700383802973324, 2.10115140257898, 3.5019189831313, ...
    4.90268653182654, 6.3034540357076, 7.70422148160049, 9.10498885604852, ...
    10.5057561452436, 11.9065233349538, 13.3072904104462, 14.7080573564048, ...
    16.1088241568413, 17.5095907949986, 18.9103572532454, 20.3111235129604, ...
    21.7118895544042, 23.1126553565776, 24.5134208970629, 25.9141861518467, ...
    27.3149510951204, 28.7157156990552, 30.1164799335463, 31.5172437659226, ...
    32.918007160614, 34.3187700787707, 35.719532477824, 37.1202943109788, ...
    38.5210555266244, 39.9218160676465, 41.3225758706231, 42.723334864877, ...
    44.1240929713558, 45.5248501013023, 46.9256061546646, 48.3263610181882, ...
    49.7271145631097, 51.1278666423533, 52.5286170870997, 53.9293657025561, ...
    55.3301122627028, 56.7308565037137, 58.1315981156439, 59.5323367318266, ...
    60.9330719152074, 62.3338031405324, 63.7345297708429, 65.1352510260352, ...
    66.5359659401756, 67.9366733025785, 69.3373715749609, 70.7380587725176, ...
    72.1387322891624, 73.5393886337675, 74.9400230196494, 76.3406287023715, ...
    77.7411958655138, 79.1417096486217, 80.5421464346171, 81.9424662991732, ...
    83.3425960440704, 84.7423855907142, 86.1414721015279, 87.5387052130272, ...
    88.9277353522959];

egLon = [0, 1.40625, 2.8125, 4.21875, 5.625, 7.03125, 8.4375, 9.84375, 11.25, ...
    12.65625, 14.0625, 15.46875, 16.875, 18.28125, 19.6875, 21.09375, 22.5, ...
    23.90625, 25.3125, 26.71875, 28.125, 29.53125, 30.9375, 32.34375, 33.75, ...
    35.15625, 36.5625, 37.96875, 39.375, 40.78125, 42.1875, 43.59375, 45, ...
    46.40625, 47.8125, 49.21875, 50.625, 52.03125, 53.4375, 54.84375, 56.25, ...
    57.65625, 59.0625, 60.46875, 61.875, 63.28125, 64.6875, 66.09375, 67.5, ...
    68.90625, 70.3125, 71.71875, 73.125, 74.53125, 75.9375, 77.34375, 78.75, ...
    80.15625, 81.5625, 82.96875, 84.375, 85.78125, 87.1875, 88.59375, 90, ...
    91.40625, 92.8125, 94.21875, 95.625, 97.03125, 98.4375, 99.84375, 101.25, ...
    102.65625, 104.0625, 105.46875, 106.875, 108.28125, 109.6875, 111.09375, ...
    112.5, 113.90625, 115.3125, 116.71875, 118.125, 119.53125, 120.9375, ...
    122.34375, 123.75, 125.15625, 126.5625, 127.96875, 129.375, 130.78125, ...
    132.1875, 133.59375, 135, 136.40625, 137.8125, 139.21875, 140.625, ...
    142.03125, 143.4375, 144.84375, 146.25, 147.65625, 149.0625, 150.46875, ...
    151.875, 153.28125, 154.6875, 156.09375, 157.5, 158.90625, 160.3125, ...
    161.71875, 163.125, 164.53125, 165.9375, 167.34375, 168.75, 170.15625, ...
    171.5625, 172.96875, 174.375, 175.78125, 177.1875, 178.59375, 180, ...
    181.40625, 182.8125, 184.21875, 185.625, 187.03125, 188.4375, 189.84375, ...
    191.25, 192.65625, 194.0625, 195.46875, 196.875, 198.28125, 199.6875, ...
    201.09375, 202.5, 203.90625, 205.3125, 206.71875, 208.125, 209.53125, ...
    210.9375, 212.34375, 213.75, 215.15625, 216.5625, 217.96875, 219.375, ...
    220.78125, 222.1875, 223.59375, 225, 226.40625, 227.8125, 229.21875, ...
    230.625, 232.03125, 233.4375, 234.84375, 236.25, 237.65625, 239.0625, ...
    240.46875, 241.875, 243.28125, 244.6875, 246.09375, 247.5, 248.90625, ...
    250.3125, 251.71875, 253.125, 254.53125, 255.9375, 257.34375, 258.75, ...
    260.15625, 261.5625, 262.96875, 264.375, 265.78125, 267.1875, 268.59375, ...
    270, 271.40625, 272.8125, 274.21875, 275.625, 277.03125, 278.4375, ...
    279.84375, 281.25, 282.65625, 284.0625, 285.46875, 286.875, 288.28125, ...
    289.6875, 291.09375, 292.5, 293.90625, 295.3125, 296.71875, 298.125, ...
    299.53125, 300.9375, 302.34375, 303.75, 305.15625, 306.5625, 307.96875, ...
    309.375, 310.78125, 312.1875, 313.59375, 315, 316.40625, 317.8125, ...
    319.21875, 320.625, 322.03125, 323.4375, 324.84375, 326.25, 327.65625, ...
    329.0625, 330.46875, 331.875, 333.28125, 334.6875, 336.09375, 337.5, ...
    338.90625, 340.3125, 341.71875, 343.125, 344.53125, 345.9375, 347.34375, ...
    348.75, 350.15625, 351.5625, 352.96875, 354.375, 355.78125, 357.1875, ...
    358.59375];

%get all dimensions and variable data
lat = get_nc_variable(fileName, 'lat');
lon = get_nc_variable(fileName, 'lon');
rvar = get_nc_variable(fileName, varName);

%convert lat and lon to vector data
latVec = reshape(lat, numel(lat), 1);
lonVec = reshape(lon, numel(lon), 1);
coordVec = double([latVec, lonVec]);
varVec = reshape(rvar, numel(rvar), 1);

%create vector of interpolation coordinates
[egLatGrid, egLonGrid] = meshgrid(egLat, egLon);
egLatVec = reshape(egLatGrid, numel(egLatGrid), 1);
egLonVec = reshape(egLonGrid, numel(egLonGrid), 1);
ICoordVec = double([egLatVec, egLonVec]);

%get rid of nans
latNan = isnan(latVec);
lonNan = isnan(lonVec);
varNan = isnan(varVec);
inputNans = latNan | lonNan | varNan;
coordVec = coordVec(~inputNans, :);

%interpolate data using Delaunay triangularization
F = scatteredInterpolant(coordVec, varVec(~inputNans), 'linear',  'none'); %or maybe nearest instead of none
var = single(F(ICoordVec));

%make sure variable values are within maximum and minimum of data
varMin = single(min(varVec(~varNan)));
varMax = single(max(varVec(~varNan)));
var(var > varMax) = varMax;
var(var < varMin) = varMin;

%convert variable data to 3D matrix with (time, lat, lon) coords, singleton
%time dimension
var = reshape(var, length(egLat), length(egLon));
var = shiftdim(var, -1);

%compute output file name
breaks = strfind(fileName, '_');
outFile = [varName, 'Unrotated', fileName(breaks(1):end)];

%create new nc file with correct variable data, delete associated old data
system(['cp ', fileName, ' ', outFile]);
newVarName = [varName, 'Unrotated'];
delete_nc_variable(outFile, varName);
overwrite_nc_variable(outFile, 'lat', egLat, 'lat');
overwrite_nc_variable(outFile, 'lon', egLon, 'lon');
overwrite_nc_variable(outFile, newVarName, var, newVarName);

varsToDelete = 'rlat,rlon,lat_vertices,lon_vertices,rlat_vertices,rlon_vertices,lat_bnds,lon_bnds,rlat_bnds,rlon_bnds';
delete_nc_variable(outFile, varsToDelete);
end