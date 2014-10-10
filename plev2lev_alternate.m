%plev2lev_alternate(plevVar, plev) converts an nc variable from using plev
%as a dimension to using lev as a dimension instead using
%scatteredInterpolant interpolation method.  This is different from
%plev2lev in that it makes the assumption that the plev coordinates are in
%a implicit hybrid sigma coordinates and that the max and min plev values
%correspond to the max and min lev values.  The plev values are scaled
%linearly accordingly and the variable values are then linearly
%interpolated.
function levVar = plev2lev_alternate(plevVar, plev)

%lev dimension values from b30 example file
egLev = [3.54463800000001, 7.38881350000001, 13.967214, 23.944625, ...
    37.2302900000001, 53.1146050000002, 70.0591500000003, 85.4391150000003, ...
    100.514695, 118.250335, 139.115395, 163.66207, 192.539935, 226.513265, ...
    266.481155, 313.501265000001, 368.817980000002, 433.895225000001, ...
    510.455255000002, 600.524200000003, 696.796290000003, 787.702060000003, ...
    867.160760000001, 929.648875000002, 970.554830000001, 992.5561];

%transform coordinates to hybrid sigma coordiantes corresponding to egLev
lev = plev*(max(egLev) - min(egLev))/(max(plev) - min(plev));
lev = lev + min(egLev);

%determine scatter points at which we have plevVar values
'gridding current values'
[curLon, curLat, curLev] = ndgrid(1:size(plevVar, 1), 1:size(plevVar, 2), lev);
'reshaping current value grid'
curLonVec = reshape(curLon, numel(curLon), 1);
curLatVec = reshape(curLat, numel(curLat), 1);
curLevVec = reshape(curLev, numel(curLev), 1);
plevVarVec = reshape(plevVar, numel(plevVar), 1);

%get rid of nans in current function values
varNan = isnan(plevVarVec);
curLonVec = curLonVec(~varNan, :);
curLatVec = curLatVec(~varNan, :);
curLevVec = curLevVec(~varNan, :);

%generate interpolant using Delaunay triangulation, linear interpolation,
%linear extrapolation
'computing interpolant'
F = scatteredInterpolant([curLonVec, curLatVec, curLevVec], plevVarVec, 'linear',  'linear');

%determine points at which we want plevVar interpolation values
'gridding interpolation points'
[lon, lat, egLevGrid] = ndgrid(1:size(plevVar, 1), 1:size(plevVar, 2), egLev);
'reshaping interpolation points'
lonVec = reshape(lon, numel(lon), 1);
latVec = reshape(lat, numel(lat), 1);
egLevVec = reshape(egLevGrid, numel(egLevGrid), 1);

%determine function values at interpolation points
'interpolating'
interpVar = F([lonVec, latVec, egLevVec]);

%reshape computed values into the output variable
'reshaping interpolated variable'
levVar = reshape(interpVar, size(plevVar, 1), size(plevVar, 2), length(egLev));
end