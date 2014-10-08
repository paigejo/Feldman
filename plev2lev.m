%plev2lev(ap, b, ps, plevVar) converts an nc variable from using plev as a
%dimension to using lev as a dimension instead using scatteredInterpolant
%interpolation method.  This function assumes plevVar has the correct
%lat/lon coordinates, i.e. they match with the lat/lon coordinates of ps.
function levVar = plev2lev(ap, b, ps, plevVar, plev, goalLev)

%determine pressures values at which we want plevVar values.  Note:
%goalPressure will be 3D in lon, lat, and lev
goalPressure = hybridSigma2Pressure(ap, b, ps);

%determine scatter points at which we have plevVar values
[curLon, curLat, curPlev] = ndgrid(1:size(ps, 1), 1:size(ps, 2), plev);
curLonVec = reshape(curLon, numel(curLon), 1, 1);
curLatVec = reshape(curLat, numel(curLat), 1, 1);
curPlevVec = reshape(curPlev, numel(curPlev), 1, 1);
plevVarVec = reshape(plevVar, numel(plevVar), 1, 1);

%get rid of nans in current function values
varNan = isnan(plevVarVec);
curLonVec = curLonVec(~varNan, :);
curLatVec = curLatVec(~varNan, :);
curPlevVec = curPlevVec(~varNan, :);

%generate interpolant using Delaunay triangulation, natural interpolation,
%linear extrapolation
F = scatteredInterpolant([curLonVec, curLatVec, curPlevVec], plevVarVec, 'natural',  'linear');

%determine points at which we want plevVar interpolation values
[lon, lat, ~] = ndgrid(1:size(ps, 1), 1:size(ps, 2), 1:size(goalPressure, 3));
lonVec = reshape(lon, numel(lon), 1, 1);
latVec = reshape(lat, numel(lat), 1, 1);
goalPressureVec = reshape(goalPressure, numel(goalPressure), 1, 1);

%determine function values at interpolation points
interpVar = F([lonVec, latVec, goalPressureVec]);

%reshape computed values into the output variable
levVar = reshape(interpVar, size(plevVar, 1), size(plevVar, 2), length(goalLev));
end