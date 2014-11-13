%formula given in b30 file: pressure = a(k) * p0 + b(k) * ps
%pressure(lon, lat, lev, time) = ap(lev) + b(lev) * ps(lon, lat, time) in frac units
%plev(plev) in Pa

%This function converts from lev coordinates, which are in a hybrid sigma
%pressure coordinate system, to pressure coordinates in Pascals.  Required
%are the a and b components of the hybrid sigma coordinates (a * p0, and
%b * ps).

function pressure = hybridSigma2Pressure(ap, b, ps)


pressure = repmat(ps, [1 1 length(ap)]);
for lev = 1:length(ap)
    pressure(:, :, lev) = pressure(:, :, lev) * b(lev) + ap(lev);
end
end