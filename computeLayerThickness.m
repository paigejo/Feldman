%computes atmospheric layer thickness in meters as a function of lat, lon,
%and lev. To do this, this first computes the interlayer pressures with the
%following formula:
% pressure = ap + b * ps
%where ap = a * p0, ps is surface pressure, a and b are hybrid sigma
%coefficients, and p0 is the reference pressure.  ap and b are variables
%that should be included in the nc files of variables using plev
%coordinates and have the same length as the lev dimension. bnds variables
%have the corresponding values at the layer interfaces

%Assumes height increases proportional to decrease in log of pressure:
% z = -H * log(p/p0)

function thickness = computeLayerThickness(ap_bnds, b_bnds, ps)

H = 7000; %7 km

%get pressure values 
topPressure = hybridSigma2Pressure(ap_bnds(2, :), b_bnds(2, :), ps);
bottomPressure = hybridSigma2Pressure(ap_bnds(1, :), b_bnds(1, :), ps);
pressureRatios = topPressure./bottomPressure;

thickness = -H*log(pressureRatios);
end