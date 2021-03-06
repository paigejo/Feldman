%[RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange)
%Uses MCMC integration to estimate the spectral mean for a given month of
%MODTRAN data by randomly sampling spectra from latitude and longitude.
%
%inputs
%swFile: the file to get shortwave data from.  If useSW = 0 this can be set
%to anything
%lwFile: file to get longwave data from (necessary no matter what useSW and
%useLW are set to, since lat and lon coordinates are only in LW files)
%nSamples: number of samples in the simulation
%useSW: whether or not to use shortwave data
%useLW: whether or not to use longwave data
%lwHiRes: whether or not to use high-resolution longwave data
%allSky: if 1, use all-sky data.  Else use clear-sky
%latRange: Used for specifying the latitude range of the region whose
%spectral mean is bbeing estimated.  For global mean use latRange=NaN.
%Otherwise set latRange to a 2 element vector with the first element being
%the minimum latitude and the 2nd element being the max latitude.  Latitude
%ranges from -90N (90S) to 90N.
%lonRange: Used the same way as latRange but for specifying a longitude
%range.  Longitude ranges from 0E to 360E
%If an input variable is unnecessary, set its value to NaN
%
%outputs
%RMSE: Root Mean Square Error in spectral estimate as function of sample
%maxPctDiff: max percent error by channel in spectral estimate as function
%of sample (in fraction units, not percent units)
%maxError: max absolute error by channel in spectral estimate as function
%of sample
%RMSEsamples: number of samples at each index
%trueMean: True spectral mean in the region specified
%sampleMeans: spectral estimate through time
%waveNum: the wavenumber of longwave spectra and the wavelength of
%shortwave spectra.  If useSW = useLW = 1, then this contains first
%shortwave wavelengths then longwave wavenumbers in one vector.
function [RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange)

useSW = logical(useSW);
useLW = logical(useLW);
lwHiRes = logical(lwHiRes);

badDataThreshold = 1e6; %if radiance values are higher than this, set to NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get wavenum/wavelength, lon/lat data:

waveNum = []; %SW units: nanometers.  LW units: cm^-1
nSW = 0;
nLW = 0;

%Shortwave:
if useSW
    
    %data at final several wavelength indices should be thrown out (and at
    %the fifth wavenumber?)
    swBuffer = 10;
    
    %get wavelength dimension, remove junk data at high wave lengths
    tmp = ncread(swFile, 'WAVELENGTH_LRES');
    nSW = length(tmp) - swBuffer;
    waveNum = tmp(1:(end-swBuffer));
    
end

%Longwave:
if useLW
    
    %data at final several wavenumber indices should be thrown out
    if lwHiRes
        lwBuffer = 4;
    else
        lwBuffer = 9;
    end
    
    %get wavenumber dimension, remove junk data at high wave numbers
    if lwHiRes
        tmp = ncread(lwFile, 'WAVELENGTH_HRES');
    else
        tmp = ncread(lwFile, 'WAVELENGTH_LRES');
    end
    nLW = length(tmp) - lwBuffer;
    waveNum = [waveNum; tmp(1:(end-lwBuffer))];

end

%get lon/lat data (only exists on lw file for some reason)
lat = ncread(lwFile, 'lat');
lon = ncread(lwFile, 'lon');

%%%%%%%%%%%%%%%%%%%
%get radiance data:

%Shortwave: remove junk data at high wavenumbers and over badDataThreshold
rad_SW = [];
if useSW
    
    if allSky
        rad_SW = ncread(swFile, 'RADIANCE_LRES_ALL');
    else
        rad_SW = ncread(swFile, 'RADIANCE_LRES_CLR');
    end
    %solarFlux = ncread(swFile, 'SOLAR_FLUX'); %in W/cm^2/nm
    rad_SW = rad_SW(:, :, 1:nSW);
    rad_SW(rad_SW > badDataThreshold) = NaN;
    %solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
end

%Longwave: remove junk data at high wavenumbers and over badDataThreshold,
rad_LW = [];
if useLW
    
    if allSky
        if lwHiRes
            rad_LW = ncread(lwFile, 'RADIANCE_HRES_ALL');
        else
            rad_LW = ncread(lwFile, 'RADIANCE_LRES_ALL');
        end
    else
        if lwHiRes
            rad_LW = ncread(lwFile, 'RADIANCE_HRES_CLR');
        else
            rad_LW = ncread(lwFile, 'RADIANCE_LRES_CLR');
        end
    end
    rad_LW = rad_LW(:, :, 1:nLW);
    rad_LW(rad_LW > badDataThreshold) = NaN;
end

%allocate rads variable: 3-dimensional ([lon lat channel])
rads = cat(3, rad_SW, rad_LW);
nLon = size(rads, 1);
nLat = size(rads, 2);

%convert units from W/cm^2/sr/_ to mW/m^2/sr/_ (divide by 10? No, mult by 10^7)
rads = rads*10^7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate min and max lon/lat indices, filter all data:

%calculate lon/lat range indices
if sum(isnan(lonRange)) == 0 && sum(isnan(latRange)) == 0
    
    minLon = lonRange(1);
    maxLon = lonRange(2);
    minLat = latRange(1);
    maxLat = latRange(2);
    
    [~, minLonI] = min(abs(lon - minLon));
    [~, maxLonI] = min(abs(lon - maxLon));
    [~, minLatI] = min(abs(lat - minLat));
    [~, maxLatI] = min(abs(lat - maxLat));
    
    lonRangeI = minLonI:maxLonI;
    latRangeI = minLatI:maxLatI;
    
else
    
    lonRangeI = 1:nLon;
    latRangeI = 1:nLat;
    
end

%filter data
rads = rads(lonRangeI, latRangeI, :);

%%%%%%%%%%%%%%%%%%
%perform sampling:

trueMean = shiftdim(myNanMean(rads, [1 2]), 1);
samplesPerRound = 1;
nRounds = floor(nSamples/samplesPerRound);

sampleLons = randi([1 size(rads, 1)], samplesPerRound*nRounds, 1);
sampleLats = randi([1 size(rads, 2)], samplesPerRound*nRounds, 1);
sampleSpectra = ones(length(sampleLons), size(rads, 3));
for n = 1:length(sampleLons)
    sampleSpectra(n, :) = shiftdim(rads(sampleLons(n), sampleLats(n), :), 1);
end

%calculate averages
spectraCumSum = nancumsum(sampleSpectra, 1);
sampleNan = isnan(sampleSpectra);
sampleSize = cumsum(~sampleNan, 1);

roundCumSums = spectraCumSum((1:nRounds)*samplesPerRound, :);
roundSampleSize = sampleSize((1:nRounds)*samplesPerRound, :);

sampleMeans = roundCumSums./roundSampleSize;
error = bsxfun(@minus, sampleMeans, trueMean);
maxPctDiff = max(abs(bsxfun(@rdivide, error, trueMean)), [], 2);
RMSE = sqrt(mean(error.^2, 2));
RMSEsamples = (1:nRounds)*samplesPerRound;
maxError = max(abs(error), [], 2);
end