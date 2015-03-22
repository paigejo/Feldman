%[RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange)
%Uses MCMC integration to estimate the spectral mean for a given set of
%MODTRAN data by randomly sampling spectra from latitude and longitude.
%Unlike sampleSpectra_decadal, sampleSpectra_decMult runs multiple
%simulations.
%
%inputs
%swSearchStr: argument to system 'ls' command to get shortwave files from.
%Only the files should be returned by ls (or their file paths). If useSW =
%0 this can be set to anything
%lwSearchStr: argument to system 'ls' command to get longwave data from
%nSamples: number of samples per simulation
%nSimulations: number of MCMC simulations
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
%RMSE: Root Mean Square Error in spectral estimate as function of sample and
%simulation
%maxPctDiff: max percent error by channel in spectral estimate as function
%of sample and simulation (in fraction units, not percent units)
%maxError: max absolute error by channel in spectral estimate as function
%of sample and simulation
%RMSEsamples: number of samples at each index
%trueMean: True spectral mean in the region specified
%sampleMeans: spectral estimates as function of sample and simulation
%waveNum: the wavenumber of longwave spectra and the wavelength of
%shortwave spectra.  If useSW = useLW = 1, then this contains first
%shortwave wavelengths then longwave wavenumbers in one vector.

function [RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, roundMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange)

useSW = logical(useSW);
useLW = logical(useLW);
lwHiRes = logical(lwHiRes);

%badDataThreshold = 1e6; %if radiance values are higher than this, set to NaN
badDataThreshold = Inf;

%%%%%%%%%%%%%%%%%%%%%
%get sw and lw files:

%sw files
[~, tmp] = system(['ls ', swSearchStr]);
swFiles = strsplit(tmp, sprintf('\n'));
swFile = swFiles{1};

%lw files
[~, tmp] = system(['ls ', lwSearchStr]);
lwFiles = strsplit(tmp, sprintf('\n'));
lwFile = lwFiles{1};

%calculate number of files
if useSW
    nTimeSteps = length(swFiles);
else
    nTimeSteps = length(lwFiles);
end

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

nSpectra = nLW + nSW;

%get lon/lat data (only exists on lw file for some reason)
lat = ncread(lwFile, 'lat');
lon = ncread(lwFile, 'lon');
nLat = length(lat);
nLon = length(lon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate min and max lon/lat indices:

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

nLonFilt = length(lonRangeI);
nLatFilt = length(latRangeI);


%%%%%%%%%%%%%%%%%%%%%%
%Generate data matrix:

%prealocate data matrices
dataMat = ones(nLonFilt, nLatFilt, nTimeSteps, nSpectra)*NaN;
data_SW = [];
data_LW = [];

%fill data matrix
for time = 1:nTimeSteps
    if mod(time, 12) == 0
        disp(['Timestep ', num2str(time), '/', num2str(nTimeSteps)])
    end
    
    %get shortwave data, remove junk data at high wavenumbers and over badDataThreshold
    if useSW
        swFile = swFiles{time};
        
        rad_low_SW_ALL = ncread(swFile, 'RADIANCE_LRES_ALL');
        %solarFlux = ncread(swFile, 'SOLAR_FLUX'); %in W/cm^2/nm
        data_SW = rad_low_SW_ALL(:, :, 1:nSW);
        data_SW(data_SW > badDataThreshold) = NaN;
        %solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    end
    
    %get longwave data, remove junk data at high wavenumbers and over badDataThreshold
    if useLW
        lwFile = lwFiles{time};
        
        if lwHiRes
            rad_LW_ALL = ncread(lwFile, 'RADIANCE_HRES_ALL');
        else
            rad_LW_ALL = ncread(lwFile, 'RADIANCE_LRES_ALL');
        end
        data_LW = rad_LW_ALL(:, :, 1:nLW);
        data_LW(data_LW > badDataThreshold) = NaN;
    end
    
    %combine shortwave, longwave data
    rads = cat(3, data_SW, data_LW);
    
    %filter data by lat/lon
    rads = rads(lonRangeI, latRangeI, :);
    
    %add rads to data matrix
    dataMat(:, :, time, :) = rads;
    
end

%convert units from W/cm^2/sr/_ to mW/m^2/sr/_ (divide by 10? No, mult by 10^7)
dataMat = dataMat*10^7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%perform sampling and analysis:

%get random sample indices
samplesPerRound = 1;
nRounds = floor(nSamples/samplesPerRound);
nSamples = samplesPerRound*nRounds;
RMSEsamples = (1:nRounds)*samplesPerRound;
sampleFiles = randi([1 nTimeSteps], nSamples, nSimulations);
sampleLons = randi([1 nLonFilt], nSamples, nSimulations);
sampleLats = randi([1 nLatFilt], nSamples, nSimulations);
sampleRows = (sampleFiles-1)*nLonFilt*nLatFilt + (sampleLats-1)*nLonFilt + sampleLons;

%reshape dataMat for easier sampling
dataMat = reshape(dataMat, [nLonFilt*nLatFilt*nTimeSteps, nSpectra]);

%calculate true mean
trueMean = myNanMean(dataMat, 1);

%do actual sampling and analysis
roundMeans = ones([nRounds, length(trueMean), nSimulations]);
for sim = 1:nSimulations
    
    %print progress
    if mod(sim, 50) == 0
        disp(['sim: ', num2str(sim)])
    end
    
    %get samples
    spectralSamples = dataMat(sampleRows(:, sim), :);
    
    %calculate running averages
    spectralCumSum = nancumsum(spectralSamples, 1);
    sampleNan = isnan(spectralSamples);
    sampleSize = cumsum(~sampleNan, 1);
    roundCumSums = spectralCumSum((1:nRounds)*samplesPerRound, :);
    roundSampleSize = sampleSize((1:nRounds)*samplesPerRound, :);
    roundMeans(:, :, sim) = roundCumSums./roundSampleSize;
    
end

%%%%%%%%%%%%%%%%%%%%
%Run error analysis:

%calculate error
error = bsxfun(@minus, roundMeans, trueMean);
maxPctDiff = squeeze(max(abs(bsxfun(@rdivide, error, trueMean)), [], 2));
RMSE = squeeze(sqrt(mean(error.^2, 2)));
maxError = squeeze(max(abs(error), [], 2));

end







