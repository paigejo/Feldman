%This function performs PCA on Lou or Plieades using .nc files for a all
%timesteps in sw and lw paths focusing on the WAVELENGTH_HRES and
%RADIANCE_HRES_CLR variables for the longwave, and the WAVELENGTH_LRES,
%RADIANCE_LRES_CLR, RADIANCE_LRES_ALL, and SOLAR_FLUX variables for the
%shortwave.  Only the 6 most dominant principle components are considered.
%PCA is performed on the Z-score matrix of the detrended data matrix
%(assuming a moving average trend). The component scores, the proportion of
%variance explained by the principle components, and the principle
%component matrix (each column is a component vector) are then saved to a
%file called saveName in the savePath, which is
%/global/scratch2/sd/jpaige/PCA/.  Note that the component scores are
%stored in a 4-dimensional matrix with dimensions [lon lat time
%componentNumber]. Also, this 4-dimensional matrix may contain NaNs.

%NOTE: the difference between this function and PCA.m is that a 1 year
%moving average trend is subtracted instead of a linear regression.  Hence,
%the first year of data is not used

%variables:
%{
for the long wave files, use

WAVELENGTH_HRES (200 to 2000 cm^-1, conversion to wavelength is cm^-1 = 1e4/um)

and

RADIANCE_HRES_CLR (W/cm^2/sr/cm^-1)
RADIANCE_HRES_ALL (W/cm^2/sr/cm^-1)

conversion to W/m^2/sr/um ? multiply by 10^4 to get W/m^2/sr/cm^-1 and then multiply by 1e4/wavenumber^2

OR: multiply by 10^4 to get W/m^2/sr/cm^-1 and then multiply by 1/(1e4*wavenumber^2)?

For the shortwave, use

WAVELENGTH_LRES (nm)
RADIANCE_LRES_CLR (W/cm^2/sr/nm)
RADIANCE_LRES_ALL (W/cm^2/sr/nm)
SOLAR_FLUX

To get to reflectance = pi*radiance/solar flux

NOTE: MATLAB ordering of radiance dimensions: lon, lat, wavelength
%}

function lou_PCA_moving_avg(useSW, useLW, saveName, swPath, lwPath)
useSW = logical(useSW);
useLW = logical(useLW);

%directories with spectroscopy files from $GRCRATCH
%swPath = '/global/scratch2/sd/jpaige/PCA/osse_sw/';
%lwPath = '/global/scratch2/sd/jpaige/PCA/osse_lw/';
savePath = '/global/scratch2/sd/jpaige/PCA/';

%add string and svdsecon functions to path if necessary
if isempty(strfind(path,'/global/u1/j/jpaige/git/Feldman;'))
    addpath(genpath('/global/u1/j/jpaige/git/Feldman/'));
end

% get files for each timestep (for some reason copying and pasting these
% lines and running them all at once doesn't work.  Must perform strsplit
% operations separately)
if useSW
    cd(swPath);
    [~, swFiles] = system('ls');
    swFiles = strsplit(swFiles, sprintf('\n'));
    nTime = length(swFiles);
end
if useLW
    cd(lwPath);
    [~, lwFiles] = system('ls');
    lwFiles = strsplit(lwFiles, sprintf('\n'));
    nTime = length(lwFiles);
end

%Shortwave:

if useSW
    
    %data at final several wavelength indices should be thrown out
    swBuffer = 10;
    
    %get wavenumber dimension, remove junk data at high wave numbers
    cd(swPath);
    system(['dmget ', swFiles{1}]);
    waveNumLowSW = ncread(swFiles{1}, 'WAVELENGTH_LRES');
    waveNumLowSW = waveNumLowSW(1:(end-swBuffer));
    
end

%Longwave:

if useLW
    
    %data at final several wavelength indices should be thrown out
    lwBuffer = 4;
    
    %get wavenumber dimension, remove junk data at high wave numbers
    cd(lwPath);
    system(['dmget ', lwFiles{1}]);
    waveNumHiLW = ncread(lwFiles{1}, 'WAVELENGTH_HRES');
    waveNumHiLW = waveNumHiLW(1:(end-lwBuffer));
    
end

%generate data matrix:
disp('generating data matrix')

for time = 1:length(swFiles)
    disp(['Timestep ', num2str(time), '/', num2str(length(swFiles))])
    
    %get shortwave data, remove junk data at high wavenumbers
    if useSW
        cd(swPath);
        swFile = swFiles{time};
        
        system(['dmget ', swFile]);
        rad_low_SW_CLR = ncread(swFile, 'RADIANCE_LRES_CLR');
        solarFlux = ncread(swFile, 'SOLAR_FLUX'); %in W/cm^2/nm
        rad_low_SW_CLR = rad_low_SW_CLR(:, :, 1:length(waveNumLowSW));
        solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    end
    
    %get longwave data, remove junk data at high wavenumbers
    if useLW
        cd(lwPath);
        lwFile = lwFiles{time};
        
        system(['dmget ', lwFile]);
        rad_hi_LW_CLR = ncread(lwFile, 'RADIANCE_HRES_CLR');
        rad_hi_LW_CLR = rad_hi_LW_CLR(:, :, 1:length(waveNumHiLW));
    end
    
    %allocate non-temporary variables and loop variables, if necessary.
    %dataMat starts out 4-dimensional ([lon lat time, channel]) so that
    %data for a single lat/lon grid cell is easily grouped (for trend
    %removal, stdev calculations)
    if time == 1
        ncols = 0;
        data_SW = [];
        data_LW = [];
        
        if useSW
            ncols = ncols + size(rad_low_SW_CLR, 3);
            nLon = size(rad_low_SW_CLR, 1);
            nLat = size(rad_low_SW_CLR, 2);
            data_SW = ones(size(rad_low_SW_CLR, 1), size(rad_low_SW_CLR, 2), size(rad_low_SW_CLR, 3));
        end
        if useLW
            ncols = ncols + size(rad_hi_LW_CLR, 3);
            nLon = size(rad_hi_LW_CLR, 1);
            nLat = size(rad_hi_LW_CLR, 2);
            data_LW = ones(size(rad_hi_LW_CLR, 1), size(rad_hi_LW_CLR, 2), size(rad_hi_LW_CLR, 3));
        end
        
        dataMat = ones(nLon, nLat, nTime, ncols);
        
        %modify wave number dimensions to match radiance variables
        if useLW
            waveNumHiLWSq = shiftdim(waveNumHiLW, -2);
            waveNumHiLWSq = waveNumHiLWSq.^2;
        end
    end
    
    %shortwave:
    
    if useSW
        %convert radiance to reflectance
        data_SW = rad_low_SW_CLR*pi*10^-6./solarFlux;
        
    end
    
    %longwave:
    
    if useLW
        
        %convert to radiance in meters and micrometers from radiance in centimeters
        rad_hi_LW_CLR = bsxfun(@times, rad_hi_LW_CLR*1e-4, waveNumHiLWSq); %TODO: FIX THIS
        
    end
    
    %update data (and add singleton dimension where time should be).  Time
    %must be third dimension and channel must be fourth so that matrix can
    %be reshaped to two-dimensional matrix with each column being a
    %channel, and each row corresponding to lat/lon/time with rows grouped
    %by timestep.
    dataMat(:, :, time, :) = permute(cat(3, data_SW, data_LW), [1 2 4 3]);
end

%clear memory except dataMat, savePath, and swFiles
clearvars -except dataMat useSW useLW swPath lwPath savePath saveName swFiles lwFiles

%compute zscore matrix of detrended data matrix:
disp('computing zscore/detrended matrix')

%modify data matrix so the average value for each grid cell and channel
%is zero and remove any trend (1 year moving average) in grid cell
%channel time series
for lon = 1:size(dataMat, 1)
    for lat = 1:size(dataMat, 2)
        for channel = 1:size(dataMat, 4)
            
            %calculate and subtract 1 year moving average (note that
            %dataMat(lon, lat, :, channel) is the time series for a given
            %grid cell and channel, and non-finite numbers are not taken
            %account in the moving average)
            trend = dataMat(lon, lat, :, channel);
            finite = isfinite(trend);
            finiteTrend = trend;
            finiteTrend(~finite) = 0;
            cSumVal = cumsum(finiteTrend);
            cSumNum = cumsum(finite);
            firstVal = cSumVal(1:(end-12));
            lastVal = cSumVal(13:end);
            firstNum = cSumNum(1:(end-12));
            lastNum = cSumNum(13:end);
            movingAvg = (lastVal - firstVal)./(lastNum - firstNum);
            dataMat(lon, lat, 13:end, channel) = dataMat(lon, lat, 13:end, channel) - movingAvg;
             
        end
    end
end

%clear data from startup year
dataMat = dataMat(:, :, 13:end, :);

%normalize each grid cell channel time series by its standard deviation
disp('normalizing data matrix')
dataMat = bsxfun(@rdivide, dataMat, std(dataMat, 0, 3));

%reshape matrix so it has dimensions [time*lat*lon, channel] in preparation for PCA
disp('reshaping and cleaning data matrix')
dataMat = reshape(dataMat, [size(dataMat, 1)*size(dataMat, 2)*size(dataMat, 3), size(dataMat, 4)]);

%find rows with non-finite values, remove them
goodRows = sum(isfinite(dataMat), 2) == size(dataMat, 2);
dataMat = dataMat(goodRows, :);

%do PCA, find 6 principle components, note that S is the singular values
%matrix, but contains singular values themselves not their squares
disp('performing PCA')
numComponents = 6;
[U, S, V] = svdsecon(dataMat, numComponents);
scoreMat = U*S;

%Calculate proportion variance explained by each component
totalVar = norm(dataMat, 'fro')^2;
varianceExplained = diag(S).^2/totalVar;

%get number of lon and lat values using sample variable FLNS
disp('reshaping score matrix')
if useSW
    cd(swPath);
    system(['dmget ', swFiles{1}]);
    tmp = ncread(swFiles{1}, 'FLNS');
else
    cd(lwPath);
    system(['dmget ', lwFiles{1}]);
    tmp = ncread(lwFiles{1}, 'FLNS');
end
nLon = size(tmp, 1);
nLat = size(tmp, 2);

%fill in holes in scoreMat with NaNs
lonLatScoreMat = NaN*ones(length(goodRows), size(scoreMat, 2));
lonLatScoreMat(goodRows, :) = scoreMat;

%reshape matrix so it has [lon lat time component] dimensions rather than
%[lon*lat*time component] dimensions (12 since first year of data removed)
lonLatScoreMat = reshape(lonLatScoreMat, [nLon, nLat, length(swFiles) - 12, size(lonLatScoreMat, 2)]);

%save results
disp('saving results')
cd(savePath);
save(saveName, 'lonLatScoreMat', 'varianceExplained', 'V');

end




