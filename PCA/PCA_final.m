%This function performs PCA using .nc files for a all timesteps in sw and
%lw paths focusing on the WAVELENGTH_HRES and RADIANCE_HRES_ALL variables
%for the longwave, and the WAVELENGTH_LRES, RADIANCE_LRES_ALL,
%RADIANCE_LRES_ALL, and SOLAR_FLUX variables for the shortwave.  Only the 6
%most dominant principle components are considered. PCA is performed on the
%Z-score matrix of the detrended data matrix (assuming a moving average
%trend). The component scores, the proportion of variance explained by the
%principle components, and the principle component matrix (each column is a
%component vector) are then saved to a file called saveName in the
%savePath, which is /global/scratch2/sd/jpaige/PCA/.  Note that the
%component scores are stored in a 4-dimensional matrix with dimensions [lon
%lat time componentNumber]. Also, this 4-dimensional matrix may contain
%NaNs.  The searchStr input is the search string going into the ls system
%function, and is used to determine which data files should be used in the
%analysis.

%NOTE: the difference between this function and PCA.m is that the data
%matrix is not detrended and normalized.  Also, more values are saved that
%will help with calculation of variance explained.

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

function PCA_final(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr)
useSW = logical(useSW);
useLW = logical(useLW);

%directories with spectroscopy files from $GRCRATCH
%swPath = '/global/scratch2/sd/jpaige/PCA/osse_sw/';
%lwPath = '/global/scratch2/sd/jpaige/PCA/osse_lw/';
%savePath = '/global/scratch2/sd/jpaige/PCA/';

%add string and svdsecon functions to path if necessary
if isempty(strfind(path,'/global/u1/j/jpaige/git/Feldman;'))
    addpath(genpath('/global/u1/j/jpaige/git/Feldman/'));
end

% get files for each timestep (for some reason copying and pasting these
% lines and running them all at once doesn't work.  Must perform strsplit
% operations separately)
if useSW
    cd(swPath);
    [~, swFiles] = system(['ls ', searchStr]);
    disp('Using the following shortwave files:')
    disp(swFiles)
    swFiles = strsplit(swFiles, sprintf('\n'));
    nTime = length(swFiles);
end
if useLW
    cd(lwPath);
    [~, lwFiles] = system(['ls ', searchStr]);
    disp('Using the following longwave files:')
    disp(lwFiles)
    lwFiles = strsplit(lwFiles, sprintf('\n'));
    nTime = length(lwFiles);
end

%Shortwave:

if useSW
    
    %data at final several wavelength indices should be thrown out
    swBuffer = 10;
    
    %get wavenumber dimension, remove junk data at high wave numbers
    cd(swPath);
    waveNumLowSW = ncread(swFiles{1}, 'WAVELENGTH_LRES');
    waveNumLowSW = waveNumLowSW(1:(end-swBuffer));
    
end

%Longwave:

if useLW
    
    %data at final several wavelength indices should be thrown out
    lwBuffer = 4;
    
    %get wavenumber dimension, remove junk data at high wave numbers
    cd(lwPath);
    waveNumHiLW = ncread(lwFiles{1}, 'WAVELENGTH_HRES');
    waveNumHiLW = waveNumHiLW(1:(end-lwBuffer));
    
end

%generate data matrix:
disp('generating data matrix')

if useSW
    numTimeSteps = length(swFiles);
else
    numTimeSteps = length(lwFiles);
end
for time = 1:numTimeSteps
    disp(['Timestep ', num2str(time), '/', num2str(numTimeSteps)])
    
    %get shortwave data, remove junk data at high wavenumbers
    if useSW
        cd(swPath);
        swFile = swFiles{time};
        
        rad_low_SW_ALL = ncread(swFile, 'RADIANCE_LRES_ALL');
        solarFlux = ncread(swFile, 'SOLAR_FLUX'); %in W/cm^2/nm
        rad_low_SW_ALL = rad_low_SW_ALL(:, :, 1:length(waveNumLowSW));
        solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    end
    
    %get longwave data, remove junk data at high wavenumbers
    if useLW
        cd(lwPath);
        lwFile = lwFiles{time};
        
        rad_hi_LW_ALL = ncread(lwFile, 'RADIANCE_HRES_ALL');
        rad_hi_LW_ALL = rad_hi_LW_ALL(:, :, 1:length(waveNumHiLW));
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
            ncols = ncols + size(rad_low_SW_ALL, 3);
            nLon = size(rad_low_SW_ALL, 1);
            nLat = size(rad_low_SW_ALL, 2);
            data_SW = ones(size(rad_low_SW_ALL, 1), size(rad_low_SW_ALL, 2), size(rad_low_SW_ALL, 3));
        end
        if useLW
            ncols = ncols + size(rad_hi_LW_ALL, 3);
            nLon = size(rad_hi_LW_ALL, 1);
            nLat = size(rad_hi_LW_ALL, 2);
            data_LW = ones(size(rad_hi_LW_ALL, 1), size(rad_hi_LW_ALL, 2), size(rad_hi_LW_ALL, 3));
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
        data_SW = rad_low_SW_ALL*pi*10^-6./solarFlux;
        
    end
    
    %longwave:
    
    if useLW
        
        %convert to radiance in meters and micrometers from radiance in centimeters
        %data_LW = bsxfun(@times, rad_hi_LW_ALL*1e-4, waveNumHiLWSq); %TODO: FIX THIS
        data_LW = rad_hi_LW_ALL;
        
    end
    
    %update data (and add singleton dimension where time should be).  Time
    %must be third dimension and channel must be fourth so that matrix can
    %be reshaped to two-dimensional matrix with each column being a
    %channel, and each row corresponding to lat/lon/time with rows grouped
    %by timestep.
    dataMat(:, :, time, :) = permute(cat(3, data_SW, data_LW), [1 2 4 3]);
end

%clear memory except dataMat, savePath, and swFiles
clearvars -except dataMat useSW useLW swPath lwPath savePath saveName swFiles lwFiles numTimeSteps

%compute zscore matrix of detrended data matrix:
disp('centering data matrix')

%Modify data matrix so the average value for each grid cell and channel is
%zero. Note that if less than 95% of the data for this grid cell and
%channel is finite, than none of it is taken into account in the
%calculations
for lon = 1:size(dataMat, 1)
    if mod(lon, 10) == 0
        disp(['Current lon is: ', num2str(lon)]);
    end
    
    for lat = 1:size(dataMat, 2)
        for channel = 1:size(dataMat, 4)
            
            %calculate and subtract 1 year moving average (note that
            %dataMat(lon, lat, :, channel) is the time series for a given
            %grid cell and channel, and non-finite numbers are not taken
            %account in the moving average)
            series = dataMat(lon, lat, :, channel);
            finite = isfinite(series);
            
            if sum(finite)/length(finite) < .95
                
                %if not enough data is finite, set to NaN
                dataMat(lon, lat, :, channel) = NaN;
                
            else
                
                %enough data is finite, so calculate and subtract mean
                cntr = mean(series(finite));
                dataMat(lon, lat, :, channel) = dataMat(lon, lat, :, channel) - cntr;
                
            end
             
        end
    end
end

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

%compute score matrix
disp('computing PCA scores');
scoreMat = U*S;

%get number of lon and lat values using sample variable FLNS
disp('reshaping score matrix')
if useSW
    cd(swPath);
    tmp = ncread(swFiles{1}, 'FLNS');
else
    cd(lwPath);
    tmp = ncread(lwFiles{1}, 'FLNS');
end
nLon = size(tmp, 1);
nLat = size(tmp, 2);

%fill in holes in scoreMat with NaNs
lonLatScoreMat = NaN*ones(length(goodRows), size(scoreMat, 2));
lonLatScoreMat(goodRows, :) = scoreMat;

%reshape matrix so it has [lon lat time component] dimensions rather than
%[lon*lat*time component] dimensions (12 since first year of data removed)
lonLatScoreMat = reshape(lonLatScoreMat, [nLon, nLat, numTimeSteps, numComponents]);

%Calculate SPE for each PC individually
disp('computing SPE')
SPEspacetime = ones(nLon, nLat, numTimeSteps, numComponents+1);
SPEspectrum = ones(size(dataMat, 2), numComponents+1);
for i = 1:numComponents
    
    err_mat = dataMat - U(:, i) * S(i, i) * V(:, i)'; %mxn = mx1 x 1x1 x 1xn
    SPE = err_mat.^2; %CHECK HERE IF ANY ELEMENT OF SPE IS HIGHER THAN dataMatSQ (.06 percent of PC1 predicted vals are greater than actual vals)
    
    SPEbadRows = sum(SPE, 2) > sum(dataMat.^2, 2);
    SPEbadCols = sum(SPE, 1) > sum(dataMat.^2, 1);
    
    SPEspectrum(:, i) = sum(SPE, 1);
    tmp = sum(SPE, 2);
    SPEspacetime = reshape(tmp, [nLon, nLat, numTimeSteps]);
    
end

%compute SPE, variance explained using all six PCs
err_mat = dataMat - U * S * V';
SPE = err_mat.^2;
SPEspectrum(:, numComponents+1) = sum(SPE, 1);
tmp = sum(SPE, 2);
SPEspacetime(:, :, :, numComponents+1) = reshape(tmp, [nLon, nLat, numTimeSteps]);

%compute variance explained in total and broken down by space, time (and expectrum?)
disp('computing variance explained')
totalVar = norm(dataMat, 'fro')^2;
VEtotal = diag(S).^2/totalVar;

dataMatSQ = dataMat.^2; %TODO: figure out how to rename dataMat without having to define new variable!!!!!
dataMatSQ = reshape(dataMatSQ, [nLon, nLat, numTimeSteps, size(dataMat, 2)]);
VEspace = bsxfun(@rdivide, sum(SPEspacetime, 3), sum(sum(dataMatSQ, 3), 4));
VEtime = bsxfun(@rdivide, sum(sum(SPEspacetime, 1), 2), sum(sum(sum(dataMatSQ, 1), 2), 4));
VEspectrum = bsxfun(@rdivide, SPEspectrum, squeeze(sum(sum(sum(dataMatSQ, 1), 2), 3)));

%save results
disp('saving results')
cd(savePath);
save(saveName, 'lonLatScoreMat', 'V', 'SPEspacetime', 'SPEspectrum', 'VEtotal', 'VEspace', 'VEtime', 'VEspectrum');

end









