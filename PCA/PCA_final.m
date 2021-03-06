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

function PCA_final(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr, lwHiRes, normalize)
useSW = logical(useSW);
useLW = logical(useLW);
lwHiRes = logical(lwHiRes);
normalize = logical(normalize);

%badDataThreshold = 1e6; %if radiance values are higher than this, set to NaN
badDataThreshold = Inf;

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

waveNum = []; %SW units: nanometers.  LW units: cm^-1
nSW = 0;
nLW = 0;
if useSW
    
    %data at final several wavelength indices should be thrown out (and at
    %the fifth wavenumber?)
    swBuffer = 10;
    
    %get wavelength dimension, remove junk data at high wave lengths
    cd(swPath);
    tmp = ncread(swFiles{1}, 'WAVELENGTH_LRES');
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
    cd(lwPath);
    if lwHiRes
        tmp = ncread(lwFiles{1}, 'WAVELENGTH_HRES');
    else
        tmp = ncread(lwFiles{1}, 'WAVELENGTH_LRES');
    end
    nLW = length(tmp) - lwBuffer;
    waveNum = [waveNum; tmp(1:(end-lwBuffer))];
    
end
    

%generate data matrix:
disp('generating data matrix')

if useSW
    nTimeSteps = length(swFiles);
else
    nTimeSteps = length(lwFiles);
end
for time = 1:nTimeSteps
    if mod(time, 12) == 0
        disp(['Timestep ', num2str(time), '/', num2str(nTimeSteps)])
    end
    
    %get shortwave data, remove junk data at high wavenumbers and over badDataThreshold
    if useSW
        cd(swPath);
        swFile = swFiles{time};
        
        rad_low_SW_ALL = ncread(swFile, 'RADIANCE_LRES_ALL');
        %solarFlux = ncread(swFile, 'SOLAR_FLUX'); %in W/cm^2/nm
        rad_low_SW_ALL = rad_low_SW_ALL(:, :, 1:nSW);
        rad_low_SW_ALL(rad_low_SW_ALL > badDataThreshold) = NaN;
        %solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    end
    
    %get longwave data, remove junk data at high wavenumbers and over badDataThreshold
    if useLW
        cd(lwPath);
        lwFile = lwFiles{time};
        
        if lwHiRes
            rad_LW_ALL = ncread(lwFile, 'RADIANCE_HRES_ALL');
        else
            rad_LW_ALL = ncread(lwFile, 'RADIANCE_LRES_ALL');
        end
        rad_LW_ALL = rad_LW_ALL(:, :, 1:nLW);
        rad_LW_ALL(rad_LW_ALL > badDataThreshold) = NaN;
    end
    
    %allocate non-temporary variables and loop variables, if necessary.
    %dataMat starts out 4-dimensional ([lon lat time, channel]) so that
    %data for a single lat/lon grid cell is easily grouped (for trend
    %removal, stdev calculations)
    if time == 1
        nSpectra = nLW + nSW;
        data_SW = [];
        data_LW = [];
        
        if useSW
            nLon = size(rad_low_SW_ALL, 1);
            nLat = size(rad_low_SW_ALL, 2);
            data_SW = ones(size(rad_low_SW_ALL, 1), size(rad_low_SW_ALL, 2), nSW);
        end
        if useLW
            nLon = size(rad_LW_ALL, 1);
            nLat = size(rad_LW_ALL, 2);
            data_LW = ones(size(rad_LW_ALL, 1), size(rad_LW_ALL, 2), nLW);
        end
        
        dataMat = ones(nLon, nLat, nTime, nSpectra);
    end
    
    %shortwave:
    
    if useSW
        %convert radiance to reflectance
        %data_SW = rad_low_SW_ALL*pi*10^-6./solarFlux;
        
        %Actually don't convert to reflectance, use radiance so comparable
        %with longwave data
        data_SW = rad_low_SW_ALL;
        
    end
    
    %longwave:
    
    if useLW
        
        %convert to radiance in meters and micrometers from radiance in centimeters
        %data_LW = bsxfun(@times, rad_hi_LW_ALL*1e-4, waveNumHiLWSq); %TODO: FIX THIS
        
        %Actually no need to convert
        data_LW = rad_LW_ALL;
        
    end
    
    %update data (and add singleton dimension where time should be).  Time
    %must be third dimension and channel must be fourth so that matrix can
    %be reshaped to two-dimensional matrix with each column being a
    %channel, and each row corresponding to lat/lon/time with rows grouped
    %by timestep.
    dataMat(:, :, time, :) = permute(cat(3, data_SW, data_LW), [1 2 4 3]);
end

%clear memory except variables needed
clearvars -except dataMat useSW useLW normalize swPath lwPath savePath saveName swFiles lwFiles waveNum nTimeSteps nLon nLat nSpectra nSW nLW

%If less than 95% of the data for this grid cell and channel is finite,
%than none of it is taken into account in the calculations
disp('marking bad data')
for lon = 1:size(dataMat, 1)
    if mod(lon, 10) == 0
        disp(['Current lon is: ', num2str(lon)]);
    end
    
    for lat = 1:size(dataMat, 2)
        for channel = 1:size(dataMat, 4)
            
            %dataMat(lon, lat, :, channel) is the time series for a given
            %grid cell and channel, and non-finite numbers are not taken
            %account)
            series = dataMat(lon, lat, :, channel);
            finite = isfinite(series);
            
            if sum(finite)/length(finite) < .95
                
                %if not enough data is finite, set to NaN
                dataMat(lon, lat, :, channel) = NaN;
                
            else
                
            end
             
        end
    end
end
%set observations with a NaN to NaNs
goodRows = sum(isfinite(dataMat), 4) == nSpectra;
nanRows = ones(size(goodRows));
nanRows(~goodRows) = NaN;
dataMat = bsxfun(@times, dataMat, nanRows);

%center and nomalize matrix by cell and channel if only using SW or LW (but
%not both)
if normalize && (~useSW || ~useLW)
    disp('normalizing data matrix');
    stds = nanstd(dataMat, 0, 3);
    dataMat = bsxfun(@rdivide, dataMat, stds);
    
end
%only center by gridcell-spectral value pairings for non-combined PCA
if ~useSW || ~useLW
    disp('centering data matrix');
    cntrs = nanmean(dataMat, 3);
    dataMat = bsxfun(@minus, dataMat, cntrs);
    
end

%reshape matrix so it has dimensions [time*lat*lon, channel] in preparation for PCA
disp('reshaping and cleaning data matrix')
dataMat = reshape(dataMat, [nLon*nLat*nTimeSteps, nSpectra]);
goodRows = reshape(goodRows, [nLon*nLat*nTimeSteps, 1]);

%filter out bad data
dataMat = dataMat(goodRows, :);

%if both SW and LW data is included, normalize using RMS std for SW and
%separately for LW.  Mean-center each column of the data matrix.
if normalize && useSW && useLW
    disp('normalizing data matrix')
    
    %calculate RMS variance
    swVar = nanmean(nanvar(dataMat(:, 1:nSW), 1));
    lwVar = nanmean(nanvar(dataMat(:, (nSW+1):end), 1));
    
    %normalize
    dataMat(:, 1:nSW) = dataMat(:, 1:nSW)/sqrt(swVar);
    dataMat(:, (nSW+1):end) = dataMat(:, (nSW+1):end)/sqrt(lwVar);
    
end

if useSW && useLW
    %center the data mat
    disp('centering data matrix')
    cntrs = nanmean(dataMat, 1);
    dataMat = bsxfun(@minus, dataMat, cntrs);
    
end

%{
%find rows with non-finite values, remove them
goodRows = sum(isfinite(dataMat), 2) == nSpectra;
dataMat = dataMat(goodRows, :);

%Modify data matrix so the average value for each column is zero
disp('centering data matrix')
for channel = 1:nSpectra
    
    %calculate and subtract mean of each column
    dataMat(:, channel) = dataMat(:, channel) - mean(dataMat(:, channel));
    
    %normalize column by standard deviation
    if normalize
    	dataMat(:, channel) = dataMat(:, channel)/std(dataMat(:, channel));
    end
    
end
%}

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
%fill in holes in scoreMat with NaNs
lonLatScoreMat = NaN*ones(length(goodRows), size(scoreMat, 2));
lonLatScoreMat(goodRows, :) = scoreMat;

%reshape matrix so it has [lon lat time component] dimensions rather than
%[lon*lat*time component] dimensions (12 since first year of data removed)
lonLatScoreMat = reshape(lonLatScoreMat, [nLon, nLat, nTimeSteps, numComponents]);

%Calculate SPE for each PC individually
disp('computing SPE')
SPEspacetime = ones(nLon, nLat, nTimeSteps, numComponents+1)*NaN;
avgSPEspectrum = ones(size(dataMat, 2), numComponents+1)*NaN;
for i = 1:numComponents
    
    err_mat = dataMat - U(:, i) * S(i, i) * V(:, i)'; %mxn = mx1 x 1x1 x 1xn
    SPE = err_mat.^2;
    
    %SPEbadRows = sum(SPE, 2) > sum(dataMat.^2, 2); %sum(SPEbadRows) is 0
    %except for PC6 which has 132 bad SPE rows
    %SPEbadCols = sum(SPE, 1) > sum(dataMat.^2, 1); %sum(SPEbadCols) > 0
    
    avgSPEspectrum(:, i) = myNanMean(SPE, 1);
    tmp = NaN*ones(length(goodRows), 1);
    tmp(goodRows) = nansum(SPE, 2);
    SPEspacetime(:, :, :, i) = reshape(tmp, [nLon, nLat, nTimeSteps]);
    
end

%compute SPE, variance explained using all six PCs
err_mat = dataMat - U * S * V';
SPE = err_mat.^2;
tmp = NaN*ones(length(goodRows), 1);
nData = tmp;
tmp(goodRows) = nansum(SPE, 2);
nData(goodRows) = sum(isnan(SPE), 2);
SPEspacetime(:, :, :, numComponents+1) = reshape(tmp, [nLon, nLat, nTimeSteps]);
nDataSpacetime = reshape(nData, [nLon, nLat, nTimeSteps]);
avgSPEspace = squeeze(bsxfun(@rdivide, myNanSum(SPEspacetime, 3), nansum(nDataSpacetime, 3)));
avgSPEtime = squeeze(bsxfun(@rdivide, myNanSum(SPEspacetime, [1, 2]), myNanSum(nDataSpacetime, [1, 2])));
avgSPEspectrum(:, numComponents+1) = squeeze(nanmean(SPE, 1));

%compute variance explained in total and broken down by space, time
disp('computing variance explained')
totalVar = norm(dataMat, 'fro')^2;
VEtotal = diag(S).^2/totalVar;

squareRowSums = NaN*ones(length(goodRows), 1);
squareRowSums(goodRows) = nansum(dataMat.^2, 2);
squareRowSums = reshape(squareRowSums, [nLon, nLat, nTimeSteps]);
VEspace = 1 - squeeze(bsxfun(@rdivide, nansum(SPEspacetime, 3), nansum(squareRowSums, 3)));
VEtime = 1 - squeeze(bsxfun(@rdivide, nansum(nansum(SPEspacetime, 1), 2), nansum(nansum(squareRowSums, 1), 2)));

%save results
disp('saving results')
cd(savePath);
save(saveName, 'lonLatScoreMat', 'V', 'avgSPEspace', 'avgSPEtime', 'avgSPEspectrum', 'VEtotal', 'VEspace', 'VEtime', 'waveNum');

end

%problems with VEtime (PCs 2-6 and kinda PC 7), VEspace (PCs 2-6),
%avgSPEspace (PCs 2-6), avgSPEtime (PCs 2-6)

%avgSPEspectrum looks fine







