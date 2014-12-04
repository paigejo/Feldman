%This function performs PCA using .nc files for a all timesteps in sw and
%lw paths focusing on the WAVELENGTH_HRES and RADIANCE_HRES_CLR variables
%for the longwave, and the WAVELENGTH_LRES, RADIANCE_LRES_CLR,
%RADIANCE_LRES_ALL, and SOLAR_FLUX variables for the shortwave.  Only the 6
%most dominant principle components are considered. PCA is performed on the
%Z-score matrix of the detrended data matrix (assuming a linear trend). The
%component scores, the proportion of variance explained by the principle
%components, and the principle component matrix (each column is a component
%vector) are then saved to a file called saveName in the savePath, which is
%/global/scratch2/sd/jpaige/PCA/.  Note that the component scores are
%stored in a 4-dimensional matrix with dimensions [lon lat time
%componentNumber]. Also, this 4-dimensional matrix may contain NaNs.

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

function PCA(useSW, useLW, saveName, swPath, lwPath)
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

% get files (for some reason copying and pasting these lines and running
% them all at once doesn't work.  Must perform strsplit operations
% separately)
if useSW
    cd(swPath);
    [~, swFiles] = system('ls');
    swFiles = strsplit(swFiles, sprintf('\n'));
end
if useLW
    cd(lwPath);
    [~, lwFiles] = system('ls');
    lwFiles = strsplit(lwFiles, sprintf('\n'));
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

nextRow = 1;
for fid = 1:length(swFiles)
    swFile = swFiles{fid};
    lwFile = lwFiles{fid};
    
    %get shortwave data, remove junk data at high wavenumbers
    if useSW
        cd(swPath);
        rad_low_SW_CLR = ncread(swFile, 'RADIANCE_LRES_CLR');
        rad_low_SW_CLR = rad_low_SW_CLR(:, :, 1:length(waveNumLowSW));
        solarFlux = ncread(swFile, 'SOLAR_FLUX');
        solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    end
    
    %get longwave data, remove junk data at high wavenumbers
    if useLW
        cd(lwPath);
        rad_hi_LW_CLR = ncread(lwFile, 'RADIANCE_HRES_CLR');
        rad_hi_LW_CLR = rad_hi_LW_CLR(:, :, 1:length(waveNumHiLW));
    end
    
    %allocate non-temporary variables, if necessary
    if fid == 1
        ncols = 0;
        if useSW
            nrows = length(swFiles)*size(rad_low_SW_CLR, 1)*size(rad_low_SW_CLR, 2);
            ncols = ncols + size(rad_low_SW_CLR, 3);
        end
        if useLW
            nrows = length(swFiles)*size(rad_hi_LW_CLR, 1)*size(rad_hi_LW_CLR, 2);
            ncols = ncols + size(rad_hi_LW_CLR, 3);
        end
        
        dataMat = ones(nrows, ncols);
        goodRows = ones(nrows, 1);
        data_SW = [];
        data_LW = [];
        
        %modify wave number dimensions to match radiance variables
        if useSW
            waveNumLowSWSq = shiftdim(waveNumLowSW, -2);
            waveNumLowSWSq = waveNumLowSWSq.^2;
        end
        if useLW
            waveNumHiLWSq = shiftdim(waveNumHiLW, -2);
            waveNumHiLWSq = waveNumHiLWSq.^2;
        end
    end
    
    %shortwave:
    
    if useSW
        %convert to radiance in meters and nanometers from radiance in centimeters
        rad_low_SW_CLR = bsxfun(@rdivide, rad_low_SW_CLR*1e-7, waveNumLowSWSq);
        
        %convert radiance to reflectance
        data_SW = rad_low_SW_CLR*pi./solarFlux;
        
        %reshape SW matrix so each row represents spectrum channels for a given
        %timestep and a given lat and lon
        data_SW = reshape(data_SW, [size(data_SW, 1)*size(data_SW, 2), size(data_SW, 3)]);
        
    end
    
    %longwave:
    
    if useLW
        
        %convert to radiance in meters and micrometers from radiance in centimeters
        rad_hi_LW_CLR = bsxfun(@rdivide, rad_hi_LW_CLR*1e-4, waveNumHiLWSq);
        
        %reshape LW matrix so each row represents spectrum channels for a given
        %timestep and a given lat and lon
        data_LW = reshape(rad_hi_LW_CLR, [size(rad_hi_LW_CLR, 1)*size(rad_hi_LW_CLR, 2), size(rad_hi_LW_CLR, 3)]);
        
    end
    
    %remove non-finite data
    matChunk = [data_SW, data_LW];
    goodRowChunk = sum(~isfinite(matChunk), 2) == 0;
    matChunk = matChunk(goodRowChunk, :);
    
    %update data
    numRows = size(goodRowChunk, 1);
    startRow = (fid - 1)*numRows + 1;
    endRow = startRow + numRows - 1;
    goodRows(startRow:endRow) = goodRowChunk;
    numGoodRows = size(matChunk, 1);
    dataMat(nextRow:(nextRow + numGoodRows - 1), :) = matChunk;
    
    nextRow = nextRow + numGoodRows;
end

%clear memory except dataMat, savePath, and swFiles
clearvars -except dataMat useSW useLW swPath lwPath savePath saveName swFiles lwFiles goodRows nextRow

%trim unused rows at end of dataMat
dataMat(nextRow:end, :) = [];

%compute zscore matrix of detrended data matrix:
disp('computing zscore/detrended matrix')

%compute time column
numTimesteps = length(swFiles);
row = 1:length(goodRows);
time = ceil(row/numTimesteps).';
goodRows = logical(goodRows);
time = time(goodRows);

%normalize data matrix so the average value in each column is zero and
%remove any linear trend in the columns
for col = 1:size(dataMat, 2)
    linCoeffs = polyfit(time, dataMat(:, col), 1);
    trendCol = polyval(linCoeffs, time);
    dataMat(:, col) = dataMat(:, col) - trendCol;
end

%divide each column by its standard deviation
dataMat = bsxfun(@rdivide, dataMat, std(dataMat, 0, 1));

%do PCA, find 6 principle components, note that S is the singular values
%matrix, but contains singular values themselves not their squares
disp('perform PCA')
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
lonLatScoreMat = reshape(lonLatScoreMat, [nLon, nLat, length(swFiles), size(lonLatScoreMat, 2)]);

%save results
disp('saving results')
cd(savePath);
save(saveName, 'lonLatScoreMat', 'varianceExplained', 'V');

end




