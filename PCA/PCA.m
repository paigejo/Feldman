%This script performs PCA using .nc files for a all timesteps in sw and lw
%paths focusing on the WAVELENGTH_HRES and RADIANCE_HRES_CLR variables for
%the longwave, and the WAVELENGTH_LRES, RADIANCE_LRES_CLR,
%RADIANCE_LRES_ALL, and SOLAR_FLUX variables for the shortwave.  Only the
%15 most dominant principle components are considered. PCA is performed on
%the Z-score matrix of the detrended data matrix (assuming a linear trend).
%The score matrix and the proportion of variance explained by the principle
%components are then saved to 'PCA_results.mat' in the savePath.

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

%directories with spectroscopoy files from $GRCRATCH
swPath = '/global/scratch2/sd/jpaige/PCA/osse_sw/';
lwPath = '/global/scratch2/sd/jpaige/PCA/osse_lw/';
savePath = '/global/scratch2/sd/jpaige/PCA/';

%add string and svdsecon functions to path
if isempty(strfind(path,'/global/u1/j/jpaige/git/Feldman;'))
    addpath(genpath('/global/u1/j/jpaige/git/Feldman/'));
end

% get files (for some reason copying and pasting these lines and running
% them all at once doesn't work.  Must perform strsplit operations
% separately)
cd(swPath);
[~, swFiles] = system('ls');
swFiles = strsplit(swFiles, sprintf('\n'));
cd(lwPath);
[~, lwFiles] = system('ls');
lwFiles = strsplit(lwFiles, sprintf('\n'));

%Shortwave:

%data at final several wavelength indices should be thrown out
swBuffer = 9;

%get wavenumber dimension, remove junk data at high wave numbers
cd(swPath);
waveNumLowSW = ncread(swFiles{1}, 'WAVELENGTH_LRES');
waveNumLowSW = waveNumLowSW(1:(end-swBuffer));

%Longwave:

%data at final several wavelength indices should be thrown out
lwBuffer = 4;

%get wavenumber dimension, remove junk data at high wave numbers
cd(lwPath);
waveNumHiLW = ncread(lwFiles{1}, 'WAVELENGTH_HRES');
waveNumHiLW = waveNumHiLW(1:(end-lwBuffer));

%generate data matrix:
for fid = 1:length(swFiles)
    swFile = swFiles{fid};
    lwFile = lwFiles{fid};
    
    %get shortwave data, remove junk data at high wavenumbers
    cd(swPath);
    rad_low_SW_CLR = ncread(swFile, 'RADIANCE_LRES_CLR');
    rad_low_SW_CLR = rad_low_SW_CLR(:, :, 1:length(waveNumLowSW));
    solarFlux = ncread(swFile, 'SOLAR_FLUX');
    solarFlux = solarFlux(:, :, 1:length(waveNumLowSW));
    
    %get longwave data, remove junk data at high wavenumbers
    cd(lwPath);
    rad_hi_LW_CLR = ncread(lwFile, 'RADIANCE_HRES_CLR');
    rad_hi_LW_CLR = rad_hi_LW_CLR(:, :, 1:length(waveNumHiLW));
    
    %allocate non-temporary variables, if necessary
    if fid == 1
        nrows = length(swFiles)*size(rad_low_SW_CLR, 1)*size(rad_low_SW_CLR, 2);
        ncols = size(rad_low_SW_CLR, 3) + size(rad_hi_LW_CLR, 3);
        dataMat = ones(nrows, ncols, 'single');
        
        %modify wave number dimensions to match radiance variables
        waveNumLowSWSq = shiftdim(waveNumLowSW, -2);
        waveNumHiLWSq = shiftdim(waveNumHiLW, -2);
        waveNumLowSWSq = waveNumLowSWSq.^2;
        waveNumHiLWSq = waveNumHiLWSq.^2;
    end
    
    %shortwave:
    
    %convert to radiance in meters and nanometers from radiance in centimeters
    rad_low_SW_CLR = bsxfun(@rdivide, rad_low_SW_CLR*1e-7, waveNumLowSWSq);
    
    %convert radiance to reflectance
    data_SW = rad_low_SW_CLR*pi./solarFlux;
    
    %reshape SW matrix so each row represents spectrum channels for a given
    %timestep and a given lat and lon
    data_SW = reshape(data_SW, [size(data_SW, 1)*size(data_SW, 2), size(data_SW, 3)]);
    
    %longwave:
    
    %convert to radiance in meters and micrometers from radiance in centimeters
    rad_hi_LW_CLR = bsxfun(@rdivide, rad_hi_LW_CLR*1e-4, waveNumHiLWSq);
    
    %reshape LW matrix so each row represents spectrum channels for a given
    %timestep and a given lat and lon
    data_LW = reshape(rad_hi_LW_CLR, [size(rad_hi_LW_CLR, 1)*size(rad_hi_LW_CLR, 2), size(rad_hi_LW_CLR, 3)]);
    
    %update data matrix
    numRows = size(data_LW, 1);
    startRow = (fid - 1)*numRows + 1;
    endRow = startRow + numRows - 1;
    dataMat(startRow:endRow, :) = single([data_SW, data_LW]);
end

%clear memory except dataMat, savePath, and swFiles
clearvars -except dataMat savePath swFiles

%compute zscore matrix of detrended data matrix:

%compute time column
numTimesteps = length(swFiles);
row = 1:size(dataMat, 1);
times = ceil(row/numTimesteps).';

%normalize data matrix so the average value in each column is zero and
%remove any linear trend in the columns
for col = 1:size(dataMat, 2)
    linCoeffs = polyfit(times, dataMat(:, col), 1);
    trendCol = polyval(linCoeffs, times);
    dataMat(:, col) = dataMat(:, col) - trendCol;
end

%divide each column by its standard deviation
dataMat = bsxfun(@rdivide, dataMat, std(dataMat, 0, 1));

%do PCA
numComponents = 6;
[U, S, V] = svdsecon(dataMat, numComponents);
scoreMat = U*eye(S);

% reshape score matrix to include lon, lat, time, and PC as dimensions (in that order)
cd(swPath);
lon = ncread(swFiles{1}, 'lon');
lat = ncread(swFiles{1}, 'lat');
lonLatScores = reshape(scoreMat, [length(lon), length(lat), length(swFiles), size(scoreMat, 2)]);

%Calculate proportion variance explained by each component
totalVar = norm(dataMat, 'fro')^2;
varianceExplained = S/totalVar;

cd(savePath);
save('PCA_results.mat', 'lon', 'lat', 'lonLatScores', 'varianceExplained');






