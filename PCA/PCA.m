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

For the shortwave, use

WAVELENGTH_LRES (nm)
RADIANCE_LRES_CLR (W/cm^2/sr/nm)
RADIANCE_LRES_ALL (W/cm^2/sr/nm)
SOLAR_FLUX 


To get to reflectance = pi*radiance/solar flux
%}

%directories with spectroscopoy files from $GRCRATCH
swPath = '/global/scratch2/sd/jpaige/PCA/osse_sw/';
lwPath = '/global/scratch2/sd/jpaige/PCA/osse_lw/';
savePath = '/global/scratch2/sd/jpaige/PCA/';

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

%get wavenumber dimension, convert to wavelength in nanometers
cd(swPath);
waveNumLowSW = ncread(swFiles{1}, 'WAVELENGTH_LRES');
waveNumLowSW = waveNumLowSW(1:(end-swBuffer));
waveLLowSW = 1/(waveNumLowSW * 1e7);

%Longwave:

%data at final several wavelength indices should be thrown out
lwBuffer = 4;

%get wavenumber dimension, convert to wavelength in micrometers
cd(lwPath);
waveNumHiLW = ncread(lwFiles{1}, 'WAVELENGTH_HRES');
waveNumHiLW = waveNumHiLW(1:(end-lwBuffer));
waveLHiLW = 1/(waveNumHiLW * 1e4);

%generate data matrix:
for fid = 1:length(swFiles)
    swFile = swFiles{fid};
    lwFile = lwFiles{fid};
    
    %Shortwave:
    cd(swPath);
    
    %convert to radiance in meters and nanometers from radiance in centimeters
    rad_low_SW_CLR = ncread(swFile, 'RADIANCE_LRES_CLR');
    rad_low_SW_CLR = rad_low_SW_CLR*1e7./repmat(waveNumLowSW.^2, [1, 1, size(rad_low_SW_CLR, 3), size(rad_low_SW_CLR, 4)]);
    
    %convert radiance to reflectance
    solarFlux = ncread(swFile, 'SOLAR_FLUX');
    data_SW = rad_low_SW_CLR*pi./solarFlux;
    data_SW = data_SW(:).';
    
    %Longwave:
    cd(lwPath);
    
    %convert to radiance in meters and micrometers from radiance in centimeters
    rad_hi_LW_CLR = ncread(lwFile, 'RADIANCE_HRES_CLR');
    rad_hi_LW_CLR = rad_hi_LW_CLR*1e4./repmat(waveNumHiLW.^2, [1, 1, size(rad_hi_LW_CLR, 3), size(rad_hi_LW_CLR, 4)]);
    
    %convert radiance to reflectance
    solarFlux = ncread(lwFile, 'SOLAR_FLUX');
    data_LW = rad_hi_LW_CLR*pi./solarFlux;
    data_LW = data_LW(:).';
    
    %allocate non-temporary variables, if necessary
    if fid == 1
        dataMat = ones(length(swFiles), numel(data_SW) + numel(data_LW));
    end
    
    %update data matrix
    dataMat(fid, :) = [data_SW, data_LW];
end

%compute zscore matrix of detrended data matrix:
ZScoreMat = dataMat;

%normalize data matrix so the average value in each column is zero and
%remove any linear trend in the columns
times = 1:length(swFiles);
for col = 1:size(dataMat, 2)
    linCoeffs = polyfit(times, dataMat(:, col), 1);
    trendCol = polyval(linCoeffs, times);
    ZScoreMat(:, col) = ZScoreMat(:, col) - trendCol;
end

%divide by column standard deviations
ZScoreMat = bsxfun(@rdivide, ZScoreMat, std(ZScoreMat, 0, 1));

%do PCA
numComponents = 15;
[U, S, V] = svds(ZScoreMat, numComponents);
scoreMat = U*eye(S);

%Calculate proportion variance explained by each component
totalVar = norm(ZScoreMat, 'fro')^2;
varianceExplained = S/totalVar;

cd(savePath);
save('PCA_results.mat', 'scoreMat', 'varianceExplained');






