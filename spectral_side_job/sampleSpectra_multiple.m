%[RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange)
%Uses MCMC integration to estimate the spectral mean for a given month of
%MODTRAN data by randomly sampling spectra from latitude and longitude.
%Unlike sampleSpectra, sampleSpectra_multiple runs multiple simulations.
%
%inputs
%swFile: the file to get shortwave data from.  If useSW = 0 this can be set
%to anything
%lwFile: file to get longwave data from (necessary no matter what useSW and
%useLW are set to, since lat and lon coordinates are only in LW files)
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

function [RMSE, maxPctDiff, maxError, RMSEsamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange)

for sim = 1:nSimulations
    
    %initialize variables
    if sim == 1
        
        [RMSE_tmp, maxPctDiff_tmp, maxError_tmp, RMSEsamples, trueMean, sampleMeans_tmp, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
        
        RMSE = ones([length(RMSE_tmp), nSimulations]);
        sampleMeans = ones([size(sampleMeans_tmp), nSimulations]);
        maxPctDiff = ones([length(maxPctDiff_tmp) nSimulations]);
        maxError = ones([length(maxError_tmp) nSimulations]);
        
    else
        
        [RMSE_tmp, maxPctDiff_tmp, maxError_tmp, ~, ~, sampleMeans_tmp, ~] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
        
    end
    
    RMSE(:, sim) = RMSE_tmp;
    sampleMeans(:, :, sim) = sampleMeans_tmp;
    maxPctDiff(:, sim) = maxPctDiff_tmp;
    maxError(:, sim) = maxError_tmp;
    
end

end







