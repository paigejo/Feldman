%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shortwave, All-Sky Tests:

swFile = 'b30.042a.cam2.h0.2000-01.nc.procs1024.r96.out.nc';
lwFile = 'b30.042a.cam2.h0.2000-01.nc.procs512.r110.out.nc';
lwHiRes = 0;
nSamples = 100;

%Global Test:
useSW = 1;
useLW = 0;
latRange = NaN;
lonRange = NaN;
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Global Spectral Mean (SW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/nm)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Global Spectral Mean (SW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/nm)');
xlabel('Wavelength (nm)');

%Arctic Test:
useSW = 1;
useLW = 0;
latRange = [70 90];
lonRange = [0 100];
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Arctic Spectral Mean (SW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/nm)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Arctic Spectral Mean (SW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/nm)');
xlabel('Wavelength (nm)');

%Tropical West Pacific Test:
useSW = 1;
useLW = 0;
latRange = [-10 10];
lonRange = [100 150];
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Tropical West Pacific Spectral Mean (SW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/nm)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Tropical West Pacific Spectral Mean (SW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/nm)');
xlabel('Wavelength (nm)');

%%%%%%%%%%%%%%%%
%Longwave Tests:

%Longwave Global Test, All-Sky
useSW = 0;
useLW = 1;
latRange = NaN;
lonRange = NaN;
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Global Spectral Mean (LW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Global Spectral Mean (LW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');

%Arctic Test:
useSW = 0;
useLW = 1;
latRange = [70 90];
lonRange = [0 100];
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Arctic Spectral Mean (LW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Arctic Spectral Mean (LW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');

%Tropical West Pacific Test:
useSW = 0;
useLW = 1;
latRange = [-10 10];
lonRange = [100 150];
[RMSE, RMSEsamples, trueMean, sampleMeans, waveNum, lon, lat] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, latRange, lonRange);

figure;
plot(RMSEsamples, RMSE)
title('RSE in MCMC Estimation of Tropical West Pacific Spectral Mean (LW, All-Sky)');
ylabel('Radiance RSE (W/cm^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean);
title('True Tropical West Pacific Spectral Mean (LW, All-Sky)');
ylabel('Radiance (W/cm^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');





