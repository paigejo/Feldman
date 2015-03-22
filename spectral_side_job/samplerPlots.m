function samplerPlots()

%user setup parameters:
allSky = true; %if true: all-sky.  If false: clear-sky
nSamples = 5000;
nSimulations = 500;
lwHiRes = 1;
swFile = 'b30.042a.cam2.h0.2000-01.nc.procs1024.r96.out.nc';
lwFile = 'b30.042a.cam2.h0.2000-01.nc.procs512.r110.out.nc';
savePath = '/Users/paigejo/git/Feldman/spectral_side_job';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sampleSpectra plots (not-multiple simulations):

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shortwave Tests:

useSW = 1;
useLW = 0;

%Global Test:
latRange = NaN;
lonRange = NaN;
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Global Spectral Mean ', makeTitleTag()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Global Spectral Mean ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimated');

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Arctic Spectral Mean ', makeTitleTag()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Arctic Spectral Mean ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimate');

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Tropical West Pacific Spectral Mean ', makeTitleTag()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Tropical West Pacific Spectral Mean ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimate');

%%%%%%%%%%%%%%%%
%Longwave Tests:

useSW = 0;
useLW = 1;

%Longwave Global Test, All-Sky
latRange = NaN;
lonRange = NaN;

[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, RMSE)
title(['Global Spectral Mean Estimate RMSE ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Global Spectral Mean ,' makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimate');

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, RMSE)
title(['Arctic Spectral Mean RMSE ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Arctic Spectral Mean ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimate');

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra(swFile, lwFile, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

figure;
plot(numSamples, RMSE)
title(['Tropical West Pacific Spectral Mean RMSE ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Tropical West Pacific Spectral Mean ', makeTitleTag()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, sampleMeans(end, :), 'r--');
hold off;
legend('True', 'Estimate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sampleSpectra_multiple plots (multiple simulations):

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shortwave Tests:

useSW = 1;
useLW = 0;

%Global Test:
latRange = NaN;
lonRange = NaN;
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swFile');

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Global Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(numSamples, maxError)
title(['Global Pacific Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Samples');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Global Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swFile');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Arctic Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(numSamples, maxError)
title(['Arctic Pacific Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Arctic Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swFile');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(numSamples, maxPctDiff)
title(['Max Error in Tropical West Pacific Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Percent');
xlabel('Samples');

figure;
plot(numSamples, maxError)
title(['Tropical West Pacific Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Tropical West Pacific Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/nm)');
xlabel('Wavelength (nm)');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%%%%%%%%%%%%%%%%
%Longwave Tests:

useSW = 0;
useLW = 1;

%Longwave Global Test
latRange = NaN;
lonRange = NaN;

[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwFile');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(numSamples, maxError)
title(['Global Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Global Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwFile');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(numSamples, maxError)
title(['Arctic Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Arctic Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_multiple(swFile, lwFile, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwFile');

[~, worstSim] = max(maxPctDiff(end, :));
worstMean = sampleMeans(end, :, worstSim);

figure;
plot(numSamples, maxError)
title(['Tropical West Pacific Spectral Mean RMSE ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^-1)');
xlabel('Samples');

figure;
plot(waveNum, trueMean, 'g', 'LineWidth', 1.2);
title(['Tropical West Pacific Spectral Mean ', makeTitleTag_multiple()]);
ylabel('Radiance (mW/m^2/sr/cm^{-1})');
xlabel('Wavenumber (cm^{-1})');
hold on;
plot(waveNum, worstMean, 'r--');
hold off;
legend('True', 'Worst Estimate');

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%Sub/helper functions:

    function name = makeSaveName(latRange, lonRange)
        
        name = '';
        
        if useSW
            name = [name, 'sw'];
        end
        if useLW
            name = [name, 'lw'];
        end
        name = [name, '_'];
        
        if useLW 
            if lwHiRes
                name = [name, 'hres_'];
            else
                name = [name, 'lres_'];
            end
        end
        
        if allSky
            name = [name, 'all_'];
        else
            name = [name, 'clr_'];
        end
        
        if sum(~isnan(latRange)) == 0 && sum(~isnan(lonRange)) == 0
            name = [name, 'global.mat'];
        else
            name = [name, 'lat', num2str(latRange(1)), '-', num2str(latRange(2)), '_'];
            name = [name, 'lon', num2str(lonRange(1)), '-', num2str(lonRange(2)), '.mat'];
        end
        
    end

    function tag = makeTitleTag()
        
        tag = '(';
        
        if useSW && ~useLW
            tag = [tag, 'Shortwave, '];
        elseif useLW && ~useSW
            tag = [tag, 'Longwave, '];
        elseif useSW && useLW
            tag = [tag, 'Combined, '];
        end
        
        if useLW 
            if lwHiRes
                tag = [tag, 'High-Resolution, '];
            else
                tag = [tag, 'Low-Resolution, '];
            end
        end
        
        if allSky
            tag = [tag, 'All-Sky)'];
        else
            tag = [tag, 'Clear-Sky)'];
        end
        
    end

    function tag = makeTitleTag_multiple()
        
        tag = '(';
        
        if useSW && ~useLW
            tag = [tag, 'Shortwave, '];
        elseif useLW && ~useSW
            tag = [tag, 'Longwave, '];
        elseif useSW && useLW
            tag = [tag, 'Combined, '];
        end
        
        if useLW 
            if lwHiRes
                tag = [tag, 'High-Resolution, '];
            else
                tag = [tag, 'Low-Resolution, '];
            end
        end
        
        if allSky
            tag = [tag, 'All-Sky, '];
        else
            tag = [tag, 'Clear-Sky, '];
        end
        
        tag = [tag, num2str(nSimulations), ' Simulations)'];
        
    end

end


