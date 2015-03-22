function samplerData_decadal(decadeRoot, allSky, lwHiRes)
%same as samplerPlots() but for decadal averages.  decadeStr should be
%something like '2000-2009', specifying the decade the data is from.

%user setup parameters:
nSamples = 1000;
nSimulations = 200;
makePlots = false; %set to true only when testing code on desktop
%swSearchStr = '/Users/paigejo/git/Feldman/PCA/test_data/nc_files/sw/*';
%lwSearchStr = '/Users/paigejo/git/Feldman/PCA/test_data/nc_files/lw/*';
%savePath = '/Users/paigejo/git/Feldman/spectral_side_job';
swSearchStr = ['/global/scratch2/sd/jpaige/PCA/sw_files/*', decadeRoot, '*.nc'];
lwSearchStr = ['/global/scratch2/sd/jpaige/PCA/lw_files/*', decadeRoot, '*.nc'];

%other setup parameters
decadeStr = [decadeRoot, '0-', decadeRoot, '9'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sampleSpectra plots (not-multiple simulations):

%this section only has plotting, no data is saved
if makePlots
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Shortwave Tests:
    
    useSW = 1;
    useLW = 0;
    
    %Global Test:
    latRange = NaN;
    lonRange = NaN;
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    [RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decadal(swSearchStr, lwSearchStr, nSamples, useSW, useLW, lwHiRes, allSky, latRange, lonRange);
    
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
    
end

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
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swSearchStr');

if makePlots
    
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
    
end

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swSearchStr');

if makePlots
    
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
    
end

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'swSearchStr');

if makePlots
    
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
    
end

%%%%%%%%%%%%%%%%
%Longwave Tests:

useSW = 0;
useLW = 1;

%Longwave Global Test
latRange = NaN;
lonRange = NaN;

[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwSearchStr');

if makePlots
    
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

end

%Arctic Test:
latRange = [70 90];
lonRange = [0 100];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwSearchStr');

if makePlots
    
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
    
end

%Tropical West Pacific Test:
latRange = [-10 10];
lonRange = [100 150];
[RMSE, maxPctDiff, maxError, numSamples, trueMean, sampleMeans, waveNum] = sampleSpectra_decMult(swSearchStr, lwSearchStr, nSamples, nSimulations, useSW, useLW, lwHiRes, allSky, latRange, lonRange);

save(makeSaveName(latRange, lonRange), 'RMSE', 'maxPctDiff', 'maxError', 'numSamples', 'trueMean', 'sampleMeans', 'waveNum', 'lwSearchStr');

if makePlots
    
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
    
end

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
            name = [name, 'global_'];
        else
            name = [name, 'lat', num2str(latRange(1)), '-', num2str(latRange(2)), '_'];
            name = [name, 'lon', num2str(lonRange(1)), '-', num2str(lonRange(2)), '_'];
        end
        
        name = [name, decadeStr, '.mat'];
        
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
        
        tag = [tag, num2str(nSimulations), ' Simulations, ', decadeStr, ')'];
        
    end

end


