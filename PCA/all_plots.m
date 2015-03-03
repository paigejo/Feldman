%make 6 imagesc plots of degrees between PC's for each decade pairing, one
%plot for each PC:

%set these options
useSW = 1;
useLW = 1;
lwHiRes = 0;
CLRorALL = 1; % 0 is clear-sky, 1 is all-sky
doCalculations = 1; %if 1, calculates certain variables for plots.  If the calculations are already done set to 0.

%start code:

%set a few variables based on user settings
if useLW && ~useSW
    titleString = '(LW, ';
    waveLStr = 'lw';
elseif useSW && ~useLW
    titleString = '(SW, ';
    waveLStr = 'sw';
else
    titleStr = '(SW and LW, ';
    waveLStr = 'swlw';
end
if lwHiRes == 1
    resStr = 'hres';
else
    resStr = 'lres';
end
if CLRorALL
    inString = 'all_'; %either 'all_' or '' (blank for clr)
    titleString = [titleString, 'All-Sky)'];
else
    inString = '';
    titleString = [titleString, 'Clear-Sky)'];
end
nPCs = 6;
nLon = 256;
nLat = 128;
nTimeSteps = 120;
nDecades = 10;
nSW = 440;

if doCalculations
    
    for i = 0:(nDecades-1)
        
        %determine the correct file to load based on i (decade):
        %decade data specification strings
        iStr = [num2str(i), '9'];
        iMMStr = [num2str(i-1), '0'];
        
        %no data from 1999, so must modify strings for first decade
        if i == 0
            iMMStr = '00';
        end
        
        fileName = ['PCA_', waveLStr, '_final_', inString, '20', iMMStr, '-20', iStr, '.mat'];
        
        %load data, initialize variables:
        load(fileName);
        if i == 1
            nChannels = size(V, 1);
            Vs = ones(nChannels, nPCs, nDecades);
            lonLatScoreMats = ones(nLon, nLat, nTimesSteps, nPCs, nDecades);
            variancesExplained = ones(nPCs, nDecades);
            VEspace_all = ones(nLon, nLat, nPCs+1, nDecades);
            VEtime_all = ones(nTimeSteps, nPCs+1, nDecades);
        end
        
        %update variables with this decade's data
        Vs(:, :, i+1) = V;
        lonLatScoreMats(:, :, :, :, i+1) = lonLatScoreMat; %concatenate through time
        variancesExplained(:, i+1) = varianceExplained;
        VEspace_all(:, :, :, i+1) = VEspace;
        VEtime_all(:, :, i+1) = VEtime;
        
    end
    
    %calculate angles of all possible comparisons
    angles = ones(nPCs, nDecades, nDecades)*NaN;
    for i = 1:nDecades
        
        U = Vs(:, :, i);
        
        for j = 1:(i-1)
            
            V = Vs(:, :, j);
            
            for PC = 1:nPCs
                
                angles(PC, i, j) = acosd(dot(U(:, PC), V(:, PC))/sqrt(dot(U(:, PC), U(:, PC)) * dot(V(:, PC), V(:, PC))));
                
            end
        end
    end
    angles = min(angles, 180 - angles); %since we care about the linear
    %subspace spanned by the vectors, not the vectors' direction signs
    
end

clim = [0 max(angles(:))+max(angles(:))/254]; %add to max so there is slight difference between max and no data color
cmap = jet(256);
%generate imagesc plots of angles between PCs
for PC = 1:6
    
    PCAngles = squeeze(angles(PC, :, :));
    
    figure;
    a = imagesc(PCAngles, clim);
    colormap(cmap);
    colorbar();
    title(['Angles Between PC ', num2str(PC) ' of 00''s through 90''s ', titleString]);
    
    set(gca, 'XTick', 1:10);
    set(gca, 'XTickLabel', {'00''s', '10''s', '20''s', '30''s', '40''s', '50''s', '60''s', '70''s', '80''s', '90''s'})
    set(gca, 'YTick', 1:10);
    set(gca, 'YTickLabel', {'00''s', '10''s', '20''s', '30''s', '40''s', '50''s', '60''s', '70''s', '80''s', '90''s'})
    set(gca, 'TickLength', [0, 0]);
    set(a,'AlphaData',~isnan(PCAngles));
    
end

%generate imagesc plot of variances explained for each decade
figure;
imagesc(variancesExplained);
set(gca,'YDir','normal');
colorbar();
title(['Variances Explained: Principle Component Versus Decade ', titleString]);

set(gca, 'XTick', 1:nDecades);
set(gca, 'XTickLabel', {'00''s', '10''s', '20''s', '30''s', '40''s', '50''s', '60''s', '70''s', '80''s', '90''s'})
set(gca, 'YTick', 1:nPCs);
set(gca, 'TickLength', [0, 0]);

%generate error bar plot of variances explained for each decade
figure;
x = 1:nPCs;
y = mean(variancesExplained, 2);
l = y - min(variancesExplained, [], 2);
u = max(variancesExplained, [], 2) - y;
errorbar(x, y, l, u, 'b.', 'markersize', 1);
hold on
plot(x, y, 'k')
hold off
title(['Variance Explained ', titleString]);
xlabel('Principle Component');
ylabel('Fraction of Variance Explained');
set(gca, 'XTickLabel', {'', '1', '2', '3', '4', '5', '6', ''})

%generate error bar plot of cumulative variances explained for each decade
figure;
x = 1:6;
cumVarExplained = cumsum(variancesExplained, 1);
y = mean(cumVarExplained, 2);
l = y - min(cumVarExplained, [], 2);
u = max(cumVarExplained, [], 2) - y;
errorbar(x, y, l, u, 'b.', 'markersize', 1);
hold on
plot(x, y, 'k')
hold off
ylim([.6, 1]);
title(['Cumulative Variance Explained ', titleString]);
xlabel('Principle Component');
ylabel('Fraction of Variance Explained');
set(gca, 'XTickLabel', {'', '1', '2', '3', '4', '5', '6', ''})

%Calculate PC1-6 core trends as function of grid cell, decade
disp('Calculating score trend slopes')
scoreTrendSlopes = ones(nLon, nLat, nPCs, nDecades);
for dec = 1:size(lonLatScoreMats, 5)
    
    for lon = 1:size(lonLatScoreMats, 1)
        if mod(lon, 10) == 0
            disp(['Current lon is: ', num2str(lon)]);
        end
        
        for lat = 1:size(lonLatScoreMats, 2)
            
            for PC = 1:size(lonLatScoreMats, 4) %PC should go up to 6
                
                %get finite-valued score time series for this grid cell and PC
                timeSeries = squeeze(lonLatScoreMats(lon, lat, :, PC, dec));
                finite = isfinite(timeSeries);
                
                %if all data is NaN, do not calculate trend here
                if sum(finite) < 2
                    scoreTrendSlopes(lon, lat, PC) = NaN;
                    continue;
                end
                
                %enough data is finite, so calculate linear trend
                time = (1:length(timeSeries)).';
                finiteTrend = timeSeries;
                finiteTrend(~finite) = [];
                finiteTime = time;
                finiteTime(~finite) = [];
                
                %fit and save trend slope from data for this grid cell and
                %channel
                linCoeffs = polyfit(finiteTime, finiteTrend, 1);
                scoreTrendSlopes(lon, lat, PC, dec) = linCoeffs(1); %first element is linear coef, second is constant coef
                
            end
        end
    end
end
scoreTrendSlopes = permute(scoreTrendSlopes, [2 1 3 4]); %permute so that lon is x and lat is y on image plots

%Plot PC1-6 score trend slopes as function of grid cell
cmap = jet(256);
for dec = 1:nDecades
    for PC = 1:nPCs %PC should go up to 6
        
        slopesToPlot = scoreTrendSlopes(:, :, PC, dec);
        
        clim = [min(slopesToPlot(:)) max(slopesToPlot(:))];
        figure;
        a = imagesc(slopesToPlot, clim);
        colormap(cmap);
        set(gca,'YDir','normal');
        colorbar();
        title(['Slope of Trend in PC ', num2str(PC) ' Scores for Decade ', num2str(dec), ' ', titleString]);
        set(a,'AlphaData',~isnan(slopesToPlot));
        
    end
end

%Calculate if PC1-6 score trend detection time as function of grid cell assuming AR(1) model
disp('Calculating score trend detection times')
trendDetectionTimes = ones(size(scoreTrendSlopes));
for dec = 1:size(lonLatScoreMats, 5)
    
    for lon = 1:size(lonLatScoreMats, 1)
        if mod(lon, 10) == 0
            disp(['Current lon is: ', num2str(lon)]);
        end
        
        for lat = 1:size(lonLatScoreMats, 2)
            
            for PC = 1:size(lonLatScoreMats, 4) %PC should go up to 6
                
                %get finite-valued score time series for this grid cell and PC
                timeSeries = squeeze(lonLatScoreMats(lon, lat, :, PC, dec));
                finite = isfinite(timeSeries);
                
                %if too much data is NaN, do not calculate trend here
                if sum(finite) < 2
                    trendDetectionTimes(lon, lat, PC, dec) = NaN;
                    continue;
                end
                
                %otherwise, remove non-finite data
                time = (1:length(timeSeries)).';
                finiteTrend = timeSeries;
                finiteTrend(~finite) = [];
                finiteTime = time;
                finiteTime(~finite) = [];
                
                %calculate trend detection time (assuming AR(1) model)
                var_m = 0; %measurement error
                [t_detect,~] = return_detection_time_two_timeseries2( ...
                    finiteTrend,zeros(size(finiteTrend)),finiteTime,var_m);
                trendDetectionTimes(lon, lat, PC, dec) = t_detect;
                
            end
        end
    end
end
trendDetectionTimes = permute(trendDetectionTimes, [2 1 3 4]); %permute so that lon is x and lat is y on image plots

%Plot PC1-6 score trend detection times as function of grid cell and
%channel
high = log10(max(trendDetectionTimes(:))/12);
low = log10(min(trendDetectionTimes(:))/12);
range = high - low;
clim = [low-range/253 high];
cmap = jet(256);
for dec = 1:size(trendDetectionTimes, 4)
    for PC = 1:size(trendDetectionTimes, 3) %PC should go up to 6
        
        figure;
        a = imagesc(log10(trendDetectionTimes(:, :, PC, dec)/12), clim);
        colormap(cmap);
        set(gca,'YDir','normal');
        colorbar();
        title(['AR(1) Trend Detection Time (in Log10 Years) for PC ', num2str(PC) ' Scores for Decade ', num2str(dec), ' ', titleString]);
        set(a,'AlphaData',~isnan(trendDetectionTimes(:, :, PC, dec)));
        
    end
end

%Plot PC1-6 score trend detection times as function of grid cell and
%channel, but only if detection time is bigger than length of dataset (all
%other detection times plotted as black)
logDetectYears = log10(trendDetectionTimes/12);
high = max(logDetectYears(:));
low = log10(100);
clim = [low high];
significant = logDetectYears <= low;
logDetectYears(significant) = -Inf; %significant trend cells should be black so set to black (bottom color of customJet)
cmap = customJet; % the first color of this color map is black, and colors are less dark than in jet
for dec = 1:nDecades
    for PC = 1:size(logDetectYears, 3) %PC should go up to 6
        
        figure;
        a = imagesc(logDetectYears(:, :, PC, dec), clim);
        colormap(cmap)
        set(gca,'YDir','normal');
        bar = colorbar();
        title(['AR(1) Trend Detection Time (in Years) for PC ', num2str(PC) ' Scores for Decade ', num2str(dec), ': Undetectable Times ', titleString]);
        set(a,'AlphaData',~isnan(logDetectYears(:, :, PC, dec)));
        ytick = get(bar, 'YTick');
        set(bar, 'YTickLabel', 10.^ytick);
        
    end
end

%Plot PC1-6 score trend detection times as function of grid cell and
%channel, but only if detection time is smaller than length of dataset (all
%other detection times plotted as black)
logDetectYears = log10(trendDetectionTimes/12);
high = log10(100);
low = min(logDetectYears(:));
clim = [low high];
significant = logDetectYears <= high;
logDetectYears(~significant) = -Inf; %insignificant trend cells should be black so set to black (bottom color of customJet)
cmap = customJet; % the first color of this color map is black, and colors are less dark than in jet
for dec = 1:nDecades
    for PC = 1:size(logDetectYears, 3) %PC should go up to 6
        
        figure;
        a = imagesc(logDetectYears(:, :, PC, dec), clim);
        colormap(cmap)
        set(gca,'YDir','normal');
        bar = colorbar();
        title(['AR(1) Trend Detection Time (in Years) for PC ', num2str(PC) ' Scores for Decade', num2str(dec), ': Detectable Times ', titleString]);
        set(a,'AlphaData',~isnan(logDetectYears(:, :, PC, dec)));
        ytick = get(bar, 'YTick');
        set(bar, 'YTickLabel', 10.^ytick);
        
    end
end

%Variance explaned plots go here
cmap = jet(256);
for dec = 1:nDecades
    for PC = 1:(nPCs+1)
        
        
        PCstr = ['PC ', num2str(PC)];
        if PC == 7
            PCstr = 'First Six PCs';
        end
        
        %spatial plot
        tmp_space = VEspace_all(:, :, PC, dec);
        figure;
        a = imagesc(tmp_space);
        colormap(cmap);
        set(gca,'YDir','normal');
        colorbar();
        title(['Variance Explained by ', num2str(PC) ' for Decade ', num2str(dec), ' ', titleString]);
        set(a,'AlphaData',~isnan(tmp_space));
        
        %temporal plot
        tmp_time = VEtime_all(:, PC, dec);
        figure;
        plot(tmp_time);
        title(['Variance Explained by ', num2str(PC) ' for Decade ', num2str(dec), ' ', titleString]);
        xlabel('Timestep (Months)');
        
    end
end

%SPE plots go here
for i = 0:(nDecades-1)
    
    %determine the correct file to load based on i (decade):
    %decade data specification strings
    iStr = [num2str(i), '9'];
    iMMStr = [num2str(i-1), '0'];
    
    %no data from 1999, so must modify strings for first decade
    if i == 0
        iMMStr = '00';
    end
    
    fileName = ['PCA_', waveLStr, '_final_', inString, '20', iMMStr, '-20', iStr, '.mat']; 
    
    %load data:
    load(fileName, 'SPEspacetime', 'SPEspectrum');
    
    for PC = 1:(nPCs+1)
        
        %spectral error plots
        PCstr = ['PC ', num2str(PC)];
        if PC == 7
            PCstr = 'First Six PCs';
        end
        
        figure;
        avgSPEspectrum = SPEspectrum/(nLon*nLat*nTimeSteps);
        if useSW && useLW
            plot(waveNum(1:nSW), SPEspectrum(1:nSW, i+1), 'b');
            hold on;
            plot(waveNum((nSW+1):end), SPEspectrum((nSW+1):end, i+1), 'r'); %206 696 3304
            xlabel('Wavelength (nm), Wavenumber (cm^{-1})');
            hold off;
            
        elseif useSW
            plot(waveNum, SPEspectrum(:, i+1), 'b');
            xlabel('Wavelength (nm)');
            
        elseif useLW
            plot(waveNum, SPEspectrum(:, i+1), 'r');
            xlabel('Wavenumber (cm^{-1})')'
            
        end
        title(['Mean SPE For ', PCstr, ' in the ', iMMstr, '''s ', titleString]);
        
        %spatial error plots
        figure;
        avgSPEspatial = nanmean(SPEspacetime(:, :, :, PC), 3);
        imagesc(avgSPEspatial);
        colobar();
        title(['Mean SPE For ', PCstr, ' in the ', iMMstr, '''s ', titleString]);
        
    end
    
end