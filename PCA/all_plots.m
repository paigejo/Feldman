%make 6 imagesc plots of degrees between PC's for each decade pairing, one
%plot for each PC:

%set these options
useSW = 1;
useLW = 0;
lwHiRes = 0;
CLRorALL = 1; % 0 is clear-sky, 1 is all-sky
doCalculations = 0; %if 1, calculates certain variables for plots.  If the calculations are already done set to 0.

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
        iMMStr = [num2str(i), '0'];
        
        %no data from 1999, so must modify strings for first decade
        if i == 0
            iMMStr = '00';
        end
        
        fileName = ['PCA_', waveLStr, '_cellnorm_20', iMMStr, '-20', iStr, '_', resStr, '.mat'];
        
        %load data, initialize variables if necessary:
        load(fileName);
        if i == 0
            nChannels = size(V, 1);
            Vs = ones(nChannels, nPCs, nDecades)*NaN;
            lonLatScoreMats = ones(nLon, nLat, nTimeSteps, nPCs, nDecades)*NaN;
            variancesExplained = ones(nPCs, nDecades)*NaN;
            VEspace_all = ones(nLon, nLat, nPCs+1, nDecades)*NaN;
            VEtime_all = ones(nTimeSteps, nPCs+1, nDecades)*NaN;
        end
        
        %update variables with this decade's data
        Vs(:, :, i+1) = V;
        lonLatScoreMats(:, :, :, :, i+1) = lonLatScoreMat; %concatenate through time
        variancesExplained(:, i+1) = VEtotal;
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
    
    %Calculate PC1-6 core trends as function of grid cell, decade
    disp('Calculating score trend slopes')
    scoreTrendSlopes = ones(nLon, nLat, nPCs, nDecades)*NaN;
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
    
    %Calculate if PC1-6 score trend detection time as function of grid cell assuming AR(1) model
    disp('Calculating score trend detection times')
    trendDetectionTimes = ones(size(scoreTrendSlopes))*NaN;
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
    trendDetectionTimes(~isfinite(trendDetectionTimes)) = NaN; %set all infinite values to NaNs
    
    %permute lon and lat for imagesc plots
    scoreTrendSlopes = permute(scoreTrendSlopes, [2 1 3 4]);
    trendDetectionTimes = permute(trendDetectionTimes, [2 1 3 4]);
    VEspace_all = permute(VEspace_all, [2 1 3 4]);

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

%Plot PC1-6 avg score trend slopes as function of grid cell
cmap = jet(256);
for PC = 1:nPCs %PC should go up to 6
    
    avgSlopes = nanmean(scoreTrendSlopes(:, :, PC, :), 4);
    
    clim = [min(avgSlopes(:)) max(avgSlopes(:))];
    figure;
    a = imagesc(avgSlopes, clim);
    colormap(cmap);
    set(gca,'YDir','normal');
    colorbar();
    title(['Average PC ', num2str(PC) ' Score Trend Slope ', titleString]);
    set(a,'AlphaData',~isnan(avgSlopes));
    
end

%Plot PC1-6 score trend detection times as function of grid cell and
%channel
finiteTrendDetectionTimes = trendDetectionTimes(isfinite(trendDetectionTimes));
high = max(finiteTrendDetectionTimes(:))/12;
low = min(finiteTrendDetectionTimes(:))/12;
range = high - low;
clim = [low-range/253 high];
cmap = jet(256);
for PC = 1:size(finiteTrendDetectionTimes, 3) %PC should go up to 6
    
    avgTrendDetectionTimes = nanmean(finiteTrendDetectionTimes(:, :, PC, :)/12, 4);
    
    figure;
    a = imagesc(avgTrendDetectionTimes, clim);
    colormap(cmap);
    set(gca,'YDir','normal');
    colorbar();
    title(['AR(1) Average Trend Detection Time for PC ', num2str(PC) ' Scores ', titleString]);
    set(a,'AlphaData',~isnan(avgTrendDetectionTimes));
    
end

%Plot PC1-6 score trend detection times as function of grid cell and
%channel, but only if detection time is bigger than length of dataset (all
%other detection times plotted as black)
detectYears = trendDetectionTimes/12;
high = max(detectYears(:));
low = 100;
clim = [low high];
significant = detectYears <= low;
detectYears(significant) = -Inf; %significant trend cells should be black so set to black (bottom color of customJet)
cmap = customJet; % the first color of this color map is black, and colors are less dark than in jet
for PC = 1:size(detectYears, 3) %PC should go up to 6
    
    avgDetectYears = nanmean(detectYears(:, :, PC, :), 4);
    
    figure;
    a = imagesc(avgDetectYears, clim);
    colormap(cmap)
    set(gca,'YDir','normal');
    bar = colorbar();
    title(['Average AR(1) Trend Detection Time (in Years) for PC ', num2str(PC) ' Scores: Undetectable Times ', titleString]);
    set(a,'AlphaData',~isfinite(avgDetectYears));
    ytick = get(bar, 'YTick');
    
end

%Plot PC1-6 score trend detection times as function of grid cell and
%channel, but only if detection time is smaller than length of dataset (all
%other detection times plotted as black)
detectYears = trendDetectionTimes/12;
high = 100;
low = min(detectYears(:));
clim = [low high];
significant = detectYears <= high;
detectYears(~significant) = -Inf; %insignificant trend cells should be black so set to black (bottom color of customJet)
cmap = customJet; % the first color of this color map is black, and colors are less dark than in jet
for PC = 1:size(detectYears, 3) %PC should go up to 6
    
    avgDetectYears = nanmean(detectYears(:, :, PC, :), 4);
    
    figure;
    a = imagesc(avgDetectYears, clim);
    colormap(cmap)
    set(gca,'YDir','normal');
    bar = colorbar();
    title(['Average AR(1) Trend Detection Time (in Years) for PC ', num2str(PC) ' Scores: Detectable Times ', titleString]);
    set(a,'AlphaData',~isnan(avgDetectYears));
    ytick = get(bar, 'YTick');
    
end

%Variance explaned plots go here
cmap = jet(256);
for PC = 1:(nPCs+1)
    
    
    PCstr = ['PC ', num2str(PC)];
    if PC == 7
        PCstr = 'First Six PCs';
    end
    
    %spatial plot
    tmp_space = nanmean(VEspace_all(:, :, PC, :), 4);
    figure;
    a = imagesc(tmp_space);
    colormap(cmap);
    set(gca,'YDir','normal');
    colorbar();
    title(['Average Variance Explained by ', PCstr ' ', titleString]);
    set(a,'AlphaData',~isnan(tmp_space));
    
    %temporal plot
    tmp_time = nanmean(VEtime_all(:, PC, :), 3);
    figure;
    plot(tmp_time);
    title(['Average Variance Explained by ', PCstr ' ', titleString]);
    xlabel('Timestep (Months)');
    
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
    load(fileName, 'avgSPEspace', 'avgSPEspectrum', 'avgSPEtime');
    
    for PC = 1:(nPCs+1)
        
        %spectral error plots
        PCstr = ['PC ', num2str(PC)];
        if PC == 7
            PCstr = 'First Six PCs';
        end
        
        figure;
        if useSW && useLW
            plot(waveNum(1:nSW), SPEspectrum(1:nSW, i+1), 'b');
            hold on;
            plot(waveNum((nSW+1):end), SPEspectrum((nSW+1):end, i+1), 'r');
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
        imagesc(avgSPEspace(:, :, i));
        colobar();
        title(['Mean SPE For ', PCstr, ' in the ', iMMstr, '''s ', titleString]);
        
        %temporal error plots
        figure;
        plot((1:length(avgSPEtime))/12, avgSPEtime(:, PC));
        xlabel('Time (Years)');
        ylabel('SPE');
        title('SPE Through Time');
        
    end
    
end