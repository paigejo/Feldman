%this script performs all moving average PCA analysis on Lou (for all
%decades)

swPath = '/lou/s2m/dfeldman/osse_sw/';
lwPath = '/lou/s2m/dfeldman/osse_lw/';

for i = 0:9
    
    if i == 1
        continue; %since I've already calculated for this decade
    end
    
    %decade data specification strings
    iStr = [num2str(i), '9'];
    iMMStr = [num2str(i-1), '9'];
    
    %no data from 1999, so must modify strings for first decade
    if i == 0
        iMMStr = '00';
    end
    
    disp(['Running PCA code on 20', iMMStr, '-20', iStr, ' decade']);
    
    %determine inputs for PCA function based on decade
    savePath = ['/lou/s2j/jlpaige/osse_sw/b30.036a.cam2/20', iMMStr, '-20', iStr, '/'];
    saveName = ['PCA_sw_MA_20', iMMStr, '-20', iStr, '.mat'];
    useSW = 1;
    useLW = 0;
    
    if i == 0
        searchStr = 'b30.042a.cam2.h0.200*.nc';
    else
        searchStr = ['b30.042a.cam2.h0.20', iMMStr, '*.nc b30.042a.cam2.h0.20', iStr(1), '*.nc'];
    end
    
    %run PCA code
    lou_PCA_moving_avg(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr)
    
end