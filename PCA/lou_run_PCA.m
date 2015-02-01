%this script performs all moving average PCA analysis on Lou (for all
%decades)

function lou_run_PCA(useSW, useLW)

swPath = '/lou/s2m/dfeldman/osse_sw/';
lwPath = '/lou/s2m/dfeldman/osse_lw/';

dbstop if error;

commandStr = '';
if useSW
    commandStr = 'sw';
end
if useLW
    commandStr = [commandStr, 'lw'];
end
    

for i = 0:9
    
    if i < 5 %only do calculations on new decades
        continue;
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
    savePath = ['/lou/s2j/jlpaige/all_', commandStr, '/b30.042a.cam2/20', iMMStr, '-20', iStr, '/'];
    saveName = ['PCA_', commandStr, '_MA_20', iMMStr, '-20', iStr, '.mat'];
    
    if i == 0
        searchStr = 'b30.042a.cam2.h0.200*.nc';
    else
        searchStr = ['b30.042a.cam2.h0.20', iMMStr, '*.nc b30.042a.cam2.h0.20', iStr(1), '*.nc'];
    end
    
    %run PCA code
    lou_PCA_moving_avg(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr)
    
end

end