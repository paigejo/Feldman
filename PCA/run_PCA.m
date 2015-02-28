%this script performs all moving average PCA analysis on NERSC systems (for
%all decades)

function run_PCA(useSW, useLW, lwHiRes)

swPath = '/global/scratch2/sd/jpaige/PCA/sw_files/';
lwPath = '/global/scratch2/sd/jpaige/PCA/lw_files/';

dbstop if error;

commandStr = '';
if useSW
    commandStr = 'sw';
end
if useLW
    commandStr = [commandStr, 'lw'];
end
    

for i = 0:9
    
    %decade data specification strings
    iStr = [num2str(i), '9'];
    iMMStr = [num2str(i-1), '9'];
    
    %no data from 1999, so must modify strings for first decade
    if i == 0
        iMMStr = '00';
    end
    
    disp(['Running PCA code on 20', iMMStr, '-20', iStr, ' decade']);
    
    %determine inputs for PCA function based on decade
    savePath = ['/global/scratch2/sd/jpaige/PCA/all_', commandStr, '/b30.042a.cam2/20', iMMStr, '-20', iStr, '/'];
    saveName = ['PCA_', commandStr, '_final_20', iMMStr, '-20', iStr, '.mat'];
    
    if i == 0
        searchStr = 'b30.042a.cam2.h0.200*.nc';
    else
        searchStr = ['b30.042a.cam2.h0.20', iMMStr, '*.nc b30.042a.cam2.h0.20', iStr(1), '*.nc'];
    end
    
    %run PCA code
    PCA_final(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr, lwHiRes)
    
end

end