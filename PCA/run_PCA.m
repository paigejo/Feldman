%this script performs all moving average PCA analysis on NERSC systems (for
%all decades)

function run_PCA(useSW, useLW, lwHiRes, normalize)

swPath = '/global/scratch2/sd/jpaige/PCA/sw_files/';
lwPath = '/global/scratch2/sd/jpaige/PCA/lw_files/';
lwHiRes = logical(lwHiRes);

dbstop if error;

commandStr = '';
if useSW
    commandStr = 'sw';
end
if useLW
    commandStr = [commandStr, 'lw'];  
end

if lwHiRes
    resStr = '_hres';
else
    resStr = '_lres';
end

if normalize
    normStr = 'cellnorm';
else
    normStr = 'notnorm';
end

for i = 0:9
    
    %already done with first 7 decades
    if i < 8
        continue;
    end
    
    %decade data specification strings
    decadeStr = [num2str(i)];
    
    disp(['Running PCA code on 20', decadeStr, '0-20', decadeStr, '9 decade']);
    
    %determine inputs for PCA function based on decade
    savePath = ['/global/scratch2/sd/jpaige/PCA/all_', commandStr, '/b30.042a.cam2/20', decadeStr, '0-20', decadeStr, '9/'];
    saveName = ['PCA_', commandStr, '_', normStr, '_20', decadeStr, '0-20', decadeStr, '9', resStr, '.mat'];
    searchStr = ['b30.042a.cam2.h0.20', decadeStr, '*.nc'];
    
    %run PCA code
    PCA_final(useSW, useLW, saveName, savePath, swPath, lwPath, searchStr, lwHiRes, normalize)
    
end

end