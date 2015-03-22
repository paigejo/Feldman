
%save results here:
savePath = '/global/u1/j/jpaige/git/Feldman/spectral_side_job/results';
cd(savePath);

for i = 0:9
    
    decadeStr = ['20', num2str(i)];
    
    disp(['Sampling for ' decadeStr, '0-', decadeStr, '9 decade'])
    
    %generate decadal data for different values of allSky, lwHiRes
    disp('all-sky, low-res')
    samplerData_decadal(decadeStr, true, false);
    disp('all-sky, hi-res')
    samplerData_decadal(decadeStr, true, true);
    disp('clear-sky, low-res')
    samplerData_decadal(decadeStr, false, false);
    disp('clear-sky, hi-res')
    samplerData_decadal(decadeStr, false, true);
    
    
end