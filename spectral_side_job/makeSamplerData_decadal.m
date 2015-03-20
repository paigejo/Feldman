
%save results here:
savePath = '/global/u1/j/jpaige/git/Feldman/spectral_side_job/results';
cd(savePath);

for i = 0:9
    
    decadeStr = ['20', num2str(i), '0-20', num2str(i), '9'];
    
    %generate decadal data for different values of allSky, lwHiRes
    sampleData_decadal(decadeStr, true, true);
    sampleData_decadal(decadeStr, true, false);
    sampleData_decadal(decadeStr, false, true);
    sampleData_decadal(decadeStr, false, false);
    
end