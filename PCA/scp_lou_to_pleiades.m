%This function moves files from the given directory in Lou to the given
%directory in Pleiades.  This function should be called from Lou.

function scp_lou_to_pleiades(louDir, pleiadesDir)

%get files
cd(louDir);
[~, files] = system('ls');
files = strsplit(files, sprintf('\n'));

%send files one by one
for f = 1:length(files)
    file = files{f};
    
    %print progress
    disp(['sending file ', file, ' (', num2str(f), '/', num2str(length(files)), ')']);
    
    %get file data from tape and send to pleiades via scp
    system(['dmget ', file]);
    system(['scp ', file, ' pfe:', pleiadesDir]);
end