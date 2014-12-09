%This function moves files from the current directory in Lou to the given
%directory in Pleiades.  This function should be called from Lou. Remote
%str should be the string used in scp to identify a username and send
%directory to send the data to.  pw is the password string.

function scp_lou_to_remote(remoteStr, pw)

%get files
[~, files] = system('ls');
files = strsplit(files, sprintf('\n'));

%update password variable
system(['export RSYNC_PASSWORD=', pw]);

%send files one by one
for f = 1:length(files)
    file = files{f};
    
    %print progress
    disp(['sending file ', file, ' (', num2str(f), '/', num2str(length(files)), ')']);
    
    %get file data from tape and send to remote via scp
    system(['dmget ', file]);
    system(['scp ', file, ' ', remoteStr]);
end