%lou_general_command(searchStr, leftOfFile, rightOfFile) can be run on
%either Lou or Plieades.  It runs a command given by
%[leftOfFile, fileName, rightOfFile], where fileName is any file matching
%the search string searchStr.  searchStr is input to the ls system command
%and can include wildcards.  This function must be called when on Lou or
%Plieades when operating on many files because the dmget function must be
%used to ensure each file is retreived from the file system before
%operations can be run on it.
%
%e.g. lou_general_command('*b30.036a.cam2*.nc', 'cp ', ' ~/osse_sw/b30.036a.cam2/')
%The above line of code would copy all files in the current directory
%matching the '*b30.036a.cam2*.nc' search string and copy them into the
%~/osse_sw/b30.036a.cam2/ directory.

function lou_general_command(searchStr, leftOfFile, rightOfFile)

%get files
[~, files] = system(['ls ', searchStr]);
files = strsplit(files, sprintf('\n'));

%operate command on files one by one
for f = 1:length(files)
    file = files{f};
    
    %print progress
    disp(['operating on file ', num2str(f), ' out of ', num2str(length(files))]);
    
    %get file data from tape and run command on it
    system(['dmget ', file]);
    system([leftOfFile, file, rightOfFile]);
end

end