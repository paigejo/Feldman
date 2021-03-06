%This script will generate subfolders in the given directory, where each
%folder is used to store .nc files for different variables in a given time
%step.  Then, it will split up the .nc files in the given directory into
%one .nc file for each time step, putting those files in the corresponding
%folders this script created.  Returns the time directory names this script
%creates for use in combine_nc_files, which will combine all files in those
%directories, and is meant to be used afterwards (although
%rotated_to_normal_coords may need to be called before combined_nc_files if
%any variables are in rotated latitude and lonitude coordinates).

%NOTE0: Call extrapolate_seasonal_nc_data before calling this function.
%Call combine_nc_files after calling this function.

%NOTE1: This script assumes that the nc files in the given directory are
%the exact variables that will be used for RadInput for the exact time
%frame. If any files have a different time frame from the rest, an error
%will be produced.

%NOTE2: This script uses a strsplit method that was downloaded from MATLAB
%file exchange and may be different from the one built in to the current
%version of MATLAB.

%NOTE3: fixed variables (not dependent on time) are copied to each
%subdirectory created.

function timeDirectories = split_nc_files(dirPath)
cd(dirPath);

%get all files in dirPath
[~, files] = system('ls');
files = strsplit(files, sprintf('\n'));

for fid = 1:length(files)
    file = files{fid};
    ['splitting file: ', file]
    
    if fid == 1
        %Determine the time range of the files, make folders:
        
        %%%Determine time range:
        
        %%%%%get time str
        [~, ncdump] = system(['/opt/local/bin/ncdump -v time ' file]);
        split1 = strsplit(ncdump, 'time = ');
        rawTime = split1{end};
        
        %%%%%only include time str before last space
        spaces = strfind(rawTime, ' ');
        commaTimes = rawTime(1:(spaces(end)-1));
        
        %%%%%get all time strings
        times = strsplit(commaTimes, ', ');
        
        %%%make folders:
        for t = 0:(length(times)-1)
            system(['mkdir ', num2str(times{t+1})]);
        end
    end
    
    %split nc file, put pieces in the correct timestep folders:
    
    %%%get file string root (everything before time part):
    breaks = strfind(file, '_');
    fileStrRoot = file(1:breaks(end));
    
    %%%now do the splitting and putting in the correct folder:
    for t = 0:(length(times)-1)
        tIStr = num2str(t);
        tStr = num2str(times{t+1});
        
        %if sftlf file, don't split, simply copy into the time directories
        if strncmp(fileStrRoot, 'sftlf', 5)
            system(['cp ', file, ' ', tStr, '/']);
            
        else
            system(['/opt/local/bin/ncks -d time,', tIStr, ',', tIStr,' ', file, ' ', tStr, '/', fileStrRoot, tStr, '.nc']);
        end
    end
end

timeDirectories = times;
end