%This script is meant to be used after split_nc_files (although
%rotated_to_normal_coords may need to be called before this as well if any
%variables are in rotated latitude/longitude coordinates).  It takes in a
%directory that contains the specified subdirectories, whcih should be in a
%cell array.  This script combines the nc files in each of those
%subdirectories, creating a combined file for each subdirectory
%representing all variable information for the timestep corresponding to
%the subdirectory.  After this function is used, another called
%format_nc_files should be called to ensure all combined nc files are in
%the correct format for RadInput.F90.

%NOTE1: This script uses a strsplit method that was downloaded from MATLAB
%file exchange and may be different from the one built in to the current
%version of MATLAB.

%NOTE2: the snd variable gets renamed to sndLImon if its in the LImon CMOR
%table and sndOImon if its in the OImon CMOR table

function combine_nc_files(dirPath, subDirectories)
cd(dirPath);

for dir = subDirectories
    dirStr = dir{1};
    
    cd(dirStr);
    
    ['combining files in ', dirStr]
    
    %get the files in this subdirectory
    [~, files] = system('ls');
    files = strsplit(files, sprintf('\n'));
    
    %get file root name (everything after variable piece)
    file = files{1};
    breaks = strfind(file, '_');
    fileRoot = file(breaks(1):end);
    
    %find highest a variable with all necessary dimensions (e.g. clw variable)
    isCLW = strncmp(files, 'clw_', 4);
    clwFile = files{isCLW};
    
    %initialize combined file for this subdirectory
    combinedFName = ['combined', fileRoot];
    system(['cp ', clwFile, ' ', combinedFName]);
    
    %check if sic data was unrotated (in which case the original sic file
    %won't be added to the combined file, only the unrotated will be)
    if sum(strncmp(files, 'sicUnrotated_', 13)) == 1
        unrotated = 1;
        
    else
        unrotated = 0;
    end
    
    %add information from all files to it:
    for file = files
        f = file{1};
        
        %skip sic file if sicUnrotated file was generated
        if strncmp(f, 'sic_', 4) && unrotated
            continue;
            
        end
        
        %combine file
        system(['/opt/local/bin/ncks -A ', f, ' ', combinedFName]);
    end
    
    %cd back to parent folder
    cd(dirPath);
end
end 