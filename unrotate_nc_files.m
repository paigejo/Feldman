%This function unrotates all 'sic' variables in the specified
%subdirectories of dirPath.  The 'sic' files are given new prefixes of
%'sicUnrotated'.  This function uses rotated_to_normal_coords.m to convert
%from rotated lat/lon coordinates to normal lat/lon coordinates.

%NOTE: If the specified subdirectories have been manually altered after
%having been created by split_nc_files aside from by using the functions in
%this folder, the function may have problems.  Specifically, if there is no
%sic file, or if there are multiple sic files, this function may produce
%unwanted results.

function unrotate_nc_files(dirPath, subDirectories)

cd(dirPath);

for dir = subDirectories
    dirStr = dir{1};
    
    cd(dirStr);
    
    %get the sic file in this subdirectory
    [~, files] = system('ls sic*');
    file = strsplit(files, sprintf('\n'));
    file = files;
    
    %unrotate sic variable
    rotated_to_normal_coords(file, 'sic');
    
    %cd back to parent folder
    cd(dirPath);
end
end