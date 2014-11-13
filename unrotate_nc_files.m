%This function unrotates all 'sic' variables in the specified
%subdirectories of dirPath.  The 'sic' files are given new prefixes of
%'sicUnrotated'.  This function uses rotated_to_normal_coords.m to convert
%from rotated lat/lon coordinates to normal lat/lon coordinates.

%NOTE: If the specified subdirectories have been manually altered after
%having been created by split_nc_files aside from by using the functions in
%this folder, the function may have problems.  Specifically, if there is no
%sic file, or if there are multiple sic files, this function may produce
%unwanted results.

%NOTE: This function assumes the subdirectories contain files with the ts
%variable that are not in rotated lat/lon coordinates.

function unrotate_nc_files(dirPath, subDirectories)

cd(dirPath);

for dir = subDirectories
    dirStr = dir{1};
    
    cd(dirStr);
    
    ['unrotating sic file in ', dirStr, ' directory']
    
    %get the sic file in this subdirectory
    [~, sicFile] = system('ls sic_*');
    sicFile = strsplit(sicFile, sprintf('\n'));
    sicFile = sicFile{1};
    
    %get a file with correct lat/lon coordinates file in this subdirectory
    %(such as ts)
    [~, tsFile] = system('ls ts_*');
    tsFile = strsplit(tsFile, sprintf('\n'));
    tsFile = tsFile{1};
    
    %get sample lat/lon data
    egLat = get_nc_variable(tsFile, 'lat');
    egLon = get_nc_variable(tsFile, 'lon');
    
    %unrotate sic variable
    rotated_to_normal_coords(sicFile, 'sic', egLat, egLon);
    
    %cd back to parent folder
    cd(dirPath);
end
end