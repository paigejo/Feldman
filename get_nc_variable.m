%Given the nc file name (assumed to be in the current directory, or the
%file path) and the variable name in the nc file, returns the variable data
%from the nc file (converted from string to double form).  Also, corrects
%for format: makes sure the order of dimensions is correct (time, lev, lat,
%lon) or (time, lat, lon) in nco ordering, or the reverse order for Matlab.

%NOTE: for the system call to ncdump to work, you might need to add the 
%line setenv('DYLD_LIBRARY_PATH',''); to your statup.m file in MATLAB's
%startup folder.

%NOTE2: I think MATLAB orders dimensions in reverse order from ncdump

%NOTE3: matlab deletes trailing singleton dimensions, so if the last
%dimension in the matlab ordering (reverse ordering) of the variable is a
%singleton one, it will be ignored.

function [varData, nDims] = get_nc_variable(fileName, varName)

varData = ncread(fileName, varName);

if strcmp(varName, 'time') || strcmp(varName, 'lev') || ...
    strcmp(varName, 'ilev') || strcmp(varName, 'lat') || ...
    strcmp(varName, 'lon') || strcmp(varName, 'plev')
    
    %then rearranging data dimensions is unnecessary
    return;
end

% ensure data is in correct format:

%get raw ncdump string
[~, ncdump] = system(['/opt/local/bin/ncdump -h ', fileName]);

%find where variable dimensions are listed (and convert to lower case)
startI = strfind(ncdump, [varName, '(']) + length(varName) + 1;
endIs = strfind(ncdump(startI:end), ')');
endI = endIs(1) - 2 + startI;
dimString = lower(ncdump(startI:endI));
dims = strsplit(dimString, ' ');
nDims = numel(dims);

%determine dimension ordering
timeLoc = strfind(dimString, 'time');
levLoc = strfind(dimString, 'lev'); %Include ilev????
latLoc = strfind(dimString, 'lat');
lonLoc = strfind(dimString, 'lon');

%permute dimensions into correct order:
if nDims == 2 && (numel(varData) ~= length(varData))
    %then only lat and lon will be included
    
    if lonLoc < latLoc
        permute(varData, [2 1]);
    end
    
elseif nDims == 3
    %Then only time, lat, and lon will be included?????????
    
    %compute inverse permutation taking dim ordering back to normal
    [~, forwardPermute] = sort([-lonLoc, -latLoc, -timeLoc]);
    [~, inversePermute] = sort(forwardPermute);
    
    %permute dimensions to necessary order
    permute(varData, inversePermute);
    
elseif nDims == 4
    %Then we have time, lev (OR plev), lat, and lon
    
    %compute inverse permutation taking dim ordering back to normal
    [~, forwardPermute] = sort([-lonLoc, -latLoc, -levLoc, -timeLoc]);
    [~, inversePermute] = sort(forwardPermute);
    
    %permute dimensions to necessary order
    permute(varData, inversePermute);
end

end