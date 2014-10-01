%This function creates the specified variable in the specified nc file.  If
%it already exists in the file, this function returns an error.  In that
%case, call overwrite_nc_variable to overwrite the original variable, which
%does not return an error if the variable didn't exist in the first place.

%NOTE1: the new data for the given variable will be double-precision.

%NOTE2: assumes that the dimensions of the variable will either be (time),
%(lat, lon), (time, lat, lon), or (time, lev, lat, lon) if the number of
%dimensions of newData are respectively: vector, 2, 3, or 4.  If newData is
%single element, no dimensions are used.

%NOTE3: if variable already exists in file, nccreate returns error, and
%this function also returns an error.

function create_nc_variable(fName, varName, varData)

%determine number of dimensions that exists and thereby infer which
%dimensions are which
dims = ndims(varData);
if numel(varData) == 1
    nccreate(fName, varName, 'Datatype', 'single');
    ncwrite(fName, varName, varData);
    return;
    
elseif sum(strcmp(varName, {'lat', 'lon', 'lev', 'time'})) > 0
    dimCell = {varName, length(varData)};
    
elseif numel(varData) == length(varData)
    dimCell = {'time', Inf};
    
elseif dims == 2
    dimCell = {'lat', size(varData, 1), 'lon', size(varData, 2)};
    
elseif dims == 3
    dimCell = {'lat', size(varData, 2), 'lon', size(varData, 3), 'time', Inf};
    
    %make 'time' last dimension (because MATLAB needs it like that for some
    %reason)
    varData = permute(varData, [2, 3, 1]);
    
elseif dims == 4
    dimCell = {'lev', size(varData, 2), 'lat', size(varData, 3), ...
        'lon', size(varData, 4), 'time', Inf};
    
    %make 'time' last dimension (because MATLAB needs it like that for some
    %reason)
    varData = permute(varData, [2, 3, 4, 1]);
end
nccreate(fName, varName, 'Dimensions', dimCell, 'Datatype', 'single');
ncwrite(fName, varName, varData);

%if we permuted the variable data, we need to permute it back with ncpdq
if dims == 3
    system(['/opt/local/bin/ncpdq -O -a time,lat,lon -v ', varName, ...
        ' ', fName, ' ', fName]);
    
elseif dims == 4
    system(['/opt/local/bin/ncpdq -O -a time,lev,lat,lon -v ', varName, ...
        ' ', fName, ' ', fName]);
    
end

end