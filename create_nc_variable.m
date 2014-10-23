%This function creates the specified variable in the specified nc file.  If
%it already exists in the file, this function returns an error.  In that
%case, call overwrite_nc_variable to overwrite the original variable, which
%does not return an error if the variable didn't exist in the first place.
%Takes in the number of dimensions of the variable in case a singleton
%time dimension has been removed.

%NOTE1: assumes that the dimensions of the variable will either be (time),
%(lat, lon), (time, lat, lon), or (time, lev, lat, lon) if the number of
%dimensions of newData are respectively: vector, 2, 3, or 4.  If newData is
%single element, no dimensions are used. The number of dimensions is
%specified by dims, which can also be in the for {'dimName1', dimSize1,
%'dimName2', dimSize2, ...} for custom dimensions.

%NOTE2: if variable already exists in file, nccreate returns error, and
%this function also returns an error.

function create_nc_variable(fName, varName, varData, dims)

%determine number of dimensions that exists and thereby infer which
%dimensions are which
if numel(varData) == 1
    nccreate(fName, varName, 'Datatype', 'single', 'Format', 'classic');
    ncwrite(fName, varName, varData);
    return;
    
elseif iscell(dims)
    dimCell = dims;
    
elseif sum(strcmp(varName, {'lat', 'lon', 'lev', 'time'})) > 0
    dimCell = {varName, length(varData)};
    
elseif numel(varData) == length(varData)
    dimCell = {'time', Inf};
    
elseif dims == 2
    dimCell = {'lon', size(varData, 1), 'lat', size(varData, 2)};
    
elseif dims == 3
    dimCell = {'lon', size(varData, 1), 'lat', size(varData, 2), 'time', Inf};
    
elseif dims == 4
    dimCell = {'lon', size(varData, 1), 'lat', size(varData, 2), ...
        'lev', size(varData, 3), 'time', Inf};
    
end

nccreate(fName, varName, 'Dimensions', dimCell, 'Datatype', 'single', ...
    'Format', 'classic');
ncwrite(fName, varName, varData);

end