%This function converts gas fraction variables stored in the dirPath
%directory from seasonal format (only 12 values are stored, one for each
%month) to a time series extending as far as the other files in the
%directory.  Call this function before calling split_nc_files.

function extrapolate_seasonal_nc_data(dirPath)
cd(dirPath);

seasonalVariables = {'cfc11', 'cfc11global', 'cfc12', 'cfc12global', ...
    'ch4', 'n2o', 'n2oglobal'};

%get all files in dirPath
[~, files] = system('ls');
files = strsplit(files, sprintf('\n'));

for fid = 1:length(files)
    file = files{fid};
    
    %get the variable of the file
    breaks = strfind(file, '_');
    varName = file(1:(breaks(1)-1));
    
    %check if it is a seasonal variable
    I = strcmp(seasonalVariables, varName);
    if sum(I) == 1
        
        %get variable and time
        variable = get_nc_variable(file, varName);
        time = get_nc_variable(file, 'time');
        
        %correct for if it is in seasonal or "one value" format
        if size(variable, 1) == 12
            
            %copy values throughout time
            variable = repmat(variable, [length(time)/12, ones(ndims(variable) - 1, 1)]);
            
            %overwrite variable
            overwrite_nc_variable(file, varName, variable, varName, 4);
            
        elseif ndims(variable) == 2 && (length(variable) ~= numel(variable))
            
            %add time dimension, since it's clearly not there
            variable = dimshift(variable, -1);
            
            %copy values throughout time
            variable = repmat(variable, [length(time), ones(ndims(variable) - 1, 1)]);
            
            %overwrite variable
            overwrite_nc_variable(file, varName, variable, varName, 4);
            
        elseif ndims(variable) == 3 && (size(variable, 1) == 1)
            
            %copy values across time (first dimension of variable)
            variable = repmat(variable, [length(time), ones(ndims(variable) - 1, 1)]);
            
            %overwrite variable
            overwrite_nc_variable(file, varName, variable, varName, 4);
            
        end
    end
end
end