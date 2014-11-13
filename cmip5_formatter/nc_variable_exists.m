%This function checks to see if the given variable named varName exists in
%the .nc file named fName.  Assumes fName is in current directory or is a
%filepath to the nc file.

function exists = nc_variable_exists(fName, varName)

%get raw ncdump string
[~, ncdump] = system(['/opt/local/bin/ncdump -h ', fName]);

%find where variable dimensions are listed (and convert to lower case)
found = strfind(ncdump, [varName, '(']);

%if not empty, variable exists, otherwise doesn't
exists = ~isempty(found);
end