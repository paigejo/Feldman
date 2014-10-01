%this function deletes the variable with name varName from the nc file, if
%it exists names fName (assumes that either fName is in the current
%directory, or fName is the path to a file

function delete_nc_variable(fName, varName)

%delete variable, if it exists
system(['/opt/local/bin/ncks -x -v ', varName, ' ', fName]);
end