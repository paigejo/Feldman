%this function deletes the variable with name varName from the nc file, if
%it exists names fName (assumes that either fName is in the current
%directory, or fName is the path to a file.
%
%NOTE: This can delete multiple files if the varName given is
%'var1,var2,...,varn' with no spaces.

function delete_nc_variable(fName, varName)

%delete variable, if it exists
system(['/opt/local/bin/ncks -C -O -x -v ', varName, ' ', fName, ' ', fName]);
end