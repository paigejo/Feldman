%This function deletes the current variable named varName in the file given
%in fName if it exists, and then recreates it with the data in newData.
%This is necessary in order to deal with difference in variable dimensions
%between the old and new variable data.  If you want to create a variable
%or overwrite it if it already exists, call this function, not
%create_nc_variable.

%NOTE1: the new data for the given variable will be double-precision.

%NOTE2: assumes that the dimensions of the variable will either be (time),
%(lat, lon), (time, lat, lon), or (time, lev, lat, lon) if the number of
%dimensions of newData are respectively: vector, 2, 3, or 4.  If newData is
%single element, no dimensions are used.

function overwrite_nc_variable(fName, varName, newData, newName)

%delete old variable
[~] = evalc('delete_nc_variable(fName, varName)');

%make new variable with new data (use evalc to not output way too much
%data.  Note: this may cause errors to go unnoticed until later
[~] = evalc('create_nc_variable(fName, newName, newData)');
end