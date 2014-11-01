%returns the names of the dimensions of the given variable in the given
%file.  Note that the order of the dimension names is in NCO ordering,
%which is the reverse of the MATLAB ordering.

function dimNames = get_nc_variable_dimensions(fileName, varName)

%get raw ncdump string
[~, ncdump] = system(['/opt/local/bin/ncdump -h ', fileName]);

%find where variable dimensions are listed (and convert to lower case)
startI = strfind(ncdump, [varName, '(']) + length(varName) + 1;
endIs = strfind(ncdump(startI:end), ')');
endI = endIs(1) - 2 + startI;
dimString = lower(ncdump(startI:endI));
dimNames = strsplit(dimString, ' ');

end