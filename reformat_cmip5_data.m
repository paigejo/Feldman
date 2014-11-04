%This script takes the user through converting from many CMIP5 output
%formats (so far, guaranteed to work for CanESM2) to the format required by
%the satellite spectroscopic simulator.  In order for the nco operators to
%work within the MATLAB environment, you might need to add the 
%line setenv('DYLD_LIBRARY_PATH',''); to your statup.m file in MATLAB's
%startup folder.  If this file does not already exist, you should create
%it.

%STEP 1: put all the variables for a given model and time period in a
%single directory.  Adjust the directory below to be your own directory.
%It is recommended that you copy the desired files into a test
%directory since they may be modified.

dataDir = '/Users/johnpaige/Desktop/cmip5/CanESM2/simvars/testDir';

%STEP 2: Make sure the variables that require time dimension have values
%for every time step (and not just a single value or values for each
%season).  The extrapolate_seasonal_nc_data.m function is meant to be used
%for this.

%extrapolate_seasonal_nc_data(dataDir);

%STEP 3: Make sure lon coordinates in sic are shifted to the correct range
%(i.e. [0, 360]).  The sic variable typically uses rotated lat/lon
%coordinates and sometimes uses strange lon ranges.

%{
cd(dataDir);
sicFile = 'sic_OImon_CanESM2_rcp85_r1i1p1_200601-210012.nc';
sicLon = get_nc_variable(sicFile, 'lon');
sicLon = sicLon + 280; %substitute shift for 280 here.
overwrite_nc_variable(sicFile, 'lon', sicLon, 'lon', {'i', size(sicLon, 1), 'j', sizeSicLon(, 2)});
%}

%STEP 4: Before beginning this step, make sure the dataDir is clear of
%anything not specifically used or output so far in this script.  The
%dataDir directory is about to get very cluttered, with one directory for
%each time step.  Run the split_nc_files.m function to create directories
%for each timestep such that variables in a directory for a timestep only
%contain data for that timestep.  Save the timeDirectories variable in case
%you have to quit and come back to this process later.

timeDirectories = split_nc_files(dataDir);
save('timeDirectories.mat', timeDirectories);

%STEP 5: Rotate sic lat/lon coordinates, if necessary.  Before running this
%operation check to make sure the lon/lat coordinates are actually rotated.
% The unrotated data will be put in a new file with 'sicUnrotated' as the
% file prefix.  Use ncview to make sure the created file is correct, then
% you may remove the original sic file if you like.

unrotate_nc_files(dataDir, timeDirectories)

%STEP 6: Combine the files in each of the newly created subdirectories into
%a single combined nc file with all of the data from each of the files for
%any given timestep.  The output file has the prefix 'combined' instead of
%a normal variable name as a prefix.

combine_nc_files(dataDir, timeDirectories);

%STEP 7: Complete the formatting process by calling the format_nc_files.m
%function.  The second line below can be uncommented to ignore the cfc11,
%cfc12, ch4, n2o, tro3, tauu, tauv, and sfcWind variables to improve
%performance, since these variables do not have very significant impact on
%spectroscopic observations. The format_nc_files function performs the
%final formatting steps on all the specified subdirectories, such as unit
%conversion, interpolation, dimension permutation, calculating variables
%based on other variables, and more.  The output file has the prefix
%'formatted' instead of a normal variable name as a prefix.  Note: this
%step can take a long time.  I recommend you try it using only one of the
%elements of timeDirectories first and seeing how long it takes.

variables = ones(1, 20);
%variables([1:3, 7:8, 16:17, 19]) = 0;
format_nc_files(dataDir, timeDirectories, variables);

%STEP 8: Use ncview to check to make sure everything is correct!






