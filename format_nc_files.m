%This function converts all files generated from the split_nc_files,
%combined_nc_files sequence to nc files with variables in the correct
%format (i.e. correct number of dimensions and units).  This function will
%generate files with the prefix 'final' for each combined file in the
%corresponding subdirectories.  Variables will be in form (time, lev, lat,
%lon), (time, lat, lon), or (lat, lon).  The input files should have only 1
%time coordinate.  If a required input file doesn't exist, a
%globally-uniform fill value is chosen (most often 0).

%correcting variable formatting:
%get_nc_variable -> ensure3D (or ensure2D) -> ensureCorrectDimensions -> unit conversion

%get_nc_variable: makes sure whatever dimensions exist are in correct order

%ensure3D: makes sure all dimensions that should exist do by copying data
%to different dimensions.  Makes sure variable is approximately the correct
%size, but not necessarily exact since interpolation may be necessary.

%ensureCorrectDimensions: interpolates data to make sure its exactly the
%right dimensions

%unit conversion: converts all units to correct units as required by
%RadInput

function format_nc_files(dirPath, subDirectories)
variableList = {'cfc11', 'cfc12', 'ch4', 'cl', 'cli', 'clw', 'n2o', ...
    'tro3', 'hus', 'hur', 'ta', 'sic', 'sftlf', 'ps', 'ts', 'tauu', ...
    'tauv', 'tas', 'sfcWind', 'snw'};

egLat = [-88.9277353522959, -87.5387052130272, -86.1414721015279, ...
    -84.7423855907142, -83.3425960440704, -81.9424662991732, ...
    -80.5421464346171, -79.1417096486217, -77.7411958655138, ...
    -76.3406287023715, -74.9400230196494, -73.5393886337675, ...
    -72.1387322891624, -70.7380587725176, -69.3373715749609, ...
    -67.9366733025785, -66.5359659401756, -65.1352510260352, ...
    -63.7345297708429, -62.3338031405324, -60.9330719152074, ...
    -59.5323367318266, -58.1315981156439, -56.7308565037137, ...
    -55.3301122627028, -53.9293657025561, -52.5286170870997, ...
    -51.1278666423533, -49.7271145631097, -48.3263610181882, ...
    -46.9256061546646, -45.5248501013023, -44.1240929713558, ...
    -42.723334864877, -41.3225758706231, -39.9218160676465, ...
    -38.5210555266244, -37.1202943109788, -35.719532477824, ...
    -34.3187700787707, -32.918007160614, -31.5172437659226, ...
    -30.1164799335463, -28.7157156990552, -27.3149510951204, ...
    -25.9141861518467, -24.5134208970629, -23.1126553565776, ...
    -21.7118895544042, -20.3111235129604, -18.9103572532454, ...
    -17.5095907949986, -16.1088241568413, -14.7080573564048, ...
    -13.3072904104462, -11.9065233349538, -10.5057561452436, ...
    -9.10498885604852, -7.70422148160049, -6.3034540357076, ...
    -4.90268653182654, -3.5019189831313, -2.10115140257898, ...
    -0.700383802973324, 0.700383802973324, 2.10115140257898, 3.5019189831313, ...
    4.90268653182654, 6.3034540357076, 7.70422148160049, 9.10498885604852, ...
    10.5057561452436, 11.9065233349538, 13.3072904104462, 14.7080573564048, ...
    16.1088241568413, 17.5095907949986, 18.9103572532454, 20.3111235129604, ...
    21.7118895544042, 23.1126553565776, 24.5134208970629, 25.9141861518467, ...
    27.3149510951204, 28.7157156990552, 30.1164799335463, 31.5172437659226, ...
    32.918007160614, 34.3187700787707, 35.719532477824, 37.1202943109788, ...
    38.5210555266244, 39.9218160676465, 41.3225758706231, 42.723334864877, ...
    44.1240929713558, 45.5248501013023, 46.9256061546646, 48.3263610181882, ...
    49.7271145631097, 51.1278666423533, 52.5286170870997, 53.9293657025561, ...
    55.3301122627028, 56.7308565037137, 58.1315981156439, 59.5323367318266, ...
    60.9330719152074, 62.3338031405324, 63.7345297708429, 65.1352510260352, ...
    66.5359659401756, 67.9366733025785, 69.3373715749609, 70.7380587725176, ...
    72.1387322891624, 73.5393886337675, 74.9400230196494, 76.3406287023715, ...
    77.7411958655138, 79.1417096486217, 80.5421464346171, 81.9424662991732, ...
    83.3425960440704, 84.7423855907142, 86.1414721015279, 87.5387052130272, ...
    88.9277353522959];

egLon = [0, 1.40625, 2.8125, 4.21875, 5.625, 7.03125, 8.4375, 9.84375, 11.25, ...
    12.65625, 14.0625, 15.46875, 16.875, 18.28125, 19.6875, 21.09375, 22.5, ...
    23.90625, 25.3125, 26.71875, 28.125, 29.53125, 30.9375, 32.34375, 33.75, ...
    35.15625, 36.5625, 37.96875, 39.375, 40.78125, 42.1875, 43.59375, 45, ...
    46.40625, 47.8125, 49.21875, 50.625, 52.03125, 53.4375, 54.84375, 56.25, ...
    57.65625, 59.0625, 60.46875, 61.875, 63.28125, 64.6875, 66.09375, 67.5, ...
    68.90625, 70.3125, 71.71875, 73.125, 74.53125, 75.9375, 77.34375, 78.75, ...
    80.15625, 81.5625, 82.96875, 84.375, 85.78125, 87.1875, 88.59375, 90, ...
    91.40625, 92.8125, 94.21875, 95.625, 97.03125, 98.4375, 99.84375, 101.25, ...
    102.65625, 104.0625, 105.46875, 106.875, 108.28125, 109.6875, 111.09375, ...
    112.5, 113.90625, 115.3125, 116.71875, 118.125, 119.53125, 120.9375, ...
    122.34375, 123.75, 125.15625, 126.5625, 127.96875, 129.375, 130.78125, ...
    132.1875, 133.59375, 135, 136.40625, 137.8125, 139.21875, 140.625, ...
    142.03125, 143.4375, 144.84375, 146.25, 147.65625, 149.0625, 150.46875, ...
    151.875, 153.28125, 154.6875, 156.09375, 157.5, 158.90625, 160.3125, ...
    161.71875, 163.125, 164.53125, 165.9375, 167.34375, 168.75, 170.15625, ...
    171.5625, 172.96875, 174.375, 175.78125, 177.1875, 178.59375, 180, ...
    181.40625, 182.8125, 184.21875, 185.625, 187.03125, 188.4375, 189.84375, ...
    191.25, 192.65625, 194.0625, 195.46875, 196.875, 198.28125, 199.6875, ...
    201.09375, 202.5, 203.90625, 205.3125, 206.71875, 208.125, 209.53125, ...
    210.9375, 212.34375, 213.75, 215.15625, 216.5625, 217.96875, 219.375, ...
    220.78125, 222.1875, 223.59375, 225, 226.40625, 227.8125, 229.21875, ...
    230.625, 232.03125, 233.4375, 234.84375, 236.25, 237.65625, 239.0625, ...
    240.46875, 241.875, 243.28125, 244.6875, 246.09375, 247.5, 248.90625, ...
    250.3125, 251.71875, 253.125, 254.53125, 255.9375, 257.34375, 258.75, ...
    260.15625, 261.5625, 262.96875, 264.375, 265.78125, 267.1875, 268.59375, ...
    270, 271.40625, 272.8125, 274.21875, 275.625, 277.03125, 278.4375, ...
    279.84375, 281.25, 282.65625, 284.0625, 285.46875, 286.875, 288.28125, ...
    289.6875, 291.09375, 292.5, 293.90625, 295.3125, 296.71875, 298.125, ...
    299.53125, 300.9375, 302.34375, 303.75, 305.15625, 306.5625, 307.96875, ...
    309.375, 310.78125, 312.1875, 313.59375, 315, 316.40625, 317.8125, ...
    319.21875, 320.625, 322.03125, 323.4375, 324.84375, 326.25, 327.65625, ...
    329.0625, 330.46875, 331.875, 333.28125, 334.6875, 336.09375, 337.5, ...
    338.90625, 340.3125, 341.71875, 343.125, 344.53125, 345.9375, 347.34375, ...
    348.75, 350.15625, 351.5625, 352.96875, 354.375, 355.78125, 357.1875, ...
    358.59375];

egLev = [3.54463800000001, 7.38881350000001, 13.967214, 23.944625, ...
    37.2302900000001, 53.1146050000002, 70.0591500000003, 85.4391150000003, ...
    100.514695, 118.250335, 139.115395, 163.66207, 192.539935, 226.513265, ...
    266.481155, 313.501265000001, 368.817980000002, 433.895225000001, ...
    510.455255000002, 600.524200000003, 696.796290000003, 787.702060000003, ...
    867.160760000001, 929.648875000002, 970.554830000001, 992.5561];

%constants for unit conversion:
goalLat = 128;
goalLon = 256;
goalLev = 26;

molMassAir = 28.97; %g/mol
molMassCFC11 = 163.1268; %g/mol
molMassCFC12 = 175.1375; %g/mol
molMassCH4 = 16.04; %g/mol
molMassN2O = 44.013; %g/mol
snowDensity = 500; %kg/m^3

%exampleFilePath is the file path to an empty (expect for dimensions) b30
%example file, which will be copied and the reformatted variables generated
%by this function will then overwrite the variables in the copied example
%file.
exampleFilePath = '~/Desktop/cmip5/b30-empty.nc';

    function variable = ensure3D(variable)
        %Make sure variable is 3-dimensional (4D including time) by copying
        %values along singleton dimensions.  Time will not be included as a
        %dimension since it is a trailing singleton dimension.
        
        if isempty(variable)
            error('ensure3D called on empty variable')
        end
        
        %get rid of singleton time dimension
        %variable = squeeze(variable);
        
        %now copy values into necessary dimensions
        dims = ndims(variable);
        if numel(variable) == 1
            variable = ones(goalLon, goalLat, goalLev)*variable;
            
        elseif numel(variable) == length(variable)
            error('variable input to ensure3D is a vector, which should not happen');
            
        elseif dims == 2
            %lev dimension must be missing
            
            %copy values along lev dimension
            variable = repmat(variable, [1, 1, goalLev]);
            
        elseif dims > 3
            error('ensure3D is being called on a variable with more than 3 dimensions');
        end
        
    end

    function variable = ensure2D(variable)
        %Make sure variable is 2D (3D including time) by copying values
        %along singleton dimensions
        
        if isempty(variable)
            error('ensure3D called on empty variable')
        end
        
        %now copy values into necessary dimensions
        dims = ndims(variable);
        if numel(variable) == 1
            variable = ones(goalLon, goalLat)*variable;
            
        elseif numel(variable) == length(variable)
            error('variable input to ensure2D is a vector, which should not happen');
            
        elseif dims > 2
            error('ensure2D is being called on a variable with more than 2 dimensions');
        end
        
    end


for dir = subDirectories
    dirName = dir{1};
    ['opening ' dirName]
    
    %open subdirectory
    cd(dirPath);
    cd(dirName);
    
    %get combined file name:
    ['finding combined file']
    [~, fileRaw] = system('find . -maxdepth 1 -name combined_\*');
    files = strsplit(fileRaw, '/');
    combinedFile = files{end};
    
    %copy combined file, calling it 'format' file
    ['creating output file']
    breaks = strfind(combinedFile, '_');
    finalFile = ['formatted', combinedFile(breaks(1):end)];
    system(['cp ', exampleFilePath, ' ', finalFile]);
    
    %get lon, lat, lev dimensions (and ilev?)
    ['reading in dimension data']
    lon = get_nc_variable(combinedFile, 'lon');
    lat = get_nc_variable(combinedFile, 'lat');
    lev = get_nc_variable(combinedFile, 'lev');
    
    %convert lev to be on 0-1000 scale instead of 0-1 scale
    if max(lev) < 1
        lev = lev*1000;
    end
    
    %get plev and variables for converting from plev
    plev = get_nc_variable(combinedFile, 'plev');
    ap = get_nc_variable(combinedFile, 'ap');
    b = get_nc_variable(combinedFile, 'b');
    ps = get_nc_variable(combinedFile, 'ps');
    
    for v = variableList
        varName = v{1};
        
        if strcmp(varName, 'cfc11')
            
            'formatting 3D variables'
            
            %if cfc11 doesn't exist, cfc11global should exist.  In that
            %case, use that variable
            if ~nc_variable_exists(combinedFile, varName)
                varName = 'cfc11global';
            end
            
            ['formatting ', varName]
            
            'reading data'
            %get data from combined file if it exists, else use fill value
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert units from molar fraction to kg/kg
            var = var * (molMassCFC11/molMassAir);
            
            'ensuring 3D'
            %make cfc11 3d (4d including time)
            var = ensure3D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            'writing to output file'
            %overwrite variable
            overwrite_nc_variable(finalFile, 'CFC11', var, 'CFC11', 4);
            
        elseif strcmp(varName, 'cfc12')
            %if cfc12 doesn't exist, cfc12global should exist.  In that
            %case, use that variable
            if ~nc_variable_exists(combinedFile, varName)
                varName = 'cfc12global';
            end
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert units from molar fraction to kg/kg
            var = var * (molMassCFC12/molMassAir);
            
            %make cfc12 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing to output file'
            overwrite_nc_variable(finalFile, 'CFC12', var, 'CFC12', 4);
            
        elseif strcmp(varName, 'ch4')
            %if ch4 doesn't exist, ch4global should exist.  In that
            %case, use that variable
            if ~nc_variable_exists(combinedFile, varName)
                varName = 'ch4global';
            end
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert units from molar fraction to kg/kg
            var = var * (molMassCH4/molMassAir);
            
            %make ch4 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'CH4', var, 'CH4', 4);
            
        elseif strcmp(varName, 'cl')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert from pct to frac if necessary
            if max(var(:)) > 2
                var = var/100;
            end
            
            %make cl 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'CLOUD', var, 'CLOUD', 4);
            
        elseif strcmp(varName, 'cli')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                cli = get_nc_variable(combinedFile, varName);
                
            else
                cli = 0;
                
            end
            
            %make sure in kg/kg not g/kg
            if max(cli(:)) > .005
                cli = cli/1000;
            end
            
            %make clwvi 3d (4d including time)
            'ensuring 3D'
            Icli = ensure3D(cli);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                Icli = ensureCorrectDimensions(Icli, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            Icli(Icli < 0) = 0;
            Icli(Icli > 1) = 1;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'CLDICE', Icli, 'CLDICE', 4);
            
        elseif strcmp(varName, 'clw')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                clw = get_nc_variable(combinedFile, varName);
                
            else
                clw = 0;
                
            end
            
            %make sure in kg/kg not g/kg
            if max(clw(:)) > .05
                clw = clw/1000;
            end
            
            %make clwvi 3d (4d including time)
            'ensuring 3D'
            Iclw = ensure3D(clw);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                Iclw = ensureCorrectDimensions(Iclw, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            Iclw(Iclw < 0) = 0;
            Iclw(Iclw > 1) = 1;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'CLDLIQ', Iclw, 'CLDLIQ', 4);
            
            %calculate variable we want (cloud ice fraction)
            'calculating cloud ice fraction'
            FICE = Icli./(Icli + Iclw);
            
            %if cli + clw close to 0, may take NaN values.  In that case,
            %set to 0
            FICE(isnan(FICE)) = 0;
            
            %ensure between 0 and 1
            FICE(FICE > 1) = 1;
            FICE(FICE < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'FICE', FICE, 'FICE', 4);
            
            %now compute ICLDWP:
            'formatting ICLDWP'
            
            %determine layer thickness as function of lon, lat, lev
            'Determining Layer Thickness'
            ap_bnds = ncread(combinedFile, 'ap_bnds');
            b_bnds = ncread(combinedFile, 'b_bnds');
            layerThickness = computeLayerThickness(ap_bnds, b_bnds, ps);
            
            %compute air densities
            'computing air densities'
            F = airDensity();
            pressure = hybridSigma2Pressure(ap, b, ps);
            pressureVec = reshape(pressure, numel(pressure), 1);
            densityVec = F(pressureVec);
            density = reshape(densityVec, size(cli, 1), size(cli, 2), size(cli, 3));
            
            %calculate ICLDWP
            ICLDWP = (clw + cli) .* density .* layerThickness;
            
            %convert from kg/m2 to g/m2 if necessary
            if max(var(:)) < 25
                var = var*1000;
            end
            
            %make clwvi 3d (4d including time)
            'ensuring 3D'
            ICLDWP = ensure3D(ICLDWP);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                ICLDWP = ensureCorrectDimensions(ICLDWP, lat, lon, lev);
            end
            
            %ensure bigger than 0
            ICLDWP(ICLDWP < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'ICLDWP', ICLDWP, 'ICLDWP', 4);
            
        elseif strcmp(varName, 'n2o')
            %if n2o doesn't exist, n2oglobal should exist.  In that
            %case, use that variable
            if ~nc_variable_exists(combinedFile, varName)
                varName = 'n2oglobal';
            end
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert units from molar fraction to kg/kg
            var = var * (molMassN2O/molMassAir);
            
            %make n2o 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, lev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writting data'
            overwrite_nc_variable(finalFile, 'N2O', var, 'N2O', 4);
            
        elseif strcmp(varName, 'tro3')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert from ppb to molar fraction
            var = var*10^(-9);
            
            %make tro3 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %convert from plev to lev
            'converting from plev to lev coordinates'
            var = plev2lev_alternate(var, plev);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, egLev);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'O3', var, 'O3', 4);
            
        elseif strcmp(varName, 'hus')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %if necessary, convert units from g/kg to kg/kg
            if(max(var(:)) > 1)
                var = var / 1000;
            end
            
            %make his 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %convert from plev to lev
            'converting from plev to lev coordinates'
            var = plev2lev_alternate(var, plev);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, egLev);
            end
            
            %ensure values between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'Q', var, 'Q', 4);
            
        elseif strcmp(varName, 'hur')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert to percentage if necessary
            if max(var(:)) < 2
                var = var * 100;
            end
            
            %make hur 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %convert from plev to lev
            'converting from plev to lev coordinates'
            var = plev2lev_alternate(var, plev);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, egLev);
            end
            
            %make sure var between 0 and 100
            var(var < 0) = 0;
            var(var > 100) = 100;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'RELHUM', var, 'RELHUM', 4);
            
        elseif strcmp(varName, 'ta')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 287; %fill value is average global surface temperature in Kelvin
                
            end
            
            %no unit conversion necessary, since already in K
            
            %make ta 3d (4d including time)
            'ensuring 3D'
            var = ensure3D(var);
            
            %convert from plev to lev
            'converting from plev to lev coordinates'
            var = plev2lev_alternate(var, plev);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, egLev);
            end
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'T', var, 'T', 4);
            
            
            
            
            
            %%%%%%%%%%%
            %Switching to 2D variables (3D including time):
            %%%%%%%%%%%
            'formatting 2D variables'
            
            
            
            
        elseif strcmp(varName, 'sic')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %make sure in fraction units not pct
            if max(var(:)) > 2
                var = var/100;
            end
            
            %make sic 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %overwrite variable
            'writing data'
            overwrite_nc_variable(finalFile, 'ICEFRAC', var, 'ICEFRAC', 3);
            
        elseif strcmp(varName, 'sftlf')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert from pct to frac if necessary
            if max(var(:)) > 2
                var = var/100;
            end
            
            %make sftlf 2d
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size, use 0 as
            %fill-value so values over ocean are set to 0
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, 0);
            end
            
            %ensure between 0 and 1
            var(var > 1) = 1;
            var(var < 0) = 0;
            
            %Write variable:
            'writting data'
            overwrite_nc_variable(finalFile, 'LANDFRAC', var, 'LANDFRAC', 2);
            
        elseif strcmp(varName, 'ps')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 101325; %this is the mean sea-level pressure in Pa
                
            end
            
            %Units don't need to be converted from Pascals, theoretically
            
            %make ps 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'PS', var, 'PS', 3);
            
        elseif strcmp(varName, 'ts')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 287; %fill value is average global surface temperature in Kelvin
                
            end
            
            %units don't need to be converted from K
            
            %make ts 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'TS', var, 'TS', 3);
            
        elseif strcmp(varName, 'tauu')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %units don't need to be converted from Pa
            
            %make tauu 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'TAUX', var, 'TAUX', 3);
            
        elseif strcmp(varName, 'tauv')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %units don't need to be converted from Pa
            
            %make tauv 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'TAUY', var, 'TAUY', 3);
            
        elseif strcmp(varName, 'tas')
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 287; %fill value is average global surface temperature in Kelvin
                
            end
            
            %units don't need to be converted from Pa
            
            %make tas 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, NaN);
            end
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'TREFHT', var, 'TREFHT', 3);
            
        elseif strcmp(varName, 'snw')
            %if snw doesn't exist, try to use lwsnl
            if ~nc_variable_exists(combinedFile, varName)
                varName = 'lwsnl';
            end
            
            ['formatting ', varName]
            
            %get data from combined file if it exists, else use fill value
            'reading data'
            if nc_variable_exists(combinedFile, varName)
                var = get_nc_variable(combinedFile, varName);
                
            else
                var = 0;
                
            end
            
            %convert from kg/m^2 to m assuming snow density fo 500kg/m^3.
            %This may be a problematic assumption, especially when using
            %lwsnl variable instead of snw
            var = var/snowDensity;
            
            %make snw 2d (3d including time)
            'ensuring 2D'
            var = ensure2D(var);
            
            %Use interpolation to make it exactly correct size, fill in
            %missing data over oceans with zeros
            if nc_variable_exists(combinedFile, varName)
                'ensuring dimensions correct size'
                var = ensureCorrectDimensions(var, lat, lon, 0);
            end
            
            %ensure bigger than 0
            var(var < 0) = 0;
            
            %Write variable:
            'writing data'
            overwrite_nc_variable(finalFile, 'SNOWHLND', var, 'SNOWHLND', 3);
            
        end
    end
end
end
