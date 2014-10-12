%returns an interpolator based on US standard atmosphere air properties
%table found at
%http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html.

function F = airDensity()

%get pressure (convert from 10^4 Pa to Pa)
pressure = [11.39 10.13 8.988 7.950 7.012 6.166 5.405 4.722 4.111 3.565 ...
    3.080 2.650 1.211 .5529 .2549 .1197 .0287 .007978 .002196 .00052 .00011];
pressure = pressure' * 10^4;

%get air density (convert from 10^-1 kg/m^3 to kg/m^3)
density = [13.47 12.25 11.12 10.07 9.093 8.194 7.364 6.601 5.900 5.258 ...
    4.671 4.135 1.948 .8891 .4008 .1841 .03996 .01027 .003097 .0008283 .0001846];
density = density'/10;

%construct interpolator
F = scatteredInterpolator(pressure, density);
end