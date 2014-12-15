%[newX, newY] = compute_mean_series(x, y) converts a series with multiple
%points per x value to a series with one point per x value, each point at a
%given x value having a y value that is the mean of the y values with that
%x value in the original series.

function [newX, newY] = compute_mean_series(x, y)

[uniqueX, IX, ~] = unique(x);
newX = 1:length(IX);
newY = newX;

for i = 1:length(IX)
    ind = IX(i);
    
    lastInd = 0;
    if i > 1
        lastInd = IX(i - 1);
    end
    
    newX(i) = uniqueX(i);
    newY(i) = mean(y( (lastInd+1):ind ));
end
end





