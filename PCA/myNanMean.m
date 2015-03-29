%myNanMean(X, dim)
%same as nanmean, but dim can be a set of dimensions.  The returned average
%is then weighted based on the number of NaNs.

function M = myNanMean(X, dim)

tot = ~isnan(X);
M = X;
for i = 1:length(dim)
    d = dim(i);
    
    tot = sum(tot, d);
    M = nansum(M, d);
    
end
M = bsxfun(@rdivide, M, tot);

end





