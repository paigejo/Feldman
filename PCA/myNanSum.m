%myNanSum(X, dim)
%same as nansum, but dim can be a set of dimensions.

function M = myNanSum(X, dim)

M = X;
for i = 1:length(dim)
    d = dim(i);
    
    M = nansum(M, d);
    
end

end