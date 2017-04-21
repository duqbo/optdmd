
function indices = match_vectors(v1,v2)
%MATCH_VECTORS Wrapper for MUNKRES
%
% Sets up a cost function so that the indices
% returned by munkres correspond to the permutation
% which minimizes the 1-norm of the difference 
% between v1(indices) and v2( (indices ~= 0) )
%
% Example:
%
%   >> indices = match_vectors(v1,v2)
%
% See also MUNKRES
v1 = v1(:);
v2 = v2(:);
costmat = abs(kron(v1.',ones(length(v2),1)) - kron(v2,ones(1,length(v1))));
[indices,cost] = munkres(costmat);

