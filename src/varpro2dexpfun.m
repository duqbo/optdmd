
function A = varpro2dexpfun(alpha,t,i)
%
% Derivatives of the matrix of exponentials
%
% Input 
% 
% alpha - vector of exponent values
% t - vector of times
% i - the desired derivative
%
% Output
%
% If Phi_i,j = exp(alpha_j t_i)
% then A = d/d(alpha_i) Phi in sparse
% format
%
m = length(t);
n = length(alpha);

if (i < 1 || i > n) 
    error('varpro2dexpfun: invalid index')
end

%A = zeros(m,n);
A = sparse(m,n);

ttemp = reshape(t,m,1);

A(:,i) = ttemp.*exp(alpha(i)*ttemp);

end
