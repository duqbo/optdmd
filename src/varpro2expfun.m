
function A = varpro2expfun(alpha,t)
%
% matrix of exponentials
%
% Input 
% 
% alpha - vector of exponent values
% t - vector of times
%
% Output
%
% A(i,j) = exp(alpha_j t_i)
%

m = length(t);
n = length(alpha);

A = zeros(m,n);

ttemp = reshape(t,m,1);
atemp = reshape(alpha,n,1);

temp = ttemp*transpose(atemp);

A = exp(temp);

end
