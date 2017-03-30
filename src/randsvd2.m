function [U,S,V] = randsvd2(A,p,q)
%
% An implementation of the random SVD
% algorithms 4.4.1 and 4.4 of
% "Finding structure with randomness"
%
% Input:
%
% A - matrix
% p - target rank p
% q - number of subspace iterations
%     to perform (the best number is
%     determined by the decay of the 
%     singular values of A. We recommend
%     taking q >= 2)
%
% Copyright Travis Askham 2017
%
% MIT License
%
[m,n] = size(A);
nr = min(max(2*p,p+5),n);
r = randn(n,nr);
Y = A*r;
[QY,~,~] = qr(Y,0);

if (q > 0)
    Ap = A';
end

% perform subspace iteration

for j = 1:q
    Y = Ap*QY;
    [QY,~,~] = qr(Y,0);
    Y = A*QY;
    [QY,~,~] = qr(Y,0);
end

AY = QY'*A;
[U,S,V] = svd(AY,0);
U = QY*U;
U = U(:,1:p);
S = S(1:p,1:p);
V = V(:,1:p);

end

