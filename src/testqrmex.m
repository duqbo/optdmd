
%% test lapack qr routines 

m = 300;
n = 200;
nrhs = 10;

A = randn(m,n) + 1i*randn(m,n);
xtrue = randn(n,nrhs) + 1i*randn(n,nrhs);

if (1 == 1)
    A = real(A);
    xtrue = real(xtrue);
end

B = A*xtrue;

tic
[AOUT,jpvt,tau] = xgeqp3_m(A);
C = xormqr_m('L','T',AOUT,tau,B);
R = triu(AOUT);
y = R\C;
x = zeros(size(y));
x(jpvt(1:length(y)),:) = y;
toc

norm(x-xtrue)
norm(A*x-B)/norm(B)
