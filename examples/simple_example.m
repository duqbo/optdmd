%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple example --- how to call code
%
% Here we fit data generated from 3 
% spatial modes, each with time dynamics 
% which are exponential in time
%
% The examples show how to call the optdmd
% wrapper with various options
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate synthetic data

% set up modes in space

x0 = 0;
x1 = 1;
nx = 200;

% space

xspace = linspace(x0,x1,nx);

% modes

f1 = sin(xspace);
f2 = cos(xspace);
f3 = tanh(xspace);

% set up time dynamics

t0 = 0;
t1 = 1;
nt = 100;

ts = linspace(t0,t1,nt);

% eigenvalues

e_1 = 1;
e_2 = -2;
e_3 = 1i;

evals = [e_1;e_2;e_3];

% create clean dynamics

xclean = f1'*exp(e_1*ts) + f2'*exp(e_2*ts) + f3'*exp(e_3*ts);

% add noise (just a little ... this problem is 
% actually pretty challenging)

sigma = 1e-3;
xdata = xclean + sigma*randn(size(xclean));

%% compute modes in various ways

% target rank

r = 3;

% 1 --- fit to unprojected data

imode = 1;
[w,e1,b] = optdmd(xdata,ts,r,imode);

% reconstructed values
x1 = w*diag(b)*exp(e1*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e1,evals);
relerr_e = norm(e1(indices)-evals)/norm(evals);

fprintf('example 1 --- fitting unprojected data\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

% 2 -- fit to projected data (routine computes
% pod modes)

imode = 2;
[w,e2,b] = optdmd(xdata,ts,r,imode);

% reconstructed values
x1 = w*diag(b)*exp(e2*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e2,evals);
relerr_e = norm(e2(indices)-evals)/norm(evals);

fprintf('example 2 --- fitting data projected by optdmd on POD modes\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

%% 3 -- fit to projected data (basis computed 
% using randomized methods)

% let's set some optimization parameters

maxiter = 30; % maximum number of iterations
tol = sigma/100; % tolerance of fit
eps_stall = 1e-9; % tolerance for detecting a stalled optimization
opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',eps_stall);

% generate a pod basis using a randomized svd
% and pass this to routine (now the routine 
% won't have to do this calculation)

[u,~,~] = randsvd2(xdata,r,3);

imode = 2;
[w,e3,b] = optdmd(xdata,ts,r,imode,opts,[],u);

% reconstructed values
x1 = w*diag(b)*exp(e3*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e3,evals);
relerr_e = norm(e3(indices)-evals)/norm(evals);

fprintf('example 3 --- fitting data projected on POD modes from randsvd2\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

%% 4 -- use your own initial guess

% set a random initial guess (not generally a good idea)

e_init = randn(3,1);

imode = 2;
[w,e4,b] = optdmd(xdata,ts,r,imode,[],e_init);

% reconstructed values
x1 = w*diag(b)*exp(e4*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e4,evals);
relerr_e = norm(e4(indices)-evals)/norm(evals);

fprintf('example 4 --- using your own initial guess\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

%% 5 -- add linear constraints

% set a random initial guess (not generally a good idea)
iseed = 8675309;
rng(iseed);
e_init = randn(3,1);

% bounds 

% note that the first ia bounds apply to the real part of
% alpha and the second ia bounds apply to the imaginary part
% For unbounded, use the appropriate choice of +/- Inf

% the below has the effect of constraining the alphas to the
% left half plane

lbc = [-Inf*ones(size(e_init)); -Inf*ones(size(e_init))];
ubc = [zeros(size(e_init)); Inf*ones(size(e_init))];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

imode = 2;
[w,e5,b] = optdmd(xdata,ts,r,imode,[],e_init,[],copts);

% reconstructed values
x1 = w*diag(b)*exp(e5*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e5,evals);
relerr_e = norm(e5(indices)-evals)/norm(evals);

fprintf('example 5 --- alphas restricted to left-half plane\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

%% 6 -- add Tikhinov regularization

% set a random initial guess (not generally a good idea)
iseed = 8675309;
rng(iseed);
e_init = randn(3,1);

% tikhinov regularization parameter
gamma = 0.05;

imode = 2;
[w,e6,b] = optdmd(xdata,ts,r,imode,[],e_init,[],[],gamma);

% reconstructed values
x1 = w*diag(b)*exp(e6*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e6,evals);
relerr_e = norm(e6(indices)-evals)/norm(evals);

fprintf('example 6 --- alphas with Tikhinov regularization\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

%% 7 -- add proximal operator

% set a random initial guess (not generally a good idea)
iseed = 8675309;
rng(iseed);
e_init = randn(3,1) + 1i*randn(3,1);

% the below has the effect of constraining the alphas to the
% left half plane (like the linear constraints above)

proxfun = @(alpha) min(real(alpha),0) + 1i*imag(alpha);

imode = 2;
[w,e7,b] = optdmd(xdata,ts,r,imode,[],e_init,[],[],[],proxfun);

% reconstructed values
x1 = w*diag(b)*exp(e7*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e7,evals);
relerr_e = norm(e7(indices)-evals)/norm(evals);

fprintf('example 7 --- alphas restricted to left-half plane w/ prox\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)


%% 8 -- show constrained vs regularized vs optimal

figure()
set(groot, 'defaultLineMarkerSize',15)
scatter(real(evals),imag(evals),'bo')
hold on
scatter(real(e4),imag(e4),'rd')
scatter(real(e5),imag(e5),'k+')
scatter(real(e6),imag(e6),'mx')
legend({'true eigenvalues','opt dmd','opt dmd + constraint',...
    'opt dmd w/ regularizer'},'Location','NorthWest')
title('compare eigenvalues')

