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

e1 = 1;
e2 = -2;
e3 = 1i;

evals = [e1;e2;e3];

% create clean dynamics

xclean = f1'*exp(e1*ts) + f2'*exp(e2*ts) + f3'*exp(e3*ts);

% add noise (just a little ... this problem is 
% actually pretty challenging)

sigma = 1e-3;
xdata = xclean + sigma*randn(size(xclean));

%% compute modes in various ways

% target rank

r = 3;

% 1 --- fit to unprojected data

imode = 1;
[w,e,b] = optdmd(xdata,ts,r,imode);

% reconstructed values
x1 = w*diag(b)*exp(e*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e,evals);
relerr_e = norm(e(indices)-evals)/norm(evals);

fprintf('example 1 --- fitting unprojected data\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

% 2 -- fit to projected data (routine computes
% pod modes)

imode = 2;
[w,e,b] = optdmd(xdata,ts,r,imode);

% reconstructed values
x1 = w*diag(b)*exp(e*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e,evals);
relerr_e = norm(e(indices)-evals)/norm(evals);

fprintf('example 2 --- fitting data projected by optdmd on POD modes\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

% 3 -- fit to projected data (basis computed 
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
[w,e,b] = optdmd(xdata,ts,r,imode,opts,[],u);

% reconstructed values
x1 = w*diag(b)*exp(e*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e,evals);
relerr_e = norm(e(indices)-evals)/norm(evals);

fprintf('example 3 --- fitting data projected on POD modes from randsvd2\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

% 4 -- use your own initial guess

% set a random initial guess (not generally a good idea)

e_init = randn(3,1);

imode = 2;
[w,e,b] = optdmd(xdata,ts,r,imode,[],e_init);

% reconstructed values
x1 = w*diag(b)*exp(e*ts);
relerr_r = norm(x1-xdata,'fro')/norm(xdata,'fro');
relerr_r_clean = norm(x1-xclean,'fro')/norm(xclean,'fro');

% compare to actual eigenvalues
indices = match_vectors(e,evals);
relerr_e = norm(e(indices)-evals)/norm(evals);

fprintf('example 4 --- using your own initial guess\n')
fprintf('relative error in reconstruction %e\n',relerr_r)
fprintf('relative error w.r.t clean data %e\n',relerr_r_clean)
fprintf('relative error of eigenvalues %e\n',relerr_e)

