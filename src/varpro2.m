function [b,alpha,niter,err,imode,alphas] = varpro2(y,t,phi,dphi, ...
    m,n,is,ia,alpha_init,opts,varargin)
%VARPRO2 Variable projection algorithm for multivariate data
%
% Attempts a fit of the columns of y as linear combinations
% of the columns of phi(alpha,t), i.e.
%
% y_k = sum_j=1^n b_jk phi_j(alpha,t)
%
% Note that phi(alpha,t) is a matrix of dimension
% m x n where m is length (t) and n is number of columns.
%
% phi_j(alpha,t) is the jth column
% y_k is the kth column of the data
%
% Input:
%
% y - M x IS matrix of data
% t - M vector of sample times
% phi(alpha,t) - M x N matrix (or sparse matrix) valued 
%              function with input alpha
% dphi(alpha,t,i) - M x N matrix (or sparse matrix) valued
%                 function of alpha: returns the derivative 
%                 of the entries of phi with respect to the 
%                 ith component of alpha
% m - integer, number of rows of data/number of sample times
% n - integer, number of columns of phi
% is - integer, number of columns of data .. number of 
%      functions to fit
% ia - integer, dimension of alpha
% alpha_init - initial guess for alpha
% opts - options structure. See varpro_opts.m for details. Can
%   be created with default values via 
%       opts = varpro_opts();
%
% varargin{1} = copts - linear constraint options structure.
%                       See varpro_lsqlinopts.m for details.
%                       allows you to enforce linear constraints.
%                       (these are treated as a projection
%                       operator applied at each step)
%
% the constrained version uses lsqlin from the optimization
% toolbox. the bounds enforced are
%
% Ac*alpha <= bc
% Ace*alpha = bce
% lbc <= alpha <= ubc
%
% For unbounded parts, set to lbc(i) = -Inf or ubc(j) = Inf
% For unneeded bounds, set input to []
% see LSQLIN for more info
%
% NOTE: Linear constraints require MATLAB R2013a
% or later
% 
% varargin{2} = gamma - tikhonov regularization term. if provided
%                       the minimization problem becomes 
%
%                min  | y - phi*b |_F^2 + | gamma alpha |_2^2 
%
%               where gamma is either a scalar or matrix.      
%
% varargin{3} = prox - proximal operator to be applied to the 
%                      vector alpha at each step (e.g. projection
%                      onto a set).
%
%
% Output:
%
% b - N x IS matrix of coefficients .. each column gives
%     the coefficients for one of the functions (columns
%     of data) corresponding to the best fit
% alpha - N vector of values of alpha for best fit
% niter - number of iterations of the Marquardt algorithm
% err - the error for each iteration of the algorithm
% imode - failure mode
%            imode = 0, normal execution, tolerance reached
%            imode = 1, maxiter reached before tolerance
%            imode = 4, failed to find new search direction
%                       at step niter
%
% Example:
%
%   >> [b,alpha,niter,err,imode,alphas] = varpro2(y,t,phi,dphi, ...
%    m,n,is,ia,alpha_init,opts)
%
% See also VARPRO_OPTS, VARPRO_LSQLINOPTS, LSQLIN

%
% Copyright (c) 2018 Travis Askham
%
% Available under the MIT license
%
% References: 
% - Extensions and Uses of the Variable Projection 
% Algorithm for Solving Nonlinear Least Squares Problems by 
% G. H. Golub and R. J. LeVeque ARO Report 79-3, Proceedings 
% of the 1979 Army Numerical Analsysis and Computers Conference
% - "Variable projection for nonlinear least squares problems." 
% Computational Optimization and Applications 54.3 (2013): 579-593. 
% by Dianne P. O’Leary and Bert W. Rust. 
%


% various error and warning string formats

mode8str = ['stall detected: residual reduced by less than %e' ...
    ' times residual at previous step. iteration %d' ...
    '. current residual %e\n'];
mode4str = ['failed to find appropriate step length at iteration %d' ...
    '. current residual %e\n'];
mode1str = ['failed to reach tolerance after maxiter = %d' ...
    ' iterations. current residual %e\n'];
optstr1 = ['input opts must be a structure, see varpro_opts.m.' ...
    ' Using default options ...\n'];
    
% set options, try to catch mistakes

opts_default = varpro_opts();

lambda0 = varpro2_getfield(opts,opts_default,'lambda0');
maxlam = varpro2_getfield(opts,opts_default,'maxlam');
lamup = varpro2_getfield(opts,opts_default,'lamup');
lamdown = varpro2_getfield(opts,opts_default,'lamdown');
ifmarq = varpro2_getfield(opts,opts_default,'ifmarq');
maxiter = varpro2_getfield(opts,opts_default,'maxiter');
tol = varpro2_getfield(opts,opts_default,'tol');
eps_stall = varpro2_getfield(opts,opts_default,'eps_stall');
iffulljac = varpro2_getfield(opts,opts_default,'iffulljac');
ifprint = varpro2_getfield(opts,opts_default,'ifprint');
ptf = varpro2_getfield(opts,opts_default,'ptf');

if (~isstruct(opts))
  if(ifprint == 1)
    fprintf(optstr1);
  end
  opts = opts_default;
end

% if linear constraints are on, get them

ifproxfun = 0;

if (nargin > 10 && ~isempty(varargin{1}))

  ifproxfun = 1;
  
  copts = varargin{1};
  copts_default = varpro_lsqlinopts();

  Ac = varpro2_getfield(copts,copts_default,'Ac');
  bc = varpro2_getfield(copts,copts_default,'bc');
  Ace = varpro2_getfield(copts,copts_default,'Ace');
  bce = varpro2_getfield(copts,copts_default,'bce');
  lbc = varpro2_getfield(copts,copts_default,'lbc');
  ubc = varpro2_getfield(copts,copts_default,'ubc');
  ifreal = varpro2_getfield(copts,copts_default,'ifreal');
  lsqlinopts = varpro2_getfield(copts,copts_default,'lsqlinopts');
  
  % force initial guess inside
  jpvt = 1:length(alpha_init);
  %[delta,ier] = varpro2_lsqlin(eye(length(alpha_init)),alpha_init, ...
  %    alpha_init,jpvt,Ac,bc,Ace,bce,lbc,ubc,ifreal,lsqlinopts);
  %alpha_init = alpha_init+delta;
    %  ier
  proxfun = @(alpha) varpro2_lsqlin_prox(alpha,jpvt,Ac,bc, ...
        Ace,bce,lbc,ubc,ifreal,lsqlinopts);
  alpha_init = proxfun(alpha_init);
  
end	

% if Tikhonov regularization is on, get it		 

if (nargin > 11 && ~isempty(varargin{2}))

  iftik = 1;

  gamma =  varargin{2};
  [mg,ng] = size(gamma);

  if (mg == 1 && ng == 1)
    gamma = gamma*eye(ia);
  elseif (mg ~= ia || ng ~= ia)
    error('Tikhonov regularization matrix of incorrect size');
  end
  
else
    iftik = 0;
    gamma = zeros(ia);
end

if (nargin > 12 && ~isempty(varargin{3}))
    ifproxfun = 1;
    proxfun = varargin{3};
    alpha_init = proxfun(alpha_init);
end

% initialize values

alpha = alpha_init;
alphas = zeros(length(alpha),maxiter);
if (iftik == 1)
  djacmat = zeros(m*is+ia,ia);
  rhstemp = zeros(m*is+ia,1);
else
  djacmat = zeros(m*is,ia);
  rhstemp = zeros(m*is,1);
end
err = zeros(maxiter,1);
res_scale = norm(y,'fro');
scales = zeros(ia,1);

rjac = zeros(2*ia,ia);

phimat = phi(alpha,t);
[U,S,V] = svd(phimat,'econ');
sd = diag(S);
tolrank = m*eps;
irank = sum(sd > tolrank*sd(1));
U = U(:,1:irank);
S = S(1:irank,1:irank);
V = V(:,1:irank);
b = phimat\y;
res = y - phimat*b;
errlast = 0.5*(norm(res,'fro')^2 + norm(gamma*alpha)^2);

imode = 0;

for iter = 1:maxiter
  
		   % build jacobian matrix, looping over alpha indeces
  
  for j = 1:ia
    dphitemp = dphi(alpha,t,j);
    djaca = (dphitemp - sparse(U*(sparse(U'*dphitemp))))*b;
    if (iffulljac == 1)
				% use full expression for Jacobian
      djacb = U*(S\(V'*(sparse(dphitemp'*res))));
      djacmat(1:m*is,j) = (djaca(:) + djacb(:));
    else
				% use approximate expression
      djacmat(1:m*is,j) = djaca(:);
    end
		   % the scales give the "marquardt" part of the algo.
    scales(j) = 1;
    if (ifmarq == 1)
      scales(j) = min(norm(djacmat(1:m*is,j)),1);
      scales(j) = max(scales(j),1e-6);
    end
  end

  if (iftik == 1)
     djacmat(m*is+1:end,:) = gamma;
  end

	% loop to determine lambda (lambda gives the "levenberg" part)

			% pre-compute components that don't depend on 
			% step-size parameter (lambda)
  
		  % get pivots and lapack style qr for jacobian matrix
  rhstemp(1:m*is) = res(:);
  if (iftik == 1)
     rhstemp(m*is+1:end) = -gamma*alpha;
  end

  g = djacmat'*rhstemp;
  
  [qout,djacout,jpvt] = qr(djacmat,0);
  %[djacout,jpvt,tau] = xgeqp3_m(djacmat);
  ijpvt = 1:ia;
  ijpvt(jpvt) = ijpvt;
  rjac(1:ia,:) = triu(djacout(1:ia,:));
  %rhstop = xormqr_m('L','T',djacout,tau,rhstemp); % Q'*res
  rhstop = qout'*rhstemp;
  scalespvt = scales(jpvt(1:ia)); % permute scales appropriately...
  rhs = [rhstop(1:ia); zeros(ia,1)]; % transformed right hand side
  
		  % check if current step size or shrunk version works
  
				% get step

  rjac(ia+1:2*ia,:) = lambda0*diag(scalespvt);

  delta0 = rjac\rhs;
  delta0 = delta0(ijpvt); % unscramble solution
  
				% new alpha guess
  
  if (ifproxfun == 1)
    alpha0 = proxfun(alpha + delta0);
    delta0 = alpha0-alpha; % update delta0
  else
    alpha0 = alpha + delta0;
    
  end
				% corresponding residual
  
  phimat = phi(alpha0,t);
  b0 = phimat\y;
  res0 = y-phimat*b0;
  err0 = 0.5*(norm(res0,'fro')^2 + norm(gamma*alpha0)^2);
  
				% new update rule: check
			 % predicted improvement vs actual improvement
  act_impr = errlast-err0;
  pred_impr = real(0.5*delta0'*(g));
  impr_rat = act_impr/pred_impr;
  
  if (err0 < errlast)
				% new version
				% rescale lambda based on
				% actual vs pred improvement
      
      lambda0 = lambda0*max(1.0/3.0,1-(2*impr_rat-1)^3);
      alpha = alpha0;
      errlast = err0;
      b = b0;
      res = res0;
    
  else
	% if not, increase lambda until something works
	% this makes the algorithm more and more like gradient descent
    
    for j = 1:maxlam
      
      lambda0 = lambda0*lamup;
      rjac(ia+1:2*ia,:) = lambda0*diag(scalespvt);

      delta0 = rjac\rhs;
      delta0 = delta0(ijpvt); % unscramble solution
      
      alpha0 = alpha + delta0;
      if (ifproxfun== 1)
        alpha0 = proxfun(alpha0);
        delta0 = alpha0-alpha;
      end


      phimat = phi(alpha0,t);
      b0 = phimat\y;
      res0 = y-phimat*b0;
      err0 = 0.5*(norm(res0,'fro')^2+norm(gamma*alpha0)^2);
      
      if (err0 < errlast) 
        break
      end

    end
    
    if (err0 < errlast) 
      alpha = alpha0;
      errlast = err0;
      b = b0;
      res = res0;
    else
      
				% no appropriate step length found
      
      niter = iter;
      err(niter) = errlast;
      imode = 4;
      if (ifprint == 1)
	fprintf(mode4str,iter,errlast);
      end
      return
    end
  end
  
  alphas(:,iter) = alpha;
  
  err(iter) = errlast;

  if (ifprint == 1 && mod(iter,ptf) == 0)
    fprintf('step %d err %e lambda %e\n',iter,errlast,lambda0)
  end
  
  if (errlast < tol)
    
				% tolerance met
    
    niter = iter;
    return;
  end
  
  if (iter > 1)
    if (err(iter-1)-err(iter) < eps_stall*err(iter-1))
      
				% stall detected
      
      niter = iter;
      imode = 8;
      if(ifprint == 1)
	fprintf(mode8str,eps_stall,iter,errlast);
      end
      return;
    end
  end
  
  phimat = phi(alpha,t);
  [U,S,V] = svd(phimat,'econ');
  sd = diag(S);
  irank = sum(sd > tolrank*sd(1));
  U = U(:,1:irank);
  S = S(1:irank,1:irank);
  V = V(:,1:irank);
  
end

			   % failed to meet tolerance in maxiter steps

niter = maxiter;
imode = 1;
if (ifprint == 1)
  fprintf(mode1str,maxiter,errlast);
end

end

function out = varpro2_getfield(opts,opts_default,in)
%VARPRO2_GETFIELD Get value of field from struct if it exists,
% otherwise set to default value

optstr2 = 'opts struct is missing %s field, using default\n';
optstr3 = 'opts default struct is missing %s field! bomb\n';

if (isfield(opts,in))
    out = opts.(in);
else
    fprintf(optstr2,in);
    if (isfield(opts_default,in))
        out = opts_default.(in);
    else
        error(optstr3,in);
    end
end

end

function [delta,ier] = varpro2_lsqlin(rjac,rhs,alpha,jpvt,Ac,bc, ...
				      Ace,bce,lbc,ubc,ifreal,opts)

  ia = length(alpha);
  jpvt =jpvt(1:ia);

  ier = 0;
  delta = zeros(size(alpha));
  
  if (ifreal == 1)

    if (~isreal(alpha) || ~isreal(rjac) || ~isreal(rhs))
      ier = 1;
      return
    end
    
    if (~isempty(Ac))
      bc = bc-Ac*alpha;
      Ac = Ac(:,jpvt);
    end
    
    if (~isempty(Ace))
      bce = bce-Ace*alpha;
      Ace = Ace(:,jpvt);
    end

    if (~isempty(lbc))
      lbc = lbc-alpha;
      lbc = lbc(jpvt);
    end
    
    if (~isempty(ubc))
      ubc = ubc-alpha;
      ubc = ubc(jpvt);
    end

    delta = lsqlin(rjac,rhs,Ac,bc,Ace,bce,lbc,ubc,[],opts);
  else
  
    rjacr = [real(rjac), -imag(rjac); imag(rjac), real(rjac)];
    rhsr = [real(rhs); imag(rhs)];
    alphar = [real(alpha); imag(alpha)];

    if (~isempty(Ac))
      bc = bc-Ac*alphar;
      Ac(:,1:ia) = Ac(:,jpvt);
      Ac(:,ia+1:2*ia) = Ac(:,ia+jpvt);
    end
    
    if (~isempty(Ace))
      bce = bce-Ace*alphar;
      Ace(:,1:ia) = Ace(:,jpvt);
      Ace(:,ia+1:2*ia) = Ace(:,ia+jpvt);      
    end

    if (~isempty(lbc))
      lbc = lbc-alphar;
      lbc(1:ia) = lbc(jpvt);
      lbc(ia+1:2*ia) = lbc(ia+jpvt);
    end
    
    if (~isempty(ubc))
      ubc = ubc-alphar;
      ubc(1:ia) = ubc(jpvt);
      ubc(ia+1:2*ia) = ubc(ia+jpvt);
    end
    
    deltar = lsqlin(rjacr,rhsr,Ac,bc,Ace,bce,lbc,ubc,[],opts);
    delta = deltar(1:end/2)+1i*deltar(end/2+1:end);

  end
    
end

function [alphanew,ier] = varpro2_lsqlin_prox(alpha,jpvt,Ac,bc, ...
				      Ace,bce,lbc,ubc,ifreal,opts)

  ia = length(alpha);
  jpvt =jpvt(1:ia);

  ier = 0;
  alphanew = alpha;
  
  if (ifreal == 1)

    if (~isreal(alpha) || ~isreal(rjac) || ~isreal(rhs))
      ier = 1;
      return
    end
    
    if (~isempty(Ac))
      Ac = Ac(:,jpvt);
    end
    
    if (~isempty(Ace))
      Ace = Ace(:,jpvt);
    end

    if (~isempty(lbc))
      lbc = lbc(jpvt);
    end
    
    if (~isempty(ubc))
      ubc = ubc(jpvt);
    end

    mat = eye(ia);
    
    alphanew = lsqlin(mat,alpha,Ac,bc,Ace,bce,lbc,ubc,[],opts);
  else
  
    alphar = [real(alpha); imag(alpha)];
    ia2 = 2*ia;
    mat = eye(ia2);
    
    if (~isempty(Ac))
      Ac(:,1:ia) = Ac(:,jpvt);
      Ac(:,ia+1:2*ia) = Ac(:,ia+jpvt);
    end
    
    if (~isempty(Ace))
      Ace(:,1:ia) = Ace(:,jpvt);
      Ace(:,ia+1:2*ia) = Ace(:,ia+jpvt);      
    end

    if (~isempty(lbc))
      lbc(1:ia) = lbc(jpvt);
      lbc(ia+1:2*ia) = lbc(ia+jpvt);
    end
    
    if (~isempty(ubc))
      ubc(1:ia) = ubc(jpvt);
      ubc(ia+1:2*ia) = ubc(ia+jpvt);
    end
    
    alphanewr = lsqlin(mat,alphar,Ac,bc,Ace,bce,lbc,ubc,[],opts);
    alphanew= alphanewr(1:end/2)+1i*alphanewr(end/2+1:end);

  end
    
end

