function opts = varpro_opts(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create options structure for varpro routines
%
% INPUT: 
%
% The input should be pairs of strings and 
% values for setting fields of the structure
%
% OUTPUT:
%
% The output will be a structure with the 
% given values for the specified fields and
% the default values for the rest
%
% OPTIONS STRUCTURE FORMAT
%
% The options structure has several fields,
% denoted as strings below, with default values
% in parentheses, that determine the behavior 
% of varpro2. See the descriptions below
%
% 'lambda0' (1.0) --- lambda0 is the initial 
%   value used for the regularization parameter
%   lambda in the Levenberg method (a larger
%   lambda makes it more like gradient descent)
%
% 'maxlam' (52) --- maxlam is the maximum number 
%   of steps used in the inner Levenberg loop,
%   i.e. the number of times you increase lambda
%   before quitting
%
% 'lamup' (2.0) --- lamup is the factor by which
%   you increase lambda when searching for an 
%   appropriate step
%
% 'lamdown' (2.0) --- lamdown is the factor by which
%   you decrease lambda when checking if that
%   results in an error decrease
%
% 'ifmarq' (1) --- ifmarq is a flag which determines
%   whether you use the Levenberg algorithm or the
%   Levenberg-Marquardt algorithm. ifmarq == 1 
%   results in the Levenberg-Marquardt algorithm.
%   Anything else gives the standard Levenberg
%   algorithm
%
% 'maxiter' (30) --- the maximum number of outer
%   loop iterations to use before quitting
%
% 'tol' (1.0e-6) --- the tolerance for the relative
%   error in the residual, i.e. the program will
%   terminate if 
%       norm(y-Phi(alpha)*b,'fro')/norm(y,'fro') < tol
%   is achieved.
%
% 'eps_stall' (1.0e-12) --- the tolerance for detecting 
%   a stall. If err(iter-1)-err(iter) < eps_stall*err(iter-1)
%   then a stall is detected and the program halts.
%
% 'iffulljac' (1) --- flag determines whether or not to use
%   the full expression for the Jacobian or Kaufman's 
%   approximation.
%
% Author: Travis Askham
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errstr1 = 'nargin = %d. input should be pairs of values';
errstr2 = 'input %d is not a valid field name';
errstr3 = 'input %d is non-numeric';

% default values

opts.lambda0 = 1.0;
opts.maxlam = 52;
opts.lamup = 2.0;
opts.lamdown = 2.0;
opts.ifmarq = 1;
opts.maxiter = 30;
opts.tol = 1.0e-6;
opts.eps_stall = 1.0e-12;
opts.iffulljac = 1;

% minimum values (in some reasonable sense)

optsmin.lambda0 = 0.0;
optsmin.maxlam = 0;
optsmin.lamup = 1.0;
optsmin.lamdown = 1.0;
optsmin.ifmarq = -Inf;
optsmin.maxiter = 0;
optsmin.tol = 0.0;
optsmin.eps_stall = -Inf;
optsmin.iffulljac = -Inf;

% maximum values (in some reasonable sense)

optsmax.lambda0 = 1.0e16;
optsmax.maxlam = 200;
optsmax.lamup = 1.0e16;
optsmax.lamdown = 1.0e16;
optsmax.ifmarq = Inf;
optsmax.maxiter = 1.0e12;
optsmax.tol = 1.0e16;
optsmax.eps_stall = 1.0;
optsmax.iffulljac = Inf;

% check if input comes in pairs

if(mod(nargin,2) ~= 0)
    error(errstr1,nargin);
end

% for each pair of inputs, check if the first 
% input is a valid structure field name and if 
% the second input is numeric. If so, set value of 
% that field to the second input. Otherwise, bomb
% with error message.

for i = 1:nargin/2
    if (isfield(opts,varargin{2*(i-1)+1}))
        if (isnumeric(varargin{2*i}))
            opts.(varargin{2*(i-1)+1}) = varargin{2*i};
        else
            error(errstr3,2*i)
        end
    else
        error(errstr2,2*(i-1)+1);
    end
end

varpro_opts_warn(opts,optsmin,optsmax);

end

function varpro_opts_warn(opts,optsmin,optsmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For each value of the opts structure, this 
% routine prints a warning if it is not within
% the bounds determined by optsmin and optsmax
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warnstrmin = ['option %s with value %e is less than %e, which' ...
    ' is not recommended'];
warnstrmax = ['option %s with value %e is greater than %e, which' ...
    ' is not recommended'];

names = fieldnames(opts);

for i = 1:length(names)
    s = names{i};
    optv = opts.(s);
    optminv = optsmin.(s);
    optmaxv = optsmax.(s);
    if (optv < optminv)
        warning(warnstrmin,s,optv,optminv);
    end
    if (optv > optmaxv)
        warning(warnstrmax,s,optv,optmaxv);        
    end
end
end
