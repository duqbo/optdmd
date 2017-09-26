function opts = varpro_lsqlinopts(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear constraints require MATLAB R2013a
% or later
%
% Create linear constraint options structure
% for varpro routines
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
% of the linear constraints for varpro2. The
% constraints are of the form expected by
% MATLAB's lsqlin routine. The following constraints
% are enforced
%
% Ac*alpha <= bc
% Ace*alpha = bce
% lbc <= alpha <= ubc
%
% If alpha is a complex vector of parameters
% (default behavior) the following constraints
% are enforced
%
% alphar = [real(alpha); imag(alpha)]
% Ac*alphar <= bc
% Ace*alphar = bce
% lbc <= alphar <= ubc
%
% Let ia = length(alpha)
%
% All constraints default to empty arrays.
%
% 'ifreal' (0) --- flag. if ifreal == 1 then
%               the vector alpha should always
%               be real valued (this implies that
%               the matrices Phi and dPhi and the
%               data should be real in varpro2).
%               if ifreal~=1 then the default
%               behavior for complex alpha is assumed.
%
% 'Ac' ([]) --- matrix Ac above. should have
%               length 2*ia if alpha is complex
%               (default) or length ia if alpha
%               is real. 
%
% 'bc' ([]) --- vector bc above, if Ac has m rows
%               bc should be length m
%
% 'Ace' ([]) --- matrix Ace above. should have
%               2*ia columns if alpha is complex
%               (default) or ia columns if alpha
%               is real. 
%
% 'bce' ([]) --- vector bce above, if Ac has m rows
%               bce should be length m
%
% 'lbc' ([]) --- vector lbc above, should be length
%                2*ia for complex alpha (default)
%                or length ia for real alpha.
%
% 'ubc' ([]) --- vector ubc above, should be length
%                2*ia for complex alpha (default)
%                or length ia for real alpha.
%
% 'lsqlinopts' (optimoptions('lsqlin','Algorithm','interior-point',...
%                  'Display','none')) ---
%                set the options for this interior solve.
%                Use a call to optimoptions to create, see
%                MATLAB's LSQLIN documentation for details
%
% See also LSQLIN
  
%
% Author: Travis Askham
%
% Copyright 2017 Travis Askham
%
% License: MIT License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errstr1 = 'nargin = %d. input should be pairs of values';
errstr2 = 'input %d is not a valid field name';
errstr3 = 'input %d is non-numeric';

				% default values
opts.ifreal = 0;
opts.Ac = [];
opts.bc = [];
opts.Ace = [];
opts.bce = [];
opts.lbc = [];
opts.ubc = [];
opts.lsqlinopts = optimoptions('lsqlin','Algorithm',...,
			       'interior-point','Display','none');

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

end
