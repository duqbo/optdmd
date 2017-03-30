
function varargout = XORMQR(SIDE,TRANS,A,tau,B,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A wrapper for calling lapack routines 
% for multiplying by Householder reflectors,
% as computed by XGEQP3
%
% We utilize the lapack.m file 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(A);
[mb,nb] = size(B);

if ( SIDE ~= 'L' && SIDE ~= 'R' )
    error('input SIDE must be equal to one of L or R')
end
if ( TRANS ~= 'N' && TRANS ~= 'T' && TRANS ~= 'C')
    error('input TRANS must be equal to one of N, T, or C')
end

if (nargin > 5)
    k = varargin{1};
else
    k = min(m,n);
end

if (isa(A,'double'))
    if (isreal(A))
        s = 'D';
        funstr = [s, 'ORMQR'];
        lwork = -1;
        C = lapack(funstr,SIDE,TRANS,mb,nb,k,A,m,tau, ...
            B,mb,zeros(10,1),lwork,0);
        lwork = C{11}(1);
        C = lapack(funstr,SIDE,TRANS,mb,nb,k,A,m,tau, ...
            B,mb,zeros(lwork,1),lwork,0);
    else
        s = 'Z';
        funstr = [s, 'UNMQR'];
        if (TRANS == 'T')
            TRANS = 'C';
        end
        lwork = -1;
        C = lapack(funstr,SIDE,TRANS,mb,nb,k,A,m,tau, ...
            B,mb,zeros(10,1),lwork,0);
        lwork = C{11}(1);
        C = lapack(funstr,SIDE,TRANS,mb,nb,k,A,m,tau, ...
            B,mb,zeros(2*lwork,1),lwork,0);
    end
else
    error('Data type of A unsupported at this time');
end

varargout{1} = C{9};










