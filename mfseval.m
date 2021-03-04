function U = mfseval(k,trg,cx,x,meth,opts)
% EVAL.  Sum 2D or 3D potential due to MFS sources with given coeffs, at targets
%
% U = eval(k,trg,cx,x,meth)
%  evaluates potential u at m targets trg (d-by-m, where space dimension d
%  inferred from trg), using coeff vector x in
%  the Helmholtz MFS with wavenumber k>=0, and source points cx (d-by-N).
%  Returns a row vector.
%
% Note: normalization for sources differs from math fundamental solution:
%  H_0^{(1)}(kr) in 2D,  e(ikr)/r in 3D.
%
% meth = 'd': direct summation
%        'f': fast (FMM) summation
%
% U = eval(k,trg,cx,x,meth,opts) also controls options such as...
%  opts.iprec controls FMM accuracy
%
% Just values (no derivatives yet).

if nargin<6, opts = []; end

dim = size(trg,1);
if dim<2 || dim>3, error(sprintf('dim = %d should be 2 or 3!',dim)); end
if size(cx,1)~=dim, error('cx and trg dimensions incompatible!'); end
M = size(trg,2);      % num targs
N = size(cx,2);       % num sources

% Note: in following, normalization must match the kernels in factor.m ...
if meth=='d'
  U = zeros(1,M);
  if dim==2
    for i=1:M     % loop over targs
      ri = sqrt((trg(1,i)-cx(1,:)).^2+(trg(2,i)-cx(2,:)).^2);  % row vec
      U(i) = besselh(0,k*ri) * x;    % do the sum
    end
  else
    for i=1:M     % loop over targs
      ri = sqrt((trg(1,i)-cx(1,:)).^2+(trg(2,i)-cx(2,:)).^2+(trg(3,i)-cx(3,:)).^2);  % row vec
      U(i) = (exp(1i*k*ri)./ri) * x;    % do the sum
    end
  end
  
elseif meth=='f'
  if ~isfield(opts,'iprec'), opts.iprec=4; end   % 12 digits by default
  if dim==2
    ifcharge = 1;
    ifdipole = 0;   % note dummy inputs for dipoles here...
    iffldtarg = 0;  % for now
    U = hfmm2dpart(opts.iprec,k,N,cx,ifcharge,x,ifdipole,...
                   0*x,0*cx,0,0,0,M,trg,1,iffldtarg,0);
    U = (4/1i) * U.pottarg;   % undo the hfmm2d normalization. is a row vec
  else              % 3D
    ifcharge = 1;
    ifdipole = 0;   % note dummy inputs for dipoles here...
    iffldtarg = 0;  % for now
    U = hfmm3dpart(opts.iprec,k,N,cx,ifcharge,x,ifdipole,...
                   0*x,0*cx,0,0,M,trg,1,iffldtarg);
    U = U.pottarg;   % undo the hfmm2d normalization. is a row vec
  end
end