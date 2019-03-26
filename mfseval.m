function U = mfseval(k,trg,cx,x,meth,opts)
% EVAL.  Sum 2D or 3D potential due to MFS sources with given coeffs, at targets
%
% U = eval(k,trg,cx,x,meth)
%  evaluates potential u at m targets trg (d-by-m, where space dimension d
%  inferred from trg), using coeff vector x in
%  the Helmholtz MFS with wavenumber k>=0 and source points cx (d-by-N).
%  Returns a row vector.
% meth = 'd': direct summation
%        'f': fast (FMM) summation
% opts.iprec controls FMM accuracy

if nargin<6, opts = []; end
if ~isfield(opts,'iprec'), opts.iprec=4; end

dim = size(trg,1);
if dim<2 || dim>3, error(sprintf('dim = %d should be 2 or 3!',dim)); end
if size(cx,1)~=dim, error('cx and trg dimensions incompatible!'); end
M = size(trg,2);
N = size(cx,2);

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
  % ***
  
end