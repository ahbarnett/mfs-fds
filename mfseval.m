function U = mfseval(k,trg,cx,x,meth)
% EVAL.  Sum 2D potential due to MFS sources with given coeffs, at targets
%
% U = eval(k,trg,x,meth)
%  evaluates potential u at m targets trg (2-by-m), using coeff vector x in
%  the Helmholtz MFS with wavenumber k>=0 and source points cx (2-by-N).
% meth = 'd': direct summation
%        'f': fast (FMM) summation

m = size(trg,2);
if meth=='d'
  U = zeros(1,m);
  for i=1:m     % loop over targs
    ri = sqrt((trg(1,i)-cx(1,:)).^2+(trg(2,i)-cx(2,:)).^2);  % row vec
    U(i) = besselh(0,k*ri) * x;    % do the sum
  end
  
elseif meth=='f'
  % ***
  
end