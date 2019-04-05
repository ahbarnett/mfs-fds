function X = solve(F,B,meth)
% SOLVE   solve LSQ lin sys previously factored by dense or FDS.
%
% X = solve(F,B)
%
% Inputs:
%  F = factor struct from factor.m
%  B = (M-by-p) stacked columns of RHS to solve, where A was M-by-N.
%  meth = 'l' dense direct LU,  'r' FLAM rskel, augmented by I.
%
% Outputs:
%  X = (N-by-p) stacked solution vectors.

t0=tic;
if meth=='l'
  X = F.U\(F.L\B);   % two back-subs
elseif meth=='q'
  X = F.R\(F.Q'*B);
elseif meth=='r'
  Ms = size(F.A,1); M = size(B,1);
  B = [F.lsqpar.tau*B; zeros(Ms-M,size(B,2))];
  if F.lsqpar.qr=='q'
    X = F.R\(F.R'\(F.A'*B));  % normal eqns with iterative refinement
    for i = 1:F.lsqpar.refine
      X = X + F.R\(F.R'\(F.A'*(B - F.A*X)));
    end
  else  % 's'
  end
  X = X(1:F.N,:);
end
fprintf('solve %.3g s\n',toc(t0))

