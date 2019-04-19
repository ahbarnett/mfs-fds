function X = solve(F,B,meth)
% SOLVE   solve LSQ lin sys previously factored by dense or FDS.
%
% X = solve(F,B)
%
% Inputs:
%  F = factor struct from factor.m
%  B = (M-by-p) stacked columns of RHS to solve, where A was M-by-N.
%  meth = method string, should match that used in factor.m :
%         'l' dense direct LU
%         'q' dense direct QR
%         'r' FLAM rskel, augmented by I, using lsq params stored in F
%
% Outputs:
%  X = (N-by-p) stacked solution vectors.

% Barnett fixed the LU for rect case, 4/19/19.

t0=tic;
if meth=='l'
  X = F.U\(linsolve(F.L,F.P*B,struct('LT',true)));   % perm the RHS then two back-subs. Note have to tell it the L solve is lower-tri to get O(N^2)
elseif meth=='q'
  X = F.R\(F.Q'*B);  % one mult and one back-sub, slower than LU
elseif meth=='r'
  Ms = size(F.A,1); M = size(B,1);
  B = [F.lsqpar.tau*B; zeros(Ms-M,size(B,2))];
  if F.lsqpar.qr=='q'
    solvefn = @(b)(F.R\(F.R'\(F.A'*b)));
  else  % 's'
    solvefn = @(b)(F.P*(F.R\spqr_qmult(F.Q,b,0)));
  end
  X = solvefn(B);
  for i = 1:F.lsqpar.refine  % iterative refinement
    X = X + solvefn(B - F.A*X);
  end
  X = X(1:F.N,:);
end
fprintf('solve %.3g s\n',toc(t0))

