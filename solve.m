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

if meth=='l'
  X = F.U\(F.L\B);   % two back-subs
elseif meth=='q'
  X = F.R\(F.Q'*B);
else
  % ***
  
end
