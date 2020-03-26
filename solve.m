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
  if F.lsqpar.qr=='q'
    lsfun = @(b)(qr_ls(F,b));
  else  % 's'
    lsfun = @(b)(spqr_ls(F,b));
  end
  Ms = size(F.A,1); [M,p] = size(B);
  if F.lsqpar.meth=='u'
    nc = Ms - F.N;
    lsB = zeros(F.N,p);
    lsD = [B(F.p,:); zeros(nc-M,p)];
  else  % 'o'
    nc = Ms - M - F.N;
    lsB = [B(F.p,:); zeros(F.N,p)];
    lsD = zeros(nc,p);
  end
  warning('off','FLAM:lsedc:maxIterCount')  % disable max iter warning
  X = lsedc(lsfun,F.A(nc+1:end,:),lsB,F.A(1:nc,:)/F.lsqpar.tau,lsD,F.lsqpar.tau,0,F.lsqpar.refine);
  X = X(1:F.N,:);
  X(F.q,:) = X;
end
fprintf('solve %.3g s\n',toc(t0))

function x = qr_ls(F,b)
x = F.R\(F.R'\(F.A'*b));
x = x + F.R\(F.R'\(F.A'*(b - F.A*x)));  % one step of iterative refinement for digit loss from normal equations

function x = spqr_ls(F,b)
x = spqr_qmult(F.Q,b,0);
x = x(1:size(F.R,1),:);
x = F.P*(F.R\x);
