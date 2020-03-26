function F = factor(k,rx,cx,meth,flampar,lsqpar)
% FACTOR  fast factor LSQ system, for 2D or 3D PDE kernel, MFS plus rskel FDS.
%
% F = factor(k,rx,cx,meth,flampar,lsqpar)
%
% Inputs:
%  k = wavenumber >= 0
%  rx = d-by-M surface pts,  where dim d = 2 or 3 (determined from coords)
%  cx = d-by-N source (MFS) pts
%  meth = factorization method string, should match that in solve.m :
%         'l' dense direct LU,   O(N^3)
%         'q' dense direct QR,   O(N^3) and slower than LU, but more stable?
%         'r' FLAM rskel, augmented by I.    O(N) unless space-filling.
%  flampar = FLAM parameter struct (only used if meth='r')
%  lsqpar = LSQ parameter struct (ditto)
%
% Note: normalization for MFS representation differs from fundamental solution:
%   H_0^{(1)}(kr) in 2D,  e(ikr)/r in 3D.
%
% Output:
%  F = struct containing all factor info, and lsq params if needed

dim = size(rx,1);
if dim<2 || dim>3, error(sprintf('dim = %d should be 2 or 3!',dim)); end

if meth=='l' || meth=='q'  % ---------------- dense meths
  tic; A = Kfun(rx,cx,k);
  w = whos('A'); fprintf('A (%dx%d) fill %.3g s \t (%.0f MB)\n',size(rx,2),size(cx,2),toc,w.bytes/1e6)
  tic;
  if meth=='l'
    [F.L,F.U,F.P] = lu(A,'matrix');     % F.P is the perm, as a matrix
    fprintf('dense LU %.3g s\n',toc)
  else  % 'q'
    [F.Q,F.R] = qr(A,0);      % econ
    fprintf('dense QR %.3g s\n',toc)
  end
  F.A = A;   % save it for kicks

else                       % ---------------- FDS
  Afun = @(i,j)Kfun(rx(:,i),cx(:,j),k);
  pxypts = proxypts(dim,flampar.p);  % proxy pts vs unit box
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)proxyfun(rc,rx,cx,slf,nbr,l,ctr,pxypts,k);
  t0=tic; RF = rskel(Afun,rx,cx,flampar.occ,flampar.rank_or_tol,pxyfun,flampar.opts);  % compress
  w = whos('RF'); fprintf('rskel: %.3g s \t %.0f (MB)\n',toc(t0),w.bytes/1e6)

  t0=tic; [A,p,q] = rskel_xsp(RF);  % make a big sparse mat
  w = whos('A'); fprintf('xsp: %.3g s \t %.0f (MB)\n',toc(t0),w.bytes/1e6)
  M = size(rx,2); N = size(cx,2);
  if ~isfield(lsqpar,'meth'), lsqpar.meth='u'; end  % ULS by default
  if lsqpar.meth=='u'
    % ULS: min. norm(x) s.t. A*x = b
    A = [lsqpar.tau*A; speye(N) sparse(N,size(A,2)-N)];
  else  % 'o'
    if ~isfield(lsqpar,'lambda'), lsqpar.lambda=0; end  % no regularization by default
    % OLS: min. norm(A*x - b)^2 + lambda^2*norm(x)  [s.t. to embedding identities]
    A = [lsqpar.tau*A(M+1:end,:); A(1:M,:); lsqpar.lambda*speye(N) sparse(N,size(A,2)-N)];
  end
  t0=tic;
  if lsqpar.qr=='q'
    F.R = qr(A,0);
  else  % 's'
    opts = struct('tol',0,'econ',0,'Q','Householder');
    [F.Q,F.R,F.P] = spqr(A,opts);
  end
  w = whos('F'); fprintf('factor: %.3g s \t %0.f (MB)\n',toc(t0),w.bytes/1e6)
  F.lsqpar = lsqpar;  % store for solve
  F.N = N;
  F.A = A;            % the sparse mat, needed for normal eqns & iter refine.
  F.p = p;            % row perm
  F.q = q;            % col perm
end

% ------------- helper routines for matrix elements in FLAM --------

% kernel function (see doc note about normalization)
function K = Kfun(rx,cx,k)
if size(rx,1)==2
  r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2);
  K = besselh(0,k*r);
else
  r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2 + (rx(3,:)'-cx(3,:)).^2);
  K = exp(1i*k*r)./r;
end

% proxy points for FDS -- circle/sphere around unit box
function pts = proxypts(dim,p)
if dim==2
  theta = (1:p)*2*pi/p;
  pts = 1.5*[cos(theta); sin(theta)];
else         % 3D
  error('3d proxypts not yet supported')
end

function [Kpxy,nbr] = proxyfun(rc,rx,cx,slf,nbr,l,ctr,pxypts,k)
% proxy function for FDS: outputs kernel matrix from pts <-> proxies
% for either row or col points, and outputs subselected nbrs to keep.
pts = pxypts.*l + ctr;       % scale and center the proxy circle about the box
if strcmpi(rc,'r')
  Kpxy = Kfun(rx(:,slf),pts,k);
  dr = cx(:,nbr) - ctr;
else
  Kpxy = Kfun(pts,cx(:,slf),k);
  dr = rx(:,nbr) - ctr;
end
nbr = nbr(sum((dr./l).^2) < 1.5^2);  % keep only nbr pts inside proxy
