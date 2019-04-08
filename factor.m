function F = factor(k,rx,cx,meth,flampar,lsqpar)
% FACTOR  fast factor LSQ system, for 2D or 3D PDE kernel, MFS plus rskel FDS.
%
% F = factor(k,rx,cx,meth,flampar,lsqpar)
%
% Inputs:
%  k = wavenumber >= 0
%  rx = d-by-M surface pts,  where dim d = 2 or 3 (determined from coords)
%  cx = d-by-N source (MFS) pts
%  meth = 'l' dense direct LU
%         'q' dense direct QR
%         'r' FLAM rskel, augmented by I.
%  flampar = FLAM parameter struct
%  lsqpar = LSQ parameter struct
%
% Note: normalization for MFS representation differs from fundamental solution:
%   H_0^{(1)}(kr) in 2D,  e(ikr)/r in 3D.
%
% Output:
%  F = struct containing all factor info.

dim = size(rx,1);
if dim<2 || dim>3, error(sprintf('dim = %d should be 2 or 3!',dim)); end

if meth=='l' || meth=='q'  % ---------------- dense
  tic; A = Kfun(rx,cx,k);
  w = whos('A'); fprintf('A fill %.3g s \t (%.0f MB)\n',toc,w.bytes/1e6)
  tic;
  if meth=='l'
    [F.L,F.U] = lu(A);
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

  t0=tic; A = rskel_xsp(RF);  % make a big sparse mat
  w = whos('A'); fprintf('xsp: %.3g s \t %.0f (MB)\n',toc(t0),w.bytes/1e6)
  N = size(cx,2);
  A = [lsqpar.tau*A; speye(N) sparse(N,size(A,2)-N)];  % treat as underdetermined
  t0=tic;
  if lsqpar.qr=='q'
    F.R = qr(A,0);
  else  % 's'
    opts = struct('Q','Householder');
    [F.Q,F.R,F.P] = spqr(A,opts);
  end
  w = whos('F'); fprintf('factor: %.3g s \t %0.f (MB)\n',toc(t0),w.bytes/1e6)
  F.lsqpar = lsqpar;  % store for solve
  F.N = N;
  F.A = A;
end
end  % main

% kernel function (see doc note about normalization)
function K = Kfun(rx,cx,k)
  if size(rx,1)==2
    r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2);
    K = besselh(0,k*r);
  else
    r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2 + (rx(3,:)'-cx(3,:)).^2);
    K = exp(1i*k*r)./r;
  end
end

% proxy points for FDS -- circle/sphere around unit box
function pts = proxypts(dim,p)
  if dim==2
    theta = (1:p)*2*pi/p;
    pts = 1.5*[cos(theta); sin(theta)];
  else         % 3D
    error('3d proxypts not yet supported')
  end
end

function [Kpxy,nbr] = proxyfun(rc,rx,cx,slf,nbr,l,ctr,pxypts,k)
% proxy function for FDS: outputs kernel matrix from pts <-> proxies
% for either row or col points, and outputs boolean for which pts to keep.
pts = pxypts*l + ctr';       % scale and center the proxy circle about the box
  if strcmpi(rc,'r')
    Kpxy = Kfun(rx(:,slf),pts,k);
    dx = cx(1,nbr) - ctr(1);
    dy = cx(2,nbr) - ctr(2);
  else
    Kpxy = Kfun(pts,cx(:,slf),k);
    dx = rx(1,nbr) - ctr(1);
    dy = rx(2,nbr) - ctr(2);
  end
  dist = sqrt(dx.^2 + dy.^2);
  nbr = nbr(dist/l < 1.5);  % keep only nbr pts inside proxy
end
