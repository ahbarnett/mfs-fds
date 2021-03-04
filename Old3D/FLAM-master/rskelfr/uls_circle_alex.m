% Underdetermined least squares on the unit circle, Laplace sources.

% 3/21/18 Barnett tweaked for matlab compat, geom plot, actual MFS in unit disc.
%
% Compare:
% uls_circle(1.5e4,1e4)           % this chooses delta = 5h, typ for MFS
% This gives 6-9 digits, fluctuates (note conds are ~ 1e+15)
% uls_circle(1.5e4,1e4,0.0063)    % 10h,  gets 9 digits reliably.
% But, setting beqAx=0 below gives only 5 digits, out to expected delta <= 10h.
%
% Method needs to handle highly ill-cond case, ie find a small resid,
% but only when there exists such an X, not too large.

function uls_circle(m,n,delta,occ,p,rank_or_tol,rdpiv,store)

  % set default parameters
  if nargin < 1 || isempty(m)
    m = 1.5e4;                         % m a bit larger than n in practice
  end
  if nargin < 2 || isempty(n)
    n = 1e4;
  end
  if nargin < 3 || isempty(delta)
    delta = 5.0 * (2*pi/n)                % 1e-12 is way too small in practice
  end
  if nargin < 4 || isempty(occ)
    occ = 128;
  end
  if nargin < 5 || isempty(p)
    p = 64;
  end
  if nargin < 6 || isempty(rank_or_tol)
    rank_or_tol = 1e-12;
  end
  if nargin < 7 || isempty(rdpiv)
    rdpiv = 'l';
  end
  if nargin < 8 || isempty(store)
    store = 'a';
  end

  % initialize
  theta = (1:m)*2*pi/m;
  rx = [cos(theta); sin(theta)];       % BVP on unit disc
  theta = (1:n)*2*pi/n;
  cx = (1 + delta)*[cos(theta); sin(theta)];  %  src pts outside
  fprintf('%dx%d lin sys; src pts are dist %.3g h_src from bdry\n',...
          m,n,n*delta/(2*pi));
  pic = 0; if pic
    figure; plot(rx(1,:),rx(2,:),'.-'); hold on; plot(cx(1,:),cx(2,:),'+');
    axis equal; title('uls circle test geom');
    legend('rx (colloc pts)','cx (src pts)'); drawnow
  end
  M = size(rx,2);
  N = size(cx,2);
  theta = (1:p)*2*pi/p;
  proxy = 1.5*[cos(theta); sin(theta)];

  % compress matrix using RSKELFR
  Afun = @(i,j)Afun2(i,j,rx,cx);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun2(rc,rx,cx,slf,nbr,l,ctr,proxy);
  opts = struct('rdpiv',rdpiv,'verb',1);
  F = rskelfr(Afun,rx,cx,occ,rank_or_tol,pxyfun,opts);
  w = whos('F');
  fprintf([repmat('-',1,80) '\n'])
  fprintf('mem: %6.2f (MB)\n',w.bytes/1e6)

  % compress matrix using IFMM
  opts = struct('store',store,'verb',1);
  G = ifmm(Afun,rx,cx,2*occ,1e-15,pxyfun,opts);
  w = whos('G');
  fprintf([repmat('-',1,80) '\n'])
  fprintf('mem: %6.2f (MB)\n',w.bytes/1e6)

  % test accuracy using randomized power method
  X = rand(N,1);
  X = X/norm(X);

  % NORM(A - F)/NORM(A)
  tic
  rskelfr_mv(F,X);
  t1 = toc;
  tic
  ifmm_mv(G,X,Afun);
  t2 = toc;
  [e,niter] = snorm(N,@(x)(ifmm_mv(G,x,Afun,'n') - rskelfr_mv(F,x,'n')), ...
                      @(x)(ifmm_mv(G,x,Afun,'c') - rskelfr_mv(F,x,'c')));
  e = e/snorm(N,@(x)(ifmm_mv(G,x,Afun,'n')),@(x)(ifmm_mv(G,x,Afun,'c')));
  fprintf('mv: %10.4e / %4d / %10.4e (s) / %10.4e (s)\n',e,niter,t1,t2)

  % test weak pseudoinverse apply accuracy
  X = rand(M,1);
  X = X/norm(X);
  tic
  rskelfr_sv(F,X);
  t = toc;

  % M >= N (overdetermined): NORM(I - PINV(F)*A)
  if M >= N
    [e,niter] = snorm(N,@(x)(x - rskelfr_sv(F,ifmm_mv(G,x,Afun,'n'),'n')),...
                        @(x)(x - ifmm_mv(G,rskelfr_sv(F,x,'c'),Afun,'c')));

  % M < N (underdetermined): NORM(I - A*PINV(F))
  else
    [e,niter] = snorm(M,@(x)(x - ifmm_mv(G,rskelfr_sv(F,x,'n'),Afun,'n')),...
                        @(x)(x - rskelfr_sv(F,ifmm_mv(G,x,Afun,'c'),'c')));
  end
  fprintf('sv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)

  % example application with simple estimates for residual/norm factors
  nrhs = 16;          % correct interpretation?
  %B = rand(M,nrhs);   % unfair to use random vectors here! (not in range)
  maxfreq = ceil(N/10);  % smooth, such that |Xtrue|~1e8 at delta = 30h.
  Bhat = [randn(maxfreq,nrhs); zeros(M-maxfreq,nrhs)];
  B = real(ifft(Bhat));           % overwrite w/ smooth RHS data...
  beqAx = 1; if beqAx % ...or overwrite with this w/ valid RHS vec:
    Xtrue = randn(N,nrhs); B = ifmm_mv(G,Xtrue,Afun);
  end
  X = rskelfr_sv(F,B);
  C = ifmm_mv(G,X,Afun);
  condl = snorm(M,@(x)rskelfr_mvl(F,x,'n'),@(x)rskelfr_mvl(F,x,'c')) ...
        * snorm(M,@(x)rskelfr_svl(F,x,'n'),@(x)rskelfr_svl(F,x,'c'));
  condu = snorm(N,@(x)rskelfr_mvu(F,x,'n'),@(x)rskelfr_mvu(F,x,'c')) ...
        * snorm(N,@(x)rskelfr_svu(F,x,'n'),@(x)rskelfr_svu(F,x,'c'));
  fprintf('uls: %10.4e resid nrm / %10.4e coeff nrm / %10.4e condL / %10.4e condU\n', norm(B - C)/norm(B),norm(X),condl,condu)
  if beqAx, fprintf('|Xtrue|_2 = %10.4e, |X - Xtrue|_2 = %10.4e (doesn''t need to be small for small resid)\n',norm(Xtrue), norm(X-Xtrue)), end
end
  
% kernel function
function K = Kfun(x,y)
  dx = bsxfun(@minus,x(1,:)',y(1,:));
  dy = bsxfun(@minus,x(2,:)',y(2,:));
  K = -1/(2*pi)*log(sqrt(dx.^2 + dy.^2));
end
  
% matrix entries
function A = Afun2(i,j,rx,cx)
  A = Kfun(rx(:,i),cx(:,j));
end

% proxy function
function [Kpxy,nbr] = pxyfun2(rc,rx,cx,slf,nbr,l,ctr,proxy)
  pxy = bsxfun(@plus,proxy*l,ctr');
  if strcmpi(rc,'r')
    Kpxy = Kfun(rx(:,slf),pxy);
    dx = cx(1,nbr) - ctr(1);
    dy = cx(2,nbr) - ctr(2);
  elseif strcmpi(rc,'c')
    Kpxy = Kfun(pxy,cx(:,slf));
    dx = rx(1,nbr) - ctr(1);
    dy = rx(2,nbr) - ctr(2);
  end
  dist = sqrt(dx.^2 + dy.^2);
  nbr = nbr(dist/l < 1.5);
end