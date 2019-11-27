% Consistent least squares on the unit sphere, Laplace kernel.
%
% This is essentially the 3D version of LS_CIRCLE, where the solution is now
% matched on the surface of the unit sphere and the MFS sources are placed at a
% slight offset from it.

function ls_sphere(m,n,delta,occ,p,rank_or_tol,rdpiv,store,doiter)

  % set default parameters
  if nargin < 1 || isempty(m), m = 16384; end  % number of row points
  if nargin < 2 || isempty(n), n =  8192; end  % number of col points
  if nargin < 3 || isempty(delta), delta = 1e-3; end  % MFS offset
  if nargin < 4 || isempty(occ), occ = 1024; end
  if nargin < 5 || isempty(p), p = 512; end  % number of proxy points
  if nargin < 6 || isempty(rank_or_tol), rank_or_tol = 1e-6; end
  if nargin < 7 || isempty(rdpiv), rdpiv = 'l'; end  % redundant pivoting
  if nargin < 8 || isempty(store), store = 'a'; end  % FMM storage mode
  if nargin < 9 || isempty(doiter), doiter = 1; end  % naive LSQR/CG?

  % initialize
  rx = randn(3,m); rx = (1 + delta)*rx./sqrt(sum(rx.^2));  % row points
  cx = randn(3,n); cx =             cx./sqrt(sum(cx.^2));  % col points
  M = size(rx,2);
  N = size(cx,2);
  % proxy points are quasi-uniform sampling of scaled 1.5-radius sphere
  proxy = trisphere_subdiv(p,'v'); r = randperm(size(proxy,2));
  proxy = proxy(:,r(1:p));  % reference proxy points are for unit box [-1, 1]^3

  % factor matrix using RSKELFM
  Afun = @(i,j)Afun_(i,j,rx,cx);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy);
  opts = struct('rdpiv',rdpiv,'verb',1);
  tic; F = rskelfm(Afun,rx,cx,occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('rskelfm time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % compress matrix using IFMM
  rank_or_tol = max(rank_or_tol*1e-2,1e-15);  % higher accuracy for reference
  opts = struct('store',store);
  tic; G = ifmm(Afun,rx,cx,2*occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('G'); mem = w.bytes/1e6;
  fprintf('ifmm time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % test accuracy using randomized power method
  X = rand(N,1);
  X = X/norm(X);

  % NORM(A - F)/NORM(A)
  tic; rskelfm_mv(F,X); t = toc;  % for timing
  err = snorm(N,@(x)(ifmm_mv(G,x,Afun,'n') - rskelfm_mv(F,x,'n')), ...
                @(x)(ifmm_mv(G,x,Afun,'c') - rskelfm_mv(F,x,'c')));
  err = err/snorm(N,@(x)ifmm_mv(G,x,Afun,'n'),@(x)ifmm_mv(G,x,Afun,'c'));
  fprintf('rskelfm_mv err/time: %10.4e / %10.4e (s)\n',err,t)
  tic; ifmm_mv(G,X,Afun); t = toc;
  fprintf('ifmm_mv time: %10.4e (s)\n',t)

  % test MFS solution
  B = ifmm_mv(G,X,Afun);              % random right-hand side B = A*X
  tic; Y = rskelfm_sv(F,B); t = toc;  % for timing (and concrete example)
  fprintf('ls:\n')
  % residual error: NORM(B - A*F\B)/NORM(B) <~ NORM(A - A*F\A)/NORM(A)
  err = snorm(N, ...
    @(x)(ifmm_mv(G,x,Afun,'n') - ...
         ifmm_mv(G,rskelfm_sv(F,ifmm_mv(G,x,Afun,'n'),'n'),Afun,'n')), ...
    @(x)(ifmm_mv(G,x,Afun,'c') - ...
         ifmm_mv(G,rskelfm_sv(F,ifmm_mv(G,x,Afun,'c'),'c'),Afun,'c')));
  err = err/snorm(N,@(x)ifmm_mv(G,x,Afun,'n'),@(x)ifmm_mv(G,x,Afun,'c'));
  fprintf('  resid err/time: %10.4e / %10.4e (s)\n',err,t)
  % solution error: NORM(X - F\B)/NORM(X) <= NORM(I - F\A)
  err1 = snorm(N, ...
    @(x)(x - rskelfm_sv(F,ifmm_mv(G,x,Afun,'n'))), ...
    @(x)(x - ifmm_mv(G,rskelfm_sv(F,x,'c'),Afun,'c')));
  % solution norm: NORM(F\B) <= NORM(F\A)*NORM(X)
  err2 = snorm(N, ...
    @(x)rskelfm_sv(F,ifmm_mv(G,x,Afun,'n')), ...
    @(x)ifmm_mv(G,rskelfm_sv(F,x,'c'),Afun,'c'));
  fprintf('  soln err/norm: %10.4e / %10.4e\n',err1,err2)
  % concrete example
  err1 = norm(B - ifmm_mv(G,Y,Afun))/norm(B);
  err2 = norm(X - Y);
  fprintf('  ex resid/soln err/soln norm: %10.4e / %10.4e / %10.4e\n', ...
          err1,err2,norm(Y))

  iter = nan;
  if ~isoctave()
    mv = @(x,trans)mv_lsqr(G,x,trans,Afun);

    % run LSQR
    if doiter, [~,~,~,iter] = lsqr(mv,B,1e-9,128); end

    % run LSQR with initial guess
    tic; [Z,~,~,piter] = lsqr(mv,B,1e-9,32,[],[],Y); t = toc;
    fprintf('lsqr:\n')
  else
    warning('No LSQR in Octave.')

    C = ifmm_mv(G,B,Afun,'c');
    mv = @(x)mv_cg(G,x,Afun);

    % run CG (on normal equations)
    if doiter, [~,~,~,iter] = pcg(mv,C,1e-9,128); end

    % run CG with initial guess from pseudoinverse
    tic; [Z,~,~,piter] = pcg(mv,C,1e-9,32,[],[],Y); t = toc;
    fprintf('cg:\n')
  end
  err1 = norm(B - ifmm_mv(G,Z,Afun))/norm(B);
  err2 = norm(X - Z)/norm(X);
  fprintf('  resid/soln err/time: %10.4e / %10.4e / %10.4e (s)\n',err1,err2,t)
  fprintf('  init/uninit iter: %d / %d\n',piter,iter)
end

% kernel function
function K = Kfun(x,y)
  dx = x(1,:)' - y(1,:);
  dy = x(2,:)' - y(2,:);
  dz = x(3,:)' - y(3,:);
  K = 1/(4*pi)./sqrt(dx.^2 + dy.^2 + dz.^2);
end

% matrix entries
function A = Afun_(i,j,rx,cx)
  A = Kfun(rx(:,i),cx(:,j));
end

% proxy function
function [Kpxy,nbr] = pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy)
  pxy = proxy*l + ctr';  % scale and translate reference points
  if rc == 'r'
    Kpxy = Kfun(rx(:,slf),pxy);
    dx = cx(1,nbr) - ctr(1);
    dy = cx(2,nbr) - ctr(2);
    dz = cx(3,nbr) - ctr(3);
  else
    Kpxy = Kfun(pxy,cx(:,slf));
    dx = rx(1,nbr) - ctr(1);
    dy = rx(2,nbr) - ctr(2);
    dz = rx(3,nbr) - ctr(3);
  end
  % proxy points form sphere of scaled radius 1.5 around current box
  % keep among neighbors only those within sphere
  dist = sqrt(dx.^2 + dy.^2 + dz.^2);
  nbr = nbr(dist/l < 1.5);
end

% matrix multiply for LSQR
function y = mv_lsqr(F,x,trans,Afun)
  if     strcmpi(trans,'notransp'), y = ifmm_mv(F,x,Afun,'n');
  elseif strcmpi(trans,  'transp'), y = ifmm_mv(F,x,Afun,'c');
  end
end

% matrix multiply for CG
function y = mv_cg(F,x,Afun)
  y = ifmm_mv(F,ifmm_mv(F,x,Afun,'n'),Afun,'c');
end
