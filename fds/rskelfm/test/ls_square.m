% Consistent least squares on the unit square, regularized Laplace kernel.
%
% This is like LS_CIRCLE except now we match the solution on the unit square
% with regularized Laplace sources over the same geometry. This setup is
% somewhat atypical for MFS problems but is useful for performance benchmarking
% in 2D.
%
% Inputs (defaults are used if not provided or set empty):
%
%   - M: number of row points (default: M = 16384)
%   - N: number of column points (default: N = 8192)
%   - DELTA: MFS offset (default: DELTA = 1e-3)
%   - OCC: tree occupancy parameter (default: OCC = 256)
%   - P: number of proxy points (default: P = 64)
%   - RANK_OR_TOL: local precision parameter (default: RANK_OR_TOL = 1e-6)
%   - TMAX: ID interpolation matrix entry bound (default: TMAX = 2)
%   - RRATIO: rank ratio for rectangular preprocessing (default: RRATIO = 2)
%   - RDPIV: redundant pivoting method (default: RDPIV = 'L')
%   - FASTSV: fast solve mode (default: FASTSV = 'N')
%   - STORE: FMM storage mode (default: STORE = 'A')
%   - DOITER: whether to run naive LSQR/CG (default: DOITER = 1)

function ls_square(M,N,delta,occ,p,rank_or_tol,Tmax,rratio,rdpiv,fastsv, ...
                   store,doiter)

  % set default parameters
  if nargin <  1 || isempty(M), M = 16384; end
  if nargin <  2 || isempty(N), N =  8192; end
  if nargin <  3 || isempty(delta), delta = 1e-3; end
  if nargin <  4 || isempty(occ), occ = 256; end
  if nargin <  5 || isempty(p), p = 64; end
  if nargin <  6 || isempty(rank_or_tol), rank_or_tol = 1e-6; end
  if nargin <  7 || isempty(Tmax), Tmax = 2; end
  if nargin <  8 || isempty(rratio), rratio = 2; end
  if nargin <  9 || isempty(rdpiv), rdpiv = 'l'; end
  if nargin < 10 || isempty(fastsv), fastsv = 'n'; end
  if nargin < 11 || isempty(store), store = 'a'; end
  if nargin < 12 || isempty(doiter), doiter = 1; end

  % initialize
  m = ceil(sqrt(M)); [x1,x2] = ndgrid((1:m)/m); rx = [x1(:) x2(:)]';
  r = randperm(size(rx,2)); rx = rx(:,r(1:M));  % row points
  n = ceil(sqrt(N)); [x1,x2] = ndgrid((1:n)/n); cx = [x1(:) x2(:)]';
  r = randperm(size(cx,2)); cx = cx(:,r(1:N));  % col points
  clear x1 x2
  theta = (1:p)*2*pi/p; proxy = 1.5*[cos(theta); sin(theta)];  % proxy points
  % reference proxy points are for unit box [-1, 1]^2

  % factor matrix using RSKELFM
  Afun = @(i,j)Afun_(i,j,rx,cx,delta);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,delta);
  opts = struct('Tmax',Tmax,'rratio',rratio,'rdpiv',rdpiv,'fastsv',fastsv, ...
                'verb',1);
  tic; F = rskelfm(Afun,rx,cx,occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('rskelfm time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % compress matrix using IFMM
  rank_or_tol = max(rank_or_tol*1e-2,1e-15);  % higher accuracy for reference
  opts = struct('Tmax',Tmax,'store',store);
  tic; G = ifmm(Afun,rx,cx,occ,rank_or_tol,pxyfun,opts); t = toc;
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
function K = Kfun(x,y,delta)
  dx = x(1,:)' - y(1,:);
  dy = x(2,:)' - y(2,:);
  K = -1/(2*pi)*log(sqrt(dx.^2 + dy.^2 + delta^2));  % regularized
end

% matrix entries
function A = Afun_(i,j,rx,cx,delta)
  A = Kfun(rx(:,i),cx(:,j),delta);
end

% proxy function
function [Kpxy,nbr] = pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,delta)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  if rc == 'r'
    Kpxy = Kfun(rx(:,slf),pxy,delta);
    dr = cx(:,nbr) - ctr;
  else
    Kpxy = Kfun(pxy,cx(:,slf),delta);
    dr = rx(:,nbr) - ctr;
  end
  % proxy points form ellipse of scaled "radius" 1.5 around current box
  % keep among neighbors only those within ellipse
  nbr = nbr(sum((dr./l).^2) < 1.5^2);
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
