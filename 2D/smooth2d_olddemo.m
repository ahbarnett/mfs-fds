function smooth2d_olddemo
% Test MFS with Ho's LSQ FDS and LU. Barnett 10/18/14
%addpath('~/physics/shravan/lsc2d/lsc2d_release1.01/'); % for quadr
%rmpath('/home/alex/svn/mpspack/');
%ahb 3/25/19

  curpath='../../FLAM';     % access FLAM
%  curpath = '../code/src/FLAM';
  dirs = {'compat','core','geom','hifde','hifie','ifmm','mf','misc','quad','rskel', ...
          'rskelf'};
  for s = dirs, addpath(sprintf('%s/%s',curpath,s{:})); end, clear s

% setup for Ho: (from ie_circle)
rank_or_tol = 1e-12;
p=64; theta = (1:p)*2*pi/p; proxy = 1.5*[cos(theta); sin(theta)]; clear theta
opts = []; opts.verb = 0; occ=128; % ??
%% lambda = 0;         changed to underdetermined setup--no more regularization
tau = eps^(-1/3);  % params for LSQ

v = 0; % verbosity
meth = 'l';   % lin alg solve: 'd'=dense, 'h'=Ho rskel LU, 'l'=Ho rskel LSQ
checkmv = 0;

b = 1; a = .3; w = 5; % try 25;        % smooth wobbly radial shape params
R = @(t) b*(1 + a*cos(w*t));
k = 30;            % wave #
ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

%d = 0.1; % makes A have bunch of sing vals at eps for N>400
%% Ns = 200:200:1000; %400;     % 400 is converged to 1e-14
Ns = 2.^(10:14); %% let's try something bigger
%Ns = 50000; % 12 secs for k=30, 20 sec for k=100
us = nan(numel(Ns),1); re = us; nrms = us; % data
for i=1:numel(Ns), N = Ns(i);
  M = N; %1.2*N;         %% also works for M > N
  t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t); t = quadr(t); % bdry pts
  f = -ui(t.x);  % RHS
  s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t); % MFS src pts
  
  % scale d w/ 1/N
  d = 100/N;
  
  s.t = 1i*d + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
  if v, figure(1); plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
    axis equal; title(sprintf('N=%d, d=%g',N,d)); hold off; drawnow; end
  if strcmp(meth,'d')  % mess with dense solves
    dz = repmat(t.x,[1 N]) - repmat(s.x.',[M 1]); % s = src(MFS), t = targ (bdy)
    tic, A = besselh(0,k*abs(dz)); fprintf('fill %.3g s\n',toc)     % matrix
    %c = A\f;     % solve, good
    %[q,r] = qr(A); c = r\(q'*f);
    tic, [L,U] = lu(A); fprintf('LU %.3g s\n',toc)
    tic, c = U\(L\f); fprintf('solve %.3g s\n',toc) % square, good even illcond
    % [L,U,P] = lu(A); c = U\(L\(P*f));  % also good for square
    %E = [eye(M), A; A', zeros(N)]; % aux system, bottoms at 11-12 digits
    %x = E \ [f;zeros(N,1)]; c = x(1+M:end); %[cond(A), cond(E)]
    re(i) = norm(A*c-f); % resid nrm

  else   % Kenneth Ho's FLAM
    % put rows (targs) rx and cols (src) cx in Ken notation...
    rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];
    %i = 7; j = 10; dz = abs(t.x(i)-s.x(j)); besselh(0,k*abs(dz)), Afun(i,j) % test Afun works
    t0=tic; F = rskel(@Afun,rx,cx,occ,rank_or_tol,@pxyfun,opts);
    w = whos('F'); fprintf('rskel: %.3g s,  %6.2f (MB)\n',toc(t0),w.bytes/1e6)
    % compress matrix using IFMM  ... needed for making rhs, and ls(), why?
    store = 'n'; opts = struct('store',store,'verb',0); %% changed to minimal storage
    t0=tic; G = ifmm(@Afun,rx,cx,2*occ,1e-15,@pxyfun,opts); % if want resid
    fprintf('ifmm: %.3g s\n',toc(t0))
    if checkmv & N<2e3 % check if rskel and ifmm apply vs dense A:
      dz = repmat(t.x,[1 N]) - repmat(s.x.',[M 1]); % s = src(MFS), t=targ (bdy)
      tic, Ad = besselh(0,k*abs(dz)); fprintf('dense fill %.3g s\n',toc)
      x = randn(N,1); y = Ad*x;
      fprintf('rskel apply err = %.3g\n',norm(y - rskel_mv(F,x)))
      fprintf('ifmm apply err = %.3g\n',norm(y - ifmm_mv(G,x,@Afun)))
    end
    t0=tic; A = rskel_xsp(F);  % make a big sparse mat
    t1 = toc(t0); w = whos('A');
    fprintf('xsp: %.3g sec / %6.2f (MB)\n',t1,w.bytes/1e6);
    if strcmp(meth,'h')  % try square lu sparse solve from ie_square...
      t0=tic; [L,U] = lu(A); t1=toc(t0); % square.  A is sparse, uses UMFPACK
      w = whos('L'); spmem = w.bytes;w = whos('U');spmem=(spmem + w.bytes)/1e6;
      fprintf('lu: %.3g sec / %6.2f (MB)\n',t1,spmem)
      t0=tic; c = sv(f); fprintf('sv: %.3g s\n',toc(t0))  % solve stage
    elseif strcmp(meth,'l')  % lsq from ols_square...
%%      A = [tau*A(M+1:end,:); A(1:M,:); lambda*speye(N) sparse(N,size(A,2)-N)];
%% change to underdetermined formulation
      A = [tau*A; speye(N) sparse(N,size(A,2)-N)];
      tic, r = qr(A,0); t1 = toc; w = whos('r');
      fprintf('qr: %.3g sec / %6.2f (MB)\n',t1,w.bytes/1e6)
      % set up right-hand side
      Ms = size(A,1);
%%      nC = Ms - M - N;  % cN needed by ls
      nC = Ms - N;
      X = f;  % my RHS
      t0=tic; [c,cres,niter] = ls([X; zeros(nC-M,1)]); %%[c,cres,niter] = ls(X); % is 1st arg the output vector?
      fprintf('equ-const lsq %.3g s, cres=%.3g, niter=%d\n',toc(t0),cres,niter)
    end
    Ac = ifmm_mv(G,c,@Afun); % that's why needed ifmm setup
    re(i) = norm(Ac-f); % resid nrm
    if v, figure(3); plot([real(f),real(Ac)],'+-'); title('f,Ac'); end
  end
  nrms(i) = norm(c);
  dz = x0 - s.x; us(i) = sum(c.*besselh(0,k*abs(dz))); % eval u @ x0
end
us % vals
ue =  -0.764438165507287 +     0.357188677072235i;  % converged value k=30
disp('N, resid nrm,  u err at pt,  soln 2-norm:')
format short g; disp([Ns(:), re, abs(us-ue), nrms]); format long g
if v&numel(Ns)>1, figure(2); semilogy(Ns,abs(us-ue),'+-'); xlabel('N');
hold on; semilogy(Ns,nrms,'g+-'); legend('u pt convergence','||c||_2'); end

%keyboard % for debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ho's codes need to access local vars like global ones - gets confusing:

  % kernel function       - from 2d laplace in ie_circle
  function K = Kfun(x,y,lp)
    dx = bsxfun(@minus,x(1,:)',y(1,:));
    dy = bsxfun(@minus,x(2,:)',y(2,:));
    dr = sqrt(dx.^2 + dy.^2);

%% no need for dipoles
%%    if strcmpi(lp,'s')
      K = besselh(0,k*dr);   % helm
      %K = -1/(2*pi)*log(dr);  % laplace
%%     elseif strcmpi(lp,'d')
%%      rdotn = bsxfun(@times,dx,y(1,:)) + bsxfun(@times,dy,y(2,:));
%%      K = k*besselh(1,k*dr).*rdotn./dr;  % helm
      %K = 1/(2*pi).*rdotn./dr.^2;  % laplace
%%     end
  end

  % matrix entries        - monopoles, edited from ie_circle
  function A = Afun(i,j)
%%    A = Kfun(rx(:,i),cx(:,j),'s');  % nonsymm taken from ols_square
    A = Kfun(rx(:,i),cx(:,j));
  end

  % proxy function - from ie_circle, with weights removed, why 2pi/N here?
  function [Kpxy,nbr] = pxyfun(rc,rx,cx,slf,nbr,l,ctr)
    pxy = bsxfun(@plus,proxy*l,ctr');
    if strcmpi(rc,'r')
%%      Kpxy = Kfun(rx(:,slf),pxy,'s')*(2*pi/N); % monopoles for proxy
%% note: (2*pi/N) is a quadrature weight for unit circle geometry
      Kpxy = Kfun(rx(:,slf),pxy)*(2*pi/N);
      dx = cx(1,nbr) - ctr(1);
      dy = cx(2,nbr) - ctr(2);
    elseif strcmpi(rc,'c')
%%      Kpxy = Kfun(pxy,cx(:,slf),'d')*(2*pi/N); % why flips to dipoles here?
      Kpxy = Kfun(pxy,cx(:,slf))*(2*pi/N);
      dx = rx(1,nbr) - ctr(1);            % (makes no difference)
      dy = rx(2,nbr) - ctr(2);
    end
    dist = sqrt(dx.^2 + dy.^2);
    nbr = nbr(dist/l < 1.5);
  end

  % sparse LU solve - from ie_circle, not used for ols_square
  function Y = sv(X)
    X = [X; zeros(size(A,1)-N,size(X,2))];   % pad rhs to extended sys
    if strcmpi(F.symm,'n')
      Y = U\(L\X);               % non symm back-subst to solve
    else
      Y = P*(L'\(D\(L\(P'*X))));
    end
    Y = Y(1:N,:);
  end

  % least squares solve - from ols_square. Needs R (I changed to r), A
  function x = lsfun(b)
    x = r\(r'\(A'*b));
    x = x + r\(r'\(A'*(b - A*x)));
  end

  % equality-constrained least squares solve - from ols_square. Needs nC,N,tau
%    function [Y,cres,niter] = ls(X)
%      n = size(X,2);
%      X = [X; zeros(N,n)];
%      [Y,cres,niter] = lsedc(@lsfun,A(nC+1:end,:),X,A(1:nC,:),zeros(nC,n),tau);
%      Y = Y(1:N,:);
%    end

%% underdetermined LSE solve (from uls_circle)
  function [Y,cres,niter] = ls(X)
    n = size(X,2);
    [Y,cres,niter] = lsedc(@lsfun,A(nC+1:end,:),zeros(N,n),A(1:nC,:),X,tau);
    Y = Y(1:N,:);
  end


end
%%%%%%%%%%%%%%%% now Alex funcs living outside main: %%%%%%%%

function s = quadr(s, N)
% QUADR - set up periodic trapezoid quadrature & geom for smooth closed curve
%
% s = quadr(s,N) where s is a struct containing a parametrization of the curve
%  in the form of s.Z a function from [0,2pi) to the complex plane, uses this
%  to build the set of nodes, weights, speeds, curvatures, etc, using the
%  N-node PTR on [0,2pi).
%
% s = quadr(s) where s contains at least the field s.x (column vector of N node
%  locations) generates all other fields by spectral differentiation.
%
% If s.Zp is present it is used as the function Z'; likewise s.Zpp is used for
%  Z''. The point of using these is that slightly more accurate normals and
%  curvatures may result compared to differentiation. On the other hand, Zp and
%  Zpp are not checked for plausibility.
%
% Note: curves go counter-clockwise. Node j is at Z(2pi.j/N). FFT is used so
%  it's O(N log N), or O(N) if Z, Zp and Zpp available.
%
% quadr('test') runs tests.
%
% Inputs:
%  s : struct containing either the field s.Z, a mapping from [0,2pi) to
%      the complex plane, which must be able to accept vectorized inputs,
%      or the field s.x containing column vector of the nodes.
%  N : number of nodes (not needed if s.x is present; however, if s.Z is also
%      present, N overrides the number of nodes in s.x).
% Outputs:
%  s : same struct with added column vector fields (nodes, weights, velocities,
%      curvatures, etc) needed for quadrature and Nystrom methods. Namely,
%      s.x nodes in complex plane, Z(s_j)
%      s.xp velocities Z'(s_j)
%      s.xp accelerations Z''(s_j)
%      s.t parameter values s_j for trapezoid rule, 2pi.j/N, j=1,...,N
%      s.nx outward unit normals
%      s.tang forward unit tangential vectors
%      s.sp speeds |Z'(s_j)|
%      s.w "speed weights" (2pi/N)*s.sp
%      s.cw velocity weights (ie complex speed)
%      s.cur curvatures kappa(s_j) (inverse bending radius)
%
% Example usage:
%
% s.Z = @(s) (1+0.3*cos(5*s).*exp(1i*s);                   % starfish param
% s = quadr(s,100);
% figure; plot(s.x,'k.'); hold on; plot([s.x, s.x+0.2*s.nx].', 'b-'); axis equal
%
% Now check that normals from spectral differentiation are accurate:
%
% s.Zp = @(s) -1.5*sin(5*s).*exp(1i*s) + 1i*s.Z(s);        % Z' formula
% t = quadr(s,100); norm(t.nx-s.nx)   % should be small
%
% (c) Alex Barnett 10/8/14

if strcmp(s,'test'), testquadr; return; end
if nargin>1           % use N from args
  s.t = (1:N)'*(2*pi/N);
  if isfield(s,'Z'), s.x = s.Z(s.t); end    % use formula
  if N~=length(s.x), error('N differs from length of s.x; that sucks!'); end
elseif isfield(s,'x')
  s.x = s.x(:);          % ensure col vec
  N = length(s.x);
  s.t = (1:N)'*(2*pi/N); % we don't know the actual params, but choose this
else
  error('Need to provide at least s.Z and N, or s.x. Neither found!');
end
if isfield(s,'Zp'), s.xp = s.Zp(s.t); else, s.xp = perispecdiff(s.x); end
if isfield(s,'Zpp'), s.xpp = s.Zpp(s.t); else, s.xpp = perispecdiff(s.xp); end
% Now local stuff that derives from x, xp, xpp at each node...
s.sp = abs(s.xp);
s.tang = s.xp./s.sp;
s.nx = -1i*s.tang;
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2;  % recall real(conj(a)*b) = "a dot b"
s.w = (2*pi/N)*s.sp;
s.cw = (2*pi/N)*s.xp;  % complex weights (incl complex speed)
end

function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = PERISPECDIFF(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points).
%
% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));
end

function testquadr
% not very extensive! Barnett 10/8/14
Z = @(s) (1+0.3*cos(5*s)).*exp(1i*s);                   % starfish param
s.Z = Z;
s = quadr(s,100);
s = []; s.x = Z((1:100)/100*2*pi); s = quadr(s);  % testing s.x input only
figure; plot(s.x,'k.'); hold on; plot([s.x, s.x+0.2*s.nx].', 'b-'); axis equal
% Now check that normals from spectral differentiation are accurate:
Zp = @(s) -1.5*sin(5*s).*exp(1i*s) + 1i*Z(s);        % Z' formula
s.Zp = Zp;
t = quadr(s,100); norm(t.nx-s.nx)   % should be small
s = []; s.x = 3; s.Z = Z; s = quadr(s,100) % N should override s.x
%s = []; s.x = Z((1:50)'/50*2*pi); s=quadr(s,100); should fail
end
