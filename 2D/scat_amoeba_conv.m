% driver code for 2D Helmholtz exterior scattering from smooth amoeba.
% Using MFS + FLAM FDS + FMM2D.   4/19/19

% TODO: Currently limited to 7 digits @ k=100 - why? proxy = 96 helps a bit.
%      (not refine, nor imagd, nor N, nor flampar.p nor flampar.rank_or_tol)
%      at =30,  ~8 digits.
%      at k=1 or 10, limit at 8-9 digits
%      at k=300 need flampar.p >= 200 for 1e-6 and >=300 for 1e-7, then stops.
%   To check this loss is only due to FLAM, select meth='q' below & you'll
%      see 11 digits rms at N=6e3 (& solution at pt conv to 13 digits).

clear; v = 1;        % verbosity: 0 = text only, 1 = figs, 2 = movie only

% smooth amoeba radial shape params...
rng(0);     % freeze the randomness
expt = 1;
switch expt
 case 1    %  small-scale
  k = 30; flampar.p = 96; % wavenumber; 64 num proxy pts (should dep on eps, k)
  %k = 100; flampar.p = 100;   % wavenumber
  %k = 300; flampar.p = 300; Ns = [1e4 1.5e4 2e4];   % wavenumber
  nF = 20;        % # shape Fourier modes
  amax = 0.05; atail = 0.5;    % shape mode sizes
  imagd = 0.004;     % MFS src dist, good for nF=20
  Ns = [2e3 3e3 4e3 6e3 8e3 1e4];   % convergence- bottoms at 1e-7
 case 2         % medium big (30 s solve)
  k = 100;    % wavenumber
  nF = 200;
  amax = 0.05; atail = 0.5;
  imagd = 0.0005;    % crucial src dist param:  fixed or scale w/ nF or N 
  Ns = [3e4 4e4 5e4];   % convergence study
end
aj = randn(1,nF); bj = randn(1,nF);  % # Fourier modes
ampl = min(amax,atail./(1:nF)); aj = ampl.*aj;  bj = ampl.*bj;  % mode sizes
R = smoothfourierradfunc(aj,bj);

ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x)));   % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

% linear solver choice... (try 'l' to compare FLAM loss of digits)
meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel

flampar.rank_or_tol = 1e-12;    % 1e-15 is 5x slower than 1e-12
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=300; %128;    % max pts per box for quadtree (affects speed only)

lsqpar.tau = eps^(-1/3);  % params for LSQ
lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
lsqpar.refine = 0;      % 0 or 1

us = nan*Ns;
for j = 1:numel(Ns);    % ---------------------------------- N-convergence
  N = Ns(j); fprintf('\nsolving N=%d...\n',N)

  % MFS: t = surface pt struct, s = source pt struct
  M = round(1.2*N);         % # bdry pts
  t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t);  % bdry pts
  s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);  % MFS src pts
  s.t = 1i*imagd + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
  if v==1
    figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
    axis equal; title(sprintf('N=%d, imagd=%g',N,imagd)); hold off; drawnow;
  end

  % convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
  rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];
  
  F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
  if meth~='r' && N<2e3, fprintf('eps-rank of A : %d\n',rank(F.A,1e-14*normest(F.A))), end
  rhs = -ui(t.x);     % eval uinc on bdry
  co = solve(F,rhs,meth);
  nrm = norm(co);
  ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
  rrms = norm(ur(:) - rhs)/sqrt(M);     % resid RMS
  %(sim to rel resid, since rhs has elements of size 1)

  us(j) = mfseval(k,[real(x0);imag(x0)],cx,co,'d');  % 1 targ: always use direct

  fprintf('resid rms = %.3g \t co 2-nrm = %.3g \t Re(u) = %.15g\n',rrms,nrm,real(us(j)))

  m = round(N*sqrt(2));    % indep check at this many new shifted-grid bdry pts
  clear b; b.t = (0.5:m-0.5)'/m*2*pi; b.x = exp(1i*b.t).*R(b.t);
  ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth,struct('iprec',5));
  utotb = ub + ui(b.x).';                 % physical potential on bdry
  fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m))  % u_tot on bdry = 0?
end                     % ------------------------------------

disp('Self-conv:  N       err at pt')
format short g; disp([Ns(:), abs(us(:)-us(end))]); format long g;

if v==1, dx = 0.01; xmax = 2.0; g = -xmax:dx:xmax;      % image plot grid
  [xx yy] = meshgrid(g,g);
  ii = ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu = mfseval(k,[xx;yy],cx,co,evalmeth,struct('iprec',2));
  fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf; imagesc(g,g,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  colorbar; hold on; plot(t.x,'-');
end
