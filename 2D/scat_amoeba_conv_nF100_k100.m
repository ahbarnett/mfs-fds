% driver code for 2D Helmholtz exterior scattering from smooth big amoeba.
% Using MFS + FLAM FDS + FMM2D.   4/18/19

clear; v = 1;        % verbosity: 0 = text only, 1 = figs, 2 = movie only

% smooth amoeba radial shape params...
rng(0);     % freeze the randomness
nF = 100; aj = randn(1,nF); bj = randn(1,nF);  % # Fourier modes
%sc = 0.1/sqrt(nF); aj=sc*aj; bj=sc*bj;        % flat band
sc = min(0.1,0.5./(1:nF)); aj = sc.*aj;  bj = sc.*bj;
R = smoothfourierradfunc(aj,bj);

Ns = [5e3 1e4 1.5e4 2e4 2.5e4 3e4];   % convergence study

k = 100;    % wavenumber (fix all params if want to compare to ue value below)
ti = -pi/2; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x)));   % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

% linear solver choice...
meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel

flampar.rank_or_tol = 1e-12;
flampar.p=128; %64;       % num proxy pts (should depend on eps)
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=128;    % max pts per box for quadtree

% params for LSQ
lsqpar.meth = 'u';  % 'u'=underdetermined, 'o'=overdetermined
lsqpar.tau = eps^(-1/3);  % constraint weighting
lsqpar.lambda = 1e-6;  % regularization in OLS
lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
lsqpar.refine = 0;      % 0 or 1

us = nan*Ns;
for j = 1:numel(Ns);    % ---------------------------------- N-convergence
  N = Ns(j); fprintf('\nsolving N=%d...\n',N)

  % MFS: t = surface pt struct, s = source pt struct
  M = round(1.2*N);         % # bdry pts
  t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t);  % bdry pts
  s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);  % MFS src pts
  imagd = 0.0012; %0.15/sqrt(N);   % crucial src dist param:  fixed or scale w/ nF or N
  s.t = 1i*imagd + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
  if v==1 % & j==1
    figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
    axis equal; title(sprintf('N=%d, imagd=%g',N,imagd)); hold off; drawnow;
  end

  % convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
  rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

  F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
  if meth~='r', fprintf('eps-rank of A : %d\n',rank(F.A,1e-14*normest(F.A))), end
  rhs = -ui(t.x);     % eval uinc on bdry
  co = solve(F,rhs,meth);
  nrm = norm(co);
  ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
  rrms = norm(ur(:) - rhs)/sqrt(M);

  us(j) = mfseval(k,[real(x0);imag(x0)],cx,co,'d');  % 1 targ: always use direct

  fprintf('resid rms = %.3g \t co 2-nrm = %.3g \t Re(u) = %.15g\n',rrms,nrm,real(us(j)))

  m = round(N*sqrt(2));    % check at this many new shifted-grid bdry pts
  clear b; b.t = (0.5:m-0.5)'/m*2*pi; b.x = exp(1i*b.t).*R(b.t);
  ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
  utotb = ub + ui(b.x).';                 % physical potential on bdry
  fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m))  % u_tot on bdry = 0?
end                     % ------------------------------------

disp('Self-conv:  N       err at pt')
format short g; disp([Ns(:), abs(us(:)-us(end))]); format long g;


if v==1, dx = 0.01; xmax = 2.0; g = -xmax:dx:xmax;      % image plot grid
  [xx yy] = meshgrid(g,g);
  ii = ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu = mfseval(k,[xx;yy],cx,co,evalmeth);
  fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf; imagesc(g,g,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  colorbar; hold on; plot(t.x,'-');
end

if v==2, nam = 'amoeba_spin';  % movie
  wO=VideoWriter([nam '.avi']); wO.FrameRate=20; wO.open;
  dx = 0.01; gx = -3:dx:2; gy = -2:dx:2;    % image plot grid (wave from left)
  nx = numel(gx); if mod(nx,2)==1, nx=nx-1; gx = gx(1:end-1); end % for ffmpeg
  ny = numel(gy); if mod(ny,2)==1, ny=ny-1; gy = gy(1:end-1); end
  %figure(2); % testing
  cc = jet(256); cc(1,:) = [0 0 0];   % add a black one
  nti = 2000;
  for ti = (0:nti-1)/nti*2*pi
    [xx yy] = meshgrid(gx,gy);
    zz = exp(1i*ti)*(xx+1i*yy); xx = real(zz); yy = imag(zz); % rot the grid too
    ii = ~insideradfunc(R,xx+1i*yy);
    ug = nan*ii;   % the plot grid
    xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
    ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x)));   % u_inc
    rhs = -ui(t.x);     % eval uinc on bdry
    co = solve(F,rhs,meth);
    nrm = norm(co);
    ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
    rrms = norm(ur(:) - rhs)/sqrt(M);
    fprintf('th=%.6g: resid rms = %.3g \t co 2-nrm = %.3g\n',ti,rrms,nrm)
    tic; uu = mfseval(k,[xx;yy],cx,co,evalmeth); fprintf('grid eval %.3g s\n',toc);
    uu = uu + ui([xx + 1i*yy]);   % incident on outside grid
    ug(ii) = uu;     % overwrite ext only
    %imagesc(gx,gy,real(ug)); caxis(2*[-1 1]); axis equal tight xy; colorbar;
    %drawnow;  % testing
    z = min(255,max(2,round(129 + 128*(real(ug)/2.0))));
    z(~ii) = 1;    % black for nan
    rgb = zeros(ny,nx,3); for i=1:3, rgb(:,:,i)=reshape(cc(z,i),[ny,nx]); end
    writeVideo(wO,rgb);
  end
  close(wO);     % writes AVI movie out; now encode small MP4...
  movieavi2mp4(nam);
end
