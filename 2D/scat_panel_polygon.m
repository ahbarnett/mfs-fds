% driver for Helmholtz exterior scattering from polygon.
% Barnett 5/6/19
% T

clear; v=1;   % verbosity
k = 300;       % wavenumber (for skyline keep <100 for nonres)
paneps = 1e-10;  % desired acc, really controls pmfs & p in panels
shape = 'u';     % 'r' rightangles city skyline, 'u' UK map

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.
% linear solver choice...
meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel
flampar.rank_or_tol = 1e-10;
flampar.p=200; %96; %64;  % 64 for k<50, 100 for k=100, 200 for k=300.
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=512;    % max pts per box for quadtree (might need to tune)
% note: v sens, even eps^(-1/3) can cause umfpack singular -> nan.
lsqpar.meth = 'u';  % 'u'=underdetermined, 'o'=overdetermined
lsqpar.tau = 1e5; %eps^(-1/3);  % constraint weighting  (>=1e7 causes NaNs as N gro)
lsqpar.lambda = 1e-6;  % regularization in OLS
lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
lsqpar.refine = 0;      % 0 or 1 (latter better acc)

if shape=='r' % make rand poly w/ right angles only, in [-1,1]^2
  rng(0)
  nb = 100;  % # bars
  ysc = 10/nb;
  h = 2/nb;       % horiz width of bars
  nc = (nb+1)*2;    % # corners
  c(1:2) = [-1-.2i*ysc, 1-.2i*ysc];    % bottom corners
  for b=1:nb, c(b*2+(1:2)) = 1i*rand*ysc+1-h*[b-1,b]; end
  maxL = 1.5*h;                      % max physical pan size, handles close
elseif shape=='u'    % UK map
  c = load('shapes/UK/uk_1poly_88'); c=c.zuk;
  maxL = 0.1; h = 0.1;       % only for normals plot
end
if v>1, figure(4); clf; plot([c c(1)], '.-'); title('polygon'); end

maxL = min(maxL, 1.5*pi/k);        % < 0.75 lambda scale
o=[]; o.geomeps = 1e-15;           % min pan size: ~eps and get sing->nan :(
o.qtype = 'u';                     % surf panel quadr type 'u' or 'g'

ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x)));  % u_inc
x0 = 0.5 + 1.3i;   % scatt wave soln test pt (complex notation)
uex = nan;

[t s] = build_polygon(c,'e',maxL,paneps,o);
N = numel(s.x);   % # dofs (MFS pts)
M = numel(t.x);   % # eqns (colloc pts)

if v, figure(1); clf; plot(t.x,'.'); hold on; plot(s.x, '.'); axis equal
  plot([c;c+0.2*h*t.cornors],'r-');  % show mean normals at corners
  m = (c+circshift(c,-1))/2;      % midpoints of arcs for labeling
  if v>1, for i=1:nc   % label corner MFS params
    text(real(c(i)),imag(c(i)),sprintf('\\alpha=%.3g r=%.2g Nr=(%d,%d) d=%.2g',t.angs(i), t.reffac(i), t.Nr(i,1), t.Nr(i,2), t.reldist(i)))
    text(real(m(i)),imag(m(i)),sprintf('Nb=%d',t.Nb(i)));
  end, end
  plot(x0,'+'); drawnow
end

fprintf('N=%d...  (min src ppw = %.3g; bdry ppw = %.3g)\n',N,2*pi/(k*max(abs(diff(s.x)))),2*pi/(k*max(abs(diff(t.x)))))

% convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
rhs = -ui(t.x);     % eval uinc on bdry
co = solve(F,rhs,meth);
i=1;  % dummy
nrm(i) = norm(co);
ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
rrms(i) = norm(ur(:) - rhs)/sqrt(M);

us(i) = mfseval(k,[real(x0);imag(x0)],cx,co,'d');  % 1 targ: always use direct

fprintf('pmfs=%d p=%d Np=%d: residrms=%.3g, co 2-nrm=%.3g, Re(u)=%.15g\n',t.pmfs,t.p,N/t.pmfs,rrms(i),nrm(i),real(us(i)))

o.p = 7;   % check at new bdry pts, this many (neq p) per panel
b = build_polygon(c,'e',maxL,paneps,o); m = numel(b.x);
ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
utotb = ub + ui(b.x).';                 % physical potential on bdry
brms(i) = norm(utotb)/sqrt(m);
fprintf('\tbdry rms err: %.3g \t est u err at pt = %.3g\n',brms(i),abs(us(i)-uex)) % u_tot on bdry = 0?

if v, figure(3); clf; subplot(2,1,1); semilogy(1:m, abs(utotb), '.-'); % err
  xlabel('test pt ind'); ylabel('total pot on bdry'); axis tight;
  subplot(2,1,2); semilogy(1:N, abs(co), '.-'); xlabel('src pt ind j');
  ylabel('coeff_j'); axis tight; drawnow
end

if v, dx = 0.005; gx = -2:dx:2; gy = min(imag(c))-1:dx:max(imag(c))+1; % plot
  if shape=='u', dx = 0.002; gx = -1:dx:1; gy = -1.5:dx:1.5; end  % for uk
  [xx yy] = meshgrid(gx,gy);
  ii = ~t.inside(xx+1i*yy);
  %ii = true(size(xx));
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu =  mfseval(k,[xx;yy],cx,co,evalmeth);
  fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf;
  imagesc(gx,gy,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  colorbar; hold on; plot(t.x,'-'); %plot(x0,'k+'); drawnow
end
