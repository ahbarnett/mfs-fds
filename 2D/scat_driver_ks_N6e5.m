% driver code for 2D Helmholtz exterior scattering from smooth domain.
% using MFS and FLAM FDS. Small example to test: FLAM and FMM2D.

% Barnett 3/25/19. 4/8/19: added FMM2D

clear;

b = 1; a = .3; w = 5;    % try 25;        % smooth wobbly radial shape params
R = @(t) b*(1 + a*cos(w*t));

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

N = 600000; 
ks = 30:30:3000; 
err = zeros(length(ks), 1); 
rrms = zeros(length(ks), 1); 
ds = zeros(length(ks), 1); 
for i = 1:length(ks)
    fprintf('i = %d', i); 
    k = ks(i); 
    ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
    % MFS: t = surface pt struct, s = source pt struct
    M = round(1.2*N);         % # bdry pts
    t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t);  % bdry pts
    s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);  % MFS src pts
    imagd = min(0.1, min(abs(diff(s.x)))*250);                   % Note for imagd=0.1 in this BVP, rank(A)<600
    ds(i) = imagd; 
    %imagd = 0.1;
    
    s.t = 1i*imagd + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
    % convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
    rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];
    
    % linear solver choice...
    meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel
    
    flampar.rank_or_tol = 1e-12;
    flampar.p = 2048;  % num proxy pts (should depend on eps)
    flampar.opts = []; flampar.opts.verb = 0;
    flampar.occ = 128;    % max pts per box for quadtree
    
    % params for LSQ
    lsqpar.meth = 'u';  % 'u'=underdetermined, 'o'=overdetermined
    lsqpar.tau = eps^(-1/3);  % constraint weighting
    lsqpar.lambda = 1e-6;  % regularization in OLS
    lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
    lsqpar.refine = 0;      % 0 or 1
    tic; 
    F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
    rhs = -ui(t.x);     % eval uinc on bdry
    co = solve(F,rhs,meth);
    nrm = norm(co);
    ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
    rrms(i) = norm(ur(:) - rhs)/sqrt(M);
    
    m = 137;    % check at this many new shifted-grid bdry pts
    clear b; b.t = (0.5:m-0.5)'/m*2*pi; b.x = exp(1i*b.t).*R(b.t);
    ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
    utotb = ub + ui(b.x).';                 % physical potential on bdry
    err(i) = norm(utotb)/sqrt(m); 
    toc;
end
save('scat_N6e5_ks.mat', 'err', 'rrms', 'ks', 'ds');