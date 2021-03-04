% mfs_fds_simp.m
% Simplifed code which only contains the code that we need, original code,
% please see mfs_fds_hofixed_liu_3D.m

clear;
close all;
tic;
filename = 'FDS_3D_sphere.txt';
fid = fopen(filename, 'a');
global A r tau rx cx k proxy M N % make them global variable

% ============== FDS Parameter Setup ======================================
folder='FLAM-master';     % access FLAM
curpath = pwd;
%  curpath = '../code/src/FLAM';
dirs = {'core','geom','hifde','hifie','ifmm','mf','misc','quad','rskel', ...
    'rskelf'};
for s = dirs, addpath(sprintf('%s/%s/%s',curpath,folder, s{:})); end, clear s

% setup for Ho: (from ie_circle)
rank_or_tol = 1e-10;
p = 32; 
s = ellipsoid(1.5, 1.5, 1.5);
s = setupspherequad(s,[p, p/2]); 
proxy = s.x; 

opts = [];
opts.verb = 1;
occ = 128; % ??

tau = eps^(-1/3);  % params for LSQ
v = 0; % verbosity
%  checkmv = 0;
checkmv = 1;
fprintf(fid, '------------------------------------------------------------ \n');
% ====================== END ==============================================
ds = 0.2;
for idx = 1: length(ds)
    % ====================== Object Parameter setup ===========================
    N = 1e4; %400;     % # source points
    fprintf(fid, 'the parameters are tol = %e, p = %d, occ = %d, N = %d \n', rank_or_tol, p, occ, N);
    P = 100;
    N0 = N/P;
    M0 = 1.2 * N0;
    d = ds(idx);          % parameter controls distance between boundary points to source points;
    M = 1.2 * N;      % # boundary points
    k = 3;               % gives the wavenumber of the problem;
    theta = -pi/4;   phi = pi/3;
    kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
    
    % ==================== This is just a way to setup sphere as the boundary =
    a = 0;  w = 4; pp = 0;
    
    R = @(s) 1+a*cos(w*(s-pp));
    R_t = @(s) -a*w*sin(w*(s-pp));
    
    s = source3DQuad(R, R_t, N0, P, 1, d);
    t = curve3DQuad(R, R_t, M0, P, 1);
    
    f = QuasiPerVecFilling3D(t, kx, ky, kz);
    % =========================================================================
    
    
    rx = [t.x'; t.y'; t.z'];    % a sphere with radius 1 as boundary
    cx = [s.x'; s.y'; s.z'];    % a sphere with radius 1-0.4 as source points
    
    
    fprintf(fid, 'k = %d, a = %d, w = %d, d = %4.4f\n', k, a, w, d);
    % ======================== END ============================================
    fprintf(fid, 'setting up the parameter and preparation takes %4.2f seconds. \n', toc);
    fprintf('setting up the parameter and preparation takes %4.2f seconds. \n', toc);
    % ================= Call 1. rskel function ========================
    tic;
    F = rskel(@Afun,rx,cx,occ,rank_or_tol,@pxyfun,opts);
    fprintf(fid, 'doing rskel takes %4.2f seconds. \n', toc);
    fprintf('doing rskel takes %4.2f seconds. \n', toc);
    % =================================================================
    
    store = 'n'; opts = struct('store',store,'verb',0); %% changed to minimal storage
    
    % ======================= Call 2. rskel_xsp function ==============
    tic;
    A = rskel_xsp(F);  % make a big sparse mat
    fprintf(fid, 'doing rskel_xsp takes %4.2f seconds. \n', toc);
    fprintf('doing rskel_xsp takes %4.2f seconds. \n', toc);
    % =================================================================
    
    % ====================== Call 3. QR to factorize the matrix ======
    tic;
    A = [tau*A; speye(N) sparse(N,size(A,2)-N)];
    r = qr(A,0);
    fprintf(fid, 'doing QR takes %4.2f seconds. \n', toc);
    fprintf('doing QR takes %4.2f seconds. \n', toc);
    % =================================================================
    
    % ==================== Call 4. call ls to solve given the RHS ===== 
    tic;
    b = [tau*f; zeros(size(A, 1)-M,1)];
    x = lsfun(b);
    c = x(1:N);
    res1 = norm(A*x-b);
    res2 = norm(A(1:M, 1:N) * c - tau*f);
    toc;
    fprintf(fid, 'solve the right hand side (0 iteration) takes %4.2f seconds. \n', toc);
    fprintf('solve the right hand side (0 iteration) takes %4.2f seconds. \n', toc);
     
    % ======================= test the accuracy of the program =========
    tic;
    M_test = 41;   P_test = 41;
    t_test = curve3DQuad(R, R_t, M_test, P_test, 1);
    
    u_exact = QuasiPerVecFilling3D(t_test, kx, ky, kz);
    u_test = zeros(M_test * P_test, 1);
    for ii = 1:M_test*P_test
        dx = t_test.x(ii) - s.x;
        dy = t_test.y(ii) - s.y;
        dz = t_test.z(ii) - s.z;
        dr = sqrt(dx.^2 + dy.^2 + dz.^2);
        u_test(ii) = sum(c.* ( 1/(4*pi) * exp(1i * k * dr)./dr)); % eval u @ x0
    end
    err = norm(u_exact - u_test)/sqrt(length(u_exact))
    %err = norm(u_exact - u_test)/norm(u_exact)
    
    fprintf(fid, 'the error on the boundary is %e \n', err);
    fprintf(fid, 'doing the test takes %4.2f seconds.\n', toc);
    fprintf('doing the test takes %4.2f seconds.\n', toc);
    % =========================================================================
    
    a1 = whos('A');
    a2 = whos('F');
    a3 = whos('r');
    byte = a1.bytes + a2.bytes + a3.bytes + a3.bytes;
    fprintf(fid, 'it takes %4.2f GB of RAM.\n', byte/(1e9));
    fprintf(fid, '------------------------------------------------------------ \n\n');
end