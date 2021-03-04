% mfs_fds_simp.m
% Simplifed code which only contains the code that we need, original code,
% please see mfs_fds_hofixed_liu_3D.m

clear; 
close all;
tic;
filename = 'FDS_3D_ellipsoid.txt';
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

ds = 0.1;
for jj = 1: length(ds)
    d = ds(jj);
    fprintf(fid, '------------------------------------------------------------ \n');
    % ====================== END ==============================================
    
    % ====================== Object Parameter setup ===========================
    N = 1e4; %400;     % # source points
    fprintf(fid, 'the parameters are tol = %e, p = %d, occ = %d, N = %d \n', rank_or_tol, p, occ, N);
    P = 100; 
    N0 = N/P;  
    M = 1.2 * N;      % # boundary points
    M0 = M/P;
    k = 3;               % gives the wavenumber of the problem;
    theta = -pi/4;   phi = pi/3;
    kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
    
    % ==================== This is just a way to setup sphere as the boundary =
    Rx = 0.6;
    Ry = 0.6;
    Rz = 1.0;
    
    t = getShapeEllip(Rx, Ry, Rz, M0, P);
    s = getSourceEllip(Rx, Ry, Rz, d, N0, P);
    f = QuasiPerVecFilling3D(t, kx, ky, kz);
    
    rx = [t.x'; t.y'; t.z'];    % a sphere with radius 1 as boundary
    cx = [s.x'; s.y'; s.z'];    % a sphere with radius 1-0.4 as source points
    
    
    fprintf(fid, 'k = %d, d = %4.4f, Rx = %4.4f, Ry = %4.4f, Rz = %4.4f.\n', k, d, Rx, Ry, Rz);
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
    if 1
        M_test = 40;
        
        t_test = getShapeEllip(Rx, Ry, Rz, M_test, M_test);
        
        u_exact = QuasiPerVecFilling3D(t_test, kx, ky, kz);
        u_test = zeros(M_test^2, 1);
        for ii = 1:M_test^2
            dx = t_test.x(ii) - s.x;
            dy = t_test.y(ii) - s.y;
            dz = t_test.z(ii) - s.z;
            dr = sqrt(dx.^2 + dy.^2 + dz.^2);
            u_test(ii) = sum(c.* ( 1/(4*pi) * exp(1i * k * dr)./dr)); % eval u @ x0
        end
        
        err = norm(u_exact - u_test)/sqrt(length(u_exact))
        %err = norm(u_exact - u_test)/norm(u_exact);
    end
    fprintf(fid, 'the error on the boundary is %e \n', err);
    fprintf(fid, 'doing the test takes %4.2f seconds.\n', toc);
    fprintf('doing the test takes %4.2f seconds.\n', toc);
    % =========================================================================
    
    a1 = whos('A');
    a2 = whos('F');
    a3 = whos('r');
    byte = a1.bytes + a2.bytes + a3.bytes + a3.bytes;
    fprintf(fid, 'it takes %4.2f GB of RAM.\n', byte/(1e9));
    fprintf('it takes %4.2f GB of RAM.\n', byte/(1e9));
    fprintf(fid, '------------------------------------------------------------ \n\n');
    
end