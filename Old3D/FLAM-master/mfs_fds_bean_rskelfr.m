% mfs_fds_simp.m
% Simplifed code which only contains the code that we need, original code,
% please see mfs_fds_hofixed_liu_3D.m

clear all;
close all;
tic;
filename = 'FDS_3D_results_Bean.txt';
fid = fopen(filename, 'a');
global  A r nC tau rx cx k proxy M N % make them global variable

% ============== FDS Parameter Setup ======================================
folder='FLAM-master';     % access FLAM
curpath = pwd;
%  curpath = '../code/src/FLAM';
dirs = {'core','geom','hifde','hifie','ifmm','mf','misc','quad','rskel', ...
    'rskelfr'};
for s = dirs, addpath(sprintf('%s/%s/%s',curpath,folder, s{:})); end, clear s

% setup for Ho: (from ie_circle)
rank_or_tol = 1e-11;
p = 512;
proxy = randn(3,p);
proxy = 1.5*bsxfun(@rdivide,proxy,sqrt(sum(proxy.^2)));

opts = [];
opts.verb = 1;
occ=2048; % ??

tau = eps^(-1/3);  % params for LSQ
v = 0; % verbosity
%  checkmv = 0;
checkmv = 1;
fprintf(fid, '------------------------------------------------------------ \n');
% ====================== END ==============================================
Ns = 220.^2;
ds = 0.15;
for i = 1:length(Ns)
    for j = 1:length(ds)
        N = Ns(i);
        % ====================== Object Parameter setup ===========================
        fprintf(fid, 'the parameters are tol = %e, p = %d, occ = %d, N = %d \n', rank_or_tol, p, occ, N);
        
        M = (1.2)^2 * N;      % # boundary points
        
        
        
        k = 10;               % gives the wavenumber of the problem;
        theta = -pi/4;   phi = pi/3;
        %kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta);
        kx = k;
        ky = 0;
        kz = 0;
        % ==================== This is just a way to setup sphere as the boundary =
        Rx = 0.6;
        Ry = 0.6;
        Rz = 1.0;
        alpha = 0.5;
        d = ds(j);  % 0.08 is the best for the bean shape(0.6, 0.6, 1.0, 0.5);
        
        
        t = getShapeBean(Rx, Ry, Rz, alpha, sqrt(M), sqrt(M));
        s = getSourceBean(Rx, Ry, Rz, alpha, d, sqrt(N), sqrt(N));
        f = QuasiPerVecFilling3D(t, kx, ky, kz);
        % =========================================================================
        if 0,
            %plot3(t.x, t.y, t.z, '.', s.x, s.y, s.z, '*r'); axis equal;
            X = reshape(t.x, sqrt(M), sqrt(M));
            Y = reshape(t.y, sqrt(M), sqrt(M));
            Z = reshape(t.z, sqrt(M), sqrt(M));
            F = reshape(real(f), sqrt(M), sqrt(M));
            surf(X,Y,Z, F);
            shading interp; xlabel('x'); ylabel('y'); zlabel('z');
            colormap hsv
            set(gca, 'FontSize',16);
            axis equal;
            %view(-50,30)
            camlight left;
            colorbar;
            saveas(gcf, '/Users/Yuxiang/Documents/AAResearch/dcthesis/bean.eps', 'epsc2');
            
        else
            
            rx = [t.x'; t.y'; t.z'];    % a sphere with radius 1 as boundary
            cx = [s.x'; s.y'; s.z'];    % a sphere with radius 1-0.4 as source points
            
            
            fprintf(fid, 'k = %d, d = %4.4f, Rx = %4.4f, Ry = %4.4f, Rz = %4.4f, twistness alpha = %4.4f .\n', k, d, Rx, Ry, Rz, alpha);
            % ======================== END ============================================
            fprintf(fid, 'setting up the parameter and preparation takes %4.2f seconds. \n', toc);
            fprintf('setting up the parameter and preparation takes %4.2f seconds. \n', toc);
            % ================= Call 1. rskelfr function ========================
            tic;
            F = rskelfr(@Afun,rx,cx,occ,rank_or_tol,@pxyfun,opts);
            fprintf(fid, 'doing rskelfr takes %4.2f seconds. \n', toc);
            fprintf('doing rskelfr takes %4.2f seconds. \n', toc);
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
            
            if 1,
                tic;
                b = [tau*f; zeros(size(A, 1)-M,1)];
                x = lsfun(b);
                c = x(1:N);
                res1 = norm(A*x-b);
                res2 = norm(A(1:M, 1:N) * c - tau*f);
                toc;
                fprintf(fid, 'solve the right hand side (0 iteration) takes %4.2f seconds. \n', toc);
                fprintf('solve the right hand side (0 iteration) takes %4.2f seconds. \n', toc);
                
            else
                tic;
                Ms = size(A,1);
                %      nC = Ms - M - N;  % cN needed by ls
                nC = Ms - N;
                X = f;  % my RHS
                [c,cres,niter] = ls([X; zeros(nC-M,1)]); %%[c,cres,niter] = ls(X); % is 1st arg the output vector?
                fprintf(fid, 'solve the right hand side (2 iterations) takes %4.2f seconds. \n', toc);
                fprintf('solve the right hand side (2 iterations) takes %4.2f seconds. \n', toc);
                cres
                niter
                % ==================================================================
            end
            
            
            % ======================= test the accuracy of the program =========
            tic;
            if 1
                M_test = 40;
                
                t_test = getShapeBean(Rx, Ry, Rz, alpha, M_test, M_test);
                
                u_exact = QuasiPerVecFilling3D(t_test, kx, ky, kz);
                u_test = zeros(M_test^2, 1);
                for ii = 1:M_test^2
                    dx = t_test.x(ii) - s.x;
                    dy = t_test.y(ii) - s.y;
                    dz = t_test.z(ii) - s.z;
                    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
                    u_test(ii) = sum(c.* ( 1/(4*pi) * exp(1i * k * dr)./dr)); % eval u @ x0
                end
                
                err(i,j) = norm(u_exact - u_test)/sqrt(length(u_exact));
                %err = norm(u_exact - u_test)/norm(u_exact);
            end
            fprintf(fid, 'the error on the boundary is %e \n', err(i));
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
            cc{i,j} = c;
            %save('err_40_Bean.mat', 'A','err', 'cc', 'ds', 'Ns');
        end
    end
end
