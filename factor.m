function F = factor(k,rx,cx,meth,flampar,lsqpar)
% FACTOR  fast factor LSQ system, for 2D or 3D PDE kernel, MFS plus rskel FDS.
%
% F = factor(k,rx,cx,meth,flampar,lsqpar)
% 
% Inputs:
%  k = wavenumber >= 0
%  rx = d-by-M surface pts,  where dim d = 2 or 3 (determined from coords)
%  cx = d-by-N source (MFS) pts
%  meth = 'l' dense direct LU
%         'q' dense direct QR
%         'r' FLAM rskel, augmented by I.
%  flampar = FLAM parameter struct
%  lsqpar = LSQ parameter struct
%
% Output:
%  F = struct containing all factor info.

dim = size(rx,1);
if dim<2 || dim>3, error(sprintf('dim = %d should be 2 or 3!',dim)); end

if meth=='l' || meth=='q'  % ---------------- dense
  tic;
  if dim==2
    r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2);
    A = besselh(0,k*r);
  else
    r = sqrt((rx(1,:)'-cx(1,:)).^2 + (rx(2,:)'-cx(2,:)).^2 + (rx(3,:)'-cx(3,:)).^2);
    A = exp(1i*k*r)./r;
  end
  w = whos('A'); fprintf('A fill %.3g s \t (%.0f MB)\n',toc,w.bytes/1e6)
  tic;
  if meth=='l'   % M=N only
    [F.L,F.U] = lu(A);
    fprintf('dense LU %.3g s\n',toc)
  else
    [F.Q,F.R] = qr(A,0);      % econ
    fprintf('dense QR %.3g s\n',toc)
  end
  F.A = A;   % save it for kicks

else                       % ---------------- FDS
  % ***
  
  
end

