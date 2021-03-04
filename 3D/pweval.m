function [u, un] = pweval(t, kx, ky, kz)
% [u, un] = pweval(t, kx, ky) returns the potential and normal derivatives
% evalutions for plane waves.

[M, N]= size(t.x);
u = zeros(M,N); 
u = exp( 1i*( kx * t.x + ky * t.y + kz * t.z));
u = reshape(u, M, N); 
if nargout >1,
    un = zeros(M,N); 
    un = 1i * (kx .* t.nx + ky .* t.ny + kz .* t.nz).* u;
end