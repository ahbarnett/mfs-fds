function b = QuasiPerVecFilling3D(t, kx, ky, kz) 
% 
% b = QuasiPerVecFilling3D(t, kx, ky, kz) returns the first part of RHS vector
% for 3D acoustic scattering over periodic structures; it is just the
% negative of the incident wave; 
% 
% It is used when we want to separate A, B, C, Q from the total matrix a
% and want to only consider the first part of the RHS vector
% 
% Larry Liu, 06/09/2014

u = pweval(t, kx, ky, kz); 

b = - u;  