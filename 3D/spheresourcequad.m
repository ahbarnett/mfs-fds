function [xs, ys, zs] = spheresourcequad(N, R)
% [xs, ys, zs] = spheresourcequad(N, R) returns the x, y, z coordiante of
% the sphere with radius R and N^2 points; 
% 
% By Larry Liu, May 16, 2014 

if 0,    % uniform in z-direction
	z = ((1:N)-1/2)/N * 2*R -R ;   p = (1:N)/N*2*pi; 
	[Z, P] = meshgrid(z,p); 
	zz = reshape(Z, N^2,1); 
	pp = reshape(P, N^2,1); 
	xs = sqrt(R^2 - zz.^2).* cos(pp); 
	ys = sqrt(R^2 - zz.^2).* sin(pp); 
	zs = zz; 
	
else,  % uniform in theta - direction

	t = ((1:N)-1/2)/N * pi - pi/2;   p = (1:N)/N*2*pi; 
	[T, P] = meshgrid(t,p); 
	tt = reshape(T, N^2,1); 
	pp = reshape(P, N^2,1); 
	xs = R * cos(tt).* cos(pp); 
	ys = R * cos(tt).* sin(pp); 
	zs = R * sin(tt); 
end
