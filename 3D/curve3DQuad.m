function target = curve3DQuad(f, f_t, N, NP, RR)
% target = curvequad(N, RR) returns an object which has x, y, z, nx, ny, nz
% properties, which are the coordiantes and normal derivatives of the curve
% quadrature points. dimension of target.x is N^2 x 1; 
%
% Larry Liu, May 16, 2014 


t=((1:N)-1/2)/N*pi-pi/2;    
%t = (1:N)/N*pi - pi/2; 
p = (0:NP-1)/NP*2*pi;  % key, start from 0
[T, P] = meshgrid(t,p); 
tt = reshape(T, N*NP,1); 
pp = reshape(P, N*NP,1); 

target.x=RR * f(tt).*cos(tt).*cos(pp);
target.y=RR * f(tt).*cos(tt).*sin(pp);
target.z=RR * f(tt).*sin(tt);

target.tx = (-f(tt).*sin(tt) + f_t(tt).*cos(tt)) .* cos(pp) ; 
target.ty = (-f(tt).*sin(tt) + f_t(tt).*cos(tt)) .* sin(pp) ; 
target.tz = (f(tt).*cos(tt) + f_t(tt).*sin(tt));

target.px = -f(tt).*cos(tt).*sin(pp);
target.py = f(tt).*cos(tt).*cos(pp);
target.pz = zeros(N*NP,1);

aa= -cross ([target.tx, target.ty, target.tz], [target.px, target.py, target.pz], 2); 
target.nx = aa(:,1);  
target.ny = aa(:,2); 
target.nz = aa(:,3); 
ll = sqrt(target.nx.^2 + target.ny.^2 + target.nz.^2); 
target.nx = target.nx./ll; 
target.ny = target.ny./ll; 
target.nz = target.nz./ll; 

       
