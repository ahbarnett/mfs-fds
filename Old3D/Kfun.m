% kernel function
function K = Kfun(x,y,lp,nu)
% nu is the source normal used for double layer potential for d only; 
global k
if nargin < 4
    nu = [];
end
dx = bsxfun(@minus,x(1,:)',y(1,:));
dy = bsxfun(@minus,x(2,:)',y(2,:));
dz = bsxfun(@minus,x(3,:)',y(3,:));
dr = sqrt(dx.^2 + dy.^2 + dz.^2);
if strcmpi(lp,'s')
    K = 1/(4*pi) * exp(1i * k * dr)./dr;
elseif strcmpi(lp,'d')
    rdotn = bsxfun(@times,dx,nu(1,:)) + bsxfun(@times,dy,nu(2,:)) + ...
        bsxfun(@times,dz,nu(3,:));
    K = 1/(4*pi).*rdotn./dr.^3;
end
K(dr == 0) = 0;
end