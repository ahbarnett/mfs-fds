function t = getShapeEllip(a, b, c, N1, N2)

tt=((1:N2)-1/2)/N2*2*pi;
ss=gauss(N1)*pi/2+pi/2; 
[TT, SS] = meshgrid (tt, ss);
tt = reshape(TT, N1*N2, 1); 
ss = reshape(SS, N1*N2, 1); 

t.x = a*sin(ss).*cos(tt); 
t.y = b*sin(ss).*sin(tt); 
t.z = c*cos(ss); 

t.sx = a*cos(ss).*cos(tt); 
t.sy = b*cos(ss).*sin(tt);
t.sz = -c*sin(ss); 

t.tx = -a*sin(ss).*sin(tt); 
t.ty = b*sin(ss).*cos(tt); 
t.tz = zeros(N1 * N2, 1);


t.ls = sqrt(t.sx.^2 + t.sy.^2 + t.sz.^2); 
t.lt = sqrt(t.tx.^2 + t.ty.^2 + t.tz.^2); 

aa= cross ([t.sx, t.sy, t.sz], [t.tx, t.ty, t.tz], 2); 
t.nx = aa(:,1);  
t.ny = aa(:,2); 
t.nz = aa(:,3); 
ll = sqrt(t.nx.^2 + t.ny.^2 + t.nz.^2); 
t.nx = t.nx./ll; 
t.ny = t.ny./ll; 
t.nz = t.nz./ll; 

