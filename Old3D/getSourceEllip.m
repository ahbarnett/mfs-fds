function s = getSourceEllip(a, b, c, d, N1, N2)

t = getShapeEllip(a, b, c, N1, N2); 

s.x = t.x - d * (sqrt(t.ls.^2 + t.lt.^2)) .* t.nx; 
s.y = t.y - d * (sqrt(t.ls.^2 + t.lt.^2)) .* t.ny; 
s.z = t.z - d * (sqrt(t.ls.^2 + t.lt.^2)) .* t.nz; 



