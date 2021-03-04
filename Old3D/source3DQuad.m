function source = source3DQuad(f, f_t, N, P, RR, hs)


ss = BORsource3DQuad(f, f_t, N,RR, hs); 
p = (0:P-1)/P*2*pi; 
p = reshape(p, P, 1); 
source.x = kron(ss.r, cos(p)); 
source.y = kron(ss.r, sin(p)); 
source.z = kron(ss.z, ones(P,1)); 

end