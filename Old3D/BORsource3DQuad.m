function sourceBOR = BORsource3DQuad(f, f_t, N, RR, hs)
% BORsource3DQuad(f, f_t, N,RR, hs) returns a structure that contains the locations
% of the MFS source pts in the r-z plane;
%
% Larry Liu, 06/11/2014
if 0,    % doing porportional to local speed
    t=((1:N)-1/2)/N*pi-pi/2;
    %t = (1:N)/N *pi - pi/2;
    %t = t*1.03;
    
    t = reshape(t, N,1);
    
    ss.r = RR * f(t).*cos(t);
    ss.z = RR * f(t).*sin(t);
    
    ss.tr = (-f(t).*sin(t) + f_t(t).*cos(t));
    ss.tz = (f(t).*cos(t) + f_t(t).*sin(t));
    
    
    ll = sqrt(ss.tr.^2 +ss.tz.^2);
    
    ss.nr = ss.tz./ll;
    ss.nz = -ss.tr./ll;
    
    h = hs * ll;
    
    
    sourceBOR.r = ss.r - h.* ss.nr;
    sourceBOR.z = ss.z - h.* ss.nz;
else   % using complexification
    if 0,
        
        t=((1:2*N)-1/2)/(2*N) *2*pi - pi/2;
        %t = (1:N)/N *pi - pi/2;
        %t = t*1.03;
        
        t = reshape(t, 2*N,1);
        
        ss.r = RR * f(t).*cos(t);
        ss.z = RR * f(t).*sin(t);
        
        z = ss.r + 1i * ss.z;
        zz = imagshiftcurve(z,hs);
        sourceBOR.r = real(zz(1:N));
        sourceBOR.z = imag(zz(1:N));
    else
        
        t = 1i * hs + ((1:N)-1/2)/N*pi-pi/2; 
        t = reshape(t, N, 1); 
        t = RR * exp(1i*t).*f(t);
        sourceBOR.r = real(t);
        sourceBOR.z = imag(t);
        
    end
end

