% driver for Helmholtz exterior scattering from polygon.
% Barnett 5/6/19

clear
v=1;
k = 15;
eps = 1e-10;

% make rand poly w/ right angles only, in [-1,1]^2
nb = 100;  % # bars
ysc = 10/nb;
h = 2/nb;       % horiz width of bars
nc = (nb+1)*2;    % # corners
c(1:2) = [-1-.2i*ysc, 1-.2i*ysc];    % bottom corners
for b=1:nb, c(b*2+(1:2)) = 1i*rand*ysc+1-h*[b-1,b]; end

figure; plot([c c(1)], '.-');
maxL = h;   % should also respect k (at given eps)
o=[]; o.geomeps = 1e-14;
[t s] = build_polygon(c,'e',maxL,eps,o);
figure; plot(t.x,'.'); hold on; plot(s.x, '.'); axis equal
plot([c;c+0.2*h*t.cornors],'r-');  % show mean normals at corners
for i=1:nc   % label corner MFS params
  text(real(c(i)),imag(c(i)),sprintf('\\alpha=%.3g r=%.2g Nr=%d d=%.2g',t.angs(i), t.reffac(i), t.Nr(i), t.reldist(i)))
end
