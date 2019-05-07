function [t s] = build_polygon(c,side,maxL,eps,opts)
% BUILD_POLYGON  make panelized closed polygon with bdry and MFS src panels
%
% [t s] = build_polygon(c,side,maxL,eps,opts)
%
% Inputs:
%  c - cyclic list of corners as pts in complex plane, in CCW sense
%  side - 'i' for interior BVP (ext MFS), or 'e' for ext BVP (int MFS pts)
%  maxL - maximum length of panels (crude global refinement control)
%  eps - desired accuracy
%
% Outputs:
%  t - surface struct including t.x nodes
%  s - MFS source curve struct including s.x source locations
%
% Called w/o arguments, does a self-test.

% Barnett 5/7/19

if nargin==0, test_build_polygon; return; end

nc = numel(c);           % # corners
v = circshift(c,-1)-c;   % vectors from corner to next corner in CCW sense.
v = v./abs(v);           % unit tangent vecs fwd from each corner
angs = angle(v./circshift(v,1)) + pi;    % exterior angle at each corner
cornors = v.*sqrt(circshift(v,1)./v)/1i;  % mean normals at corners
reldist = nan(nc,1); reffac = reldist; Nb = reffac; Nr = Nb;     % alloc
for i=1:nc                   % choose all MFS params, per corner
  [reffac(i) Nr(i) reldist(i)] = cornermfsrules(angs(i), eps, opts);
end

pmfs = 2 + ceil(log10(1/eps));         % MFS panel orders
p = pmfs + 2;                          % surf colloc panel orders

t.x = [];
s.x = [];
for i=1:nc          % loop over arcs: arc 1 connects corner 1 to 2, etc
  i2 = i+1; if i2>nc, i2=1; end    % ind of next corner, cyclically
  tt.Z = @(t) c(i)*(1-t) + c(i2)*t; tt.Zp = @(t) c(i2)-c(i);   % straight arc
  Nb(i) = max(2,ceil(abs(c(i2)-c(i))/maxL));        % limit max pan length
  [tt tpan] = panelarc(tt,p,Nb(i),Nr([i i2]),reffac([i i2])); % MFS pars @ corns
  o = []; o.reldist = reldist([i i2]);        % MFS params at both end corns
  if side=='e', o.reldist = -o.reldist; end
  o.cornors = cornors([i i2]);
  ss = panelarcmfs(pmfs,tt,tpan,o);
  t.x = [t.x; tt.x]; s.x = [s.x; ss.x];  % append to big col vecs
end

t.cornors = cornors;   % pass out some diagnostics
t.angs = angs; t.reffac = reffac; t.Nr = Nr; t.Nb = Nb; t.reldist = reldist;
% ***  todo t.inside = ...

%%%%
function test_build_polygon
sides = 5; maxL = 0.2; eps = 1e-10; opts.geomeps = 1e-14;
c = exp(1i*2*pi*((0:sides-1)/sides).^2);
[t s] = build_polygon(c,'e',maxL,eps,opts);
figure; plot(t.x,'.'); hold on; plot(s.x, '.'); axis equal
plot([c;c+0.2*t.cornors],'r-');
m = (c+circshift(c,-1))/2;      % midpoints of arcs
for i=1:sides
  text(real(c(i)),imag(c(i)),sprintf('\\alpha=%.3g r=%.2g Nr=%d d=%.2g',t.angs(i), t.reffac(i), t.Nr(i), t.reldist(i)))
  text(real(m(i)),imag(m(i)),sprintf('Nb=%d',t.Nb(i)));
end
