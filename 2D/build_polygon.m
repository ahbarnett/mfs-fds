function [t s] = build_polygon(c,side,maxL,eps,opts)
% BUILD_POLYGON  make panelized 2D closed polygon with bdry and MFS src panels
%
% [t s] = build_polygon(c,side,maxL,eps,opts) builds surf quadrature and MFS
%  source points for a closed polygon in 2D.
%
% Inputs:
%  c - cyclic list of corners as pts in complex plane, in CCW sense
%  side - 'i' for interior BVP (ext MFS), or 'e' for ext BVP (int MFS pts)
%  maxL - maximum length of panels (crude global refinement control)
%  eps - desired accuracy
%  opts (optional) - gets passed to panelarc, and optional fields control:
%       opts.p - overrides p for surf panel order
%
% Outputs:
%  t - surface struct including t.x nodes, t.inside(z) func, ...
%  s - MFS source curve struct including s.x source locations
%
% Called w/o arguments, does a self-test with labelled plot.

% Barnett 5/7/19

if nargin==0, test_build_polygon; return; end
if nargin<5, opts = []; end

nc = numel(c);           % # corners
v = circshift(c,-1)-c;   % vectors from corner to next corner in CCW sense.
l = abs(v);              % arc lengths
v = v./l;                % unit tangent vecs fwd from each corner
angs = angle(v./circshift(v,1)) + pi;    % exterior angle at each corner
cornors = v.*sqrt(circshift(v,1)./v)/1i;  % mean normals at corners

pmfs = 2 + ceil(log10(1/eps));         % MFS panel orders
p = pmfs + 3;                          % surf colloc panel orders
if isfield(opts,'p'), p = opts.p; end  % override p

Nb = nan(nc,1); t.reffac = Nb; t.reldist = Nb; t.Nr=nan(nc,2);     % alloc
t.x = []; s.x = [];    % get ready to append
for i=1:nc            % loop over arcs: arc 1 connects corner 1 to 2, etc
  i2 = i+1; if i2>nc, i2=1; end    % ind of next corner, cyclically
  tt.Z = @(t) c(i)*(1-t) + c(i2)*t; tt.Zp = @(t) c(i2)-c(i);   % straight arc
  Nb(i) = max(2,ceil(abs(c(i2)-c(i))/maxL));     % upper limit max pan length
  lpan = l(i)/Nb(i);    % base pan len in this arc, for stopping Nr's
  o = [];               % get MFS params at each corner of this arc...
  [reffac(1) Nr(1) o.reldist(1)] = cornermfsrules(angs(i), eps, lpan, opts);
  [reffac(2) Nr(2) o.reldist(2)] = cornermfsrules(angs(i2), eps, lpan, opts);
  t.reffac(i) = reffac(1); t.reldist(i) = o.reldist(1);  % save diagnostics
  t.Nr(i,2)=Nr(1); t.Nr(i2,1)=Nr(2);                     % save Nr's per corn
  [tt tpan] = panelarc(tt,p,Nb(i),Nr,reffac,opts);  % MFS pars @ both corns
  if side=='e', o.reldist = -o.reldist; end
  o.cornors = cornors([i i2]);
  ss = panelarcmfs(pmfs,tt,tpan,o);
  t.x = [t.x; tt.x]; s.x = [s.x; ss.x];  % append to big col vecs
end

t.inside = @(z) inpolygon(real(z),imag(z),real(c),imag(c));
t.cornors = cornors; t.p = p; t.pmfs = pmfs;   % pass out more diagnostics...
t.angs = angs; t.Nb = Nb;

%%%%
function test_build_polygon
sides = 5; maxL = 0.2; eps = 1e-10; opts.geomeps = 1e-14;
c = exp(1i*2*pi*((0:sides-1)/sides).^2);
[t s] = build_polygon(c,'e',maxL,eps,opts);
figure; plot(t.x,'.'); hold on; plot(s.x, '.'); axis equal
plot([c;c+0.2*t.cornors],'r-');
[xx yy] = meshgrid(-1.2:0.05:1.2); zz = xx(:)+1i*yy(:);   % grid of pts
ii = ~t.inside(zz); plot(zz(ii),'.','color',.8*[1 1 1]);  % show outside ones
m = (c+circshift(c,-1))/2;      % midpoints of arcs for labeling
for i=1:sides
  text(real(c(i)),imag(c(i)),sprintf('\\alpha=%.3g r=%.2g Nr=(%d,%d) d=%.2g',t.angs(i), t.reffac(i), t.Nr(i,1), t.Nr(i,2), t.reldist(i)))
  text(real(m(i)),imag(m(i)),sprintf('Nb=%d',t.Nb(i)));
end
