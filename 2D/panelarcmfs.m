function s = panelarcmfs(p,t,tpan,opts)
% PANELARCMFS.  set up MFS point panels for one arc, matching possible corners
%
% s = panelarcmfs(p,t,tpan,opts) generates MFS source points given a
%  boundary panelization of one open arc.
%
% Inputs:
%  p = number of sources per panel
%  t = arc struct including t.Z the parameterization from [0,1) to
%      the complex plane.
%  tpan = parameter start points of surface panels, in [0,1)
%  opts - optional struct with:
%         opts.reldist - relative distance of sources from panels
%                        (negative for interior). If two values are given,
%                        it will interpolate between these from one end to other
%         opts.cornors - average outward unit normals at the two corners, helps
%                        make better MFS sources
%         opts.meth = 'd' disconnected straight MFS panels
%                     'c' connected straight MFS panels
% Outputs:
%  s = struct with fields: x giving complex locations of sources
%
% Called w/o args, does self-test.
%
% Notes: no fancy relative coordinates for now; too much refinement leads to
% coincident points to machine precision.

% Barnett 5/6/19

if nargin==0, test_panelarcmfs; return; end
if nargin<4, opts = []; end
if ~isfield(opts,'reldist'), opts.reldist = [-0.4]; end
if numel(opts.reldist)==1, opts.reldist = opts.reldist*[1 1]; end
if ~isfield(opts,'meth'), opts.meth='c'; end
if ~isfield(opts,'cornors')     % fake some normals at the arc's two ends
  cn(1)=diff(t.Z(tpan([1 2]))); cn(2)=diff(t.Z([tpan(end) 1])); cn
  opts.cornors = -1i*cn./abs(cn);
end
Np = numel(tpan);           % # panels
tpan = [tpan(:).' 1.0];     % make row vec, with final endpt, length Np+1

zpan = t.Z(tpan(1:end-1));  % get panel start locs - only time t.Z used for now
zpane = [zpan(2:end) t.Z(1)];  % ends of each surf panel
v = zpane-zpan;             % surf panel start-to-end vectors, as row of C #s
n = -1i*(v(2:end)+v(1:end-1))/2;   % crude est outward normals at interior juncs
n = [opts.cornors(1) n opts.cornors(2)];   % append normals at two arc ends
n = n./abs(n);
mtpan = (tpan(1:end-1) + tpan(2:end))/2;      % mean param values of each pan
if opts.meth=='d'
  reld = opts.reldist(1)*(1.0-mtpan) + opts.reldist(2)*mtpan;   % lin interp
  mfspans = [zpan;zpane] + ([v;v]/1i).*reld;       % transl in normal dir
elseif opts.meth=='c'
  l = abs(v);    % surf panel lengths
  lpts = [l(1) (l(1:end-1)+l(2:end))/2 l(end)];  % typ pan lens at Np+1 brkpts
  reld = opts.reldist(1)*(1.0-tpan) + opts.reldist(2)*tpan;   % lin interp
  mfsends = [zpan zpane(end)] + reld.*n.*lpts;
  mfspans = nan(2,Np);   % build MFS source panel starts and ends as C #s...
  mfspans(1,:) = mfsends(1:end-1);  % connect each pan from one end pt to next
  mfspans(2,:) = mfsends(2:end);
else
  error(sprintf('unknown opts.meth: %s\n',opts.meth))
end
  
s.x = nan(p*Np,1);  % now build list of src pts, via unif linear panels...
for i=1:Np
  s.x((i-1)*p+(1:p)) = mfspans(1,i) + (mfspans(2,i)-mfspans(1,i))*(0.5:p-0.5)/p;
end

%%%%%%%
function test_panelarcmfs   % just plot vis test fow now
a = 5;
t.Z = @(t) exp(1i*a*(t-2)); t.Zp = @(t) 1i*a*exp(1i*a*(t-2));  % Z, Z' in [0,1]
[t tpan zpan] = panelarc(t,12,3,[5 10],[3 2]);
p = 10;
for m = 'dc', o.meth = m;      % loop over methods
  o.reldist = [-0.5 -0.2];      % <0 for interior src
  s = panelarcmfs(p,t,tpan,o);
  zpan = t.Z(tpan);
  figure; n=1e3; plot(t.Z((0:n)/n),'-'); hold on; plot(zpan,'+');
  plot(s.x,'r.'); axis equal tight;
  title(sprintf('test panelarcmfs, meth=%s',o.meth))
end
