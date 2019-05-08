function [s tpan zpan] = panelarc(s,p,Nb,Nrs,reffacs,opts)
% PANELARC.  set up boundary point panels covering one arc connecting 2 corners
%
% s = panelarc(s,p,Nb,Nrs,reffacs,opts)
%
% Inputs:
%  s: curve struct with s.Z complex param funcs on [0,1], and s.Zp its complex
%     deriv, if weights are to be needed, and s.Zpp its 2nd deriv, if curvature
%     needed.
%  p: nodes per panel
%  Nb: number of equally-spaced "base" panels before refinement
%  Nrs: 2-element vector giving # panel refinements on corner at 0 and at 1.
%  reffacs: 2-element vector giving refinement ratios at corner at 0 and at 1.
%  opts: (optional) struct with optional fields:
%        opts.qtype = 'u' (uniform on panels) or 'g' (Gauss-L on panels).
%
% Outputs:
%   s: struct with new boundary discretization fields including
%         x - nodes as complex numbers
%         w - weights (speed weights), real numbers
%   tpan: parameter values of start of each panel
%   zpan: 2-by-Np complex coords of panel starts and ends
%
% Called w/o arguments, does a simple self-test.
%
% Notes: no fancy relative coordinates for now; Nrs too high leads to
% coincident points to machine precision.

% Barnett 5/6/19

if nargin==0, test_panelarc; return; end
if nargin<2 || isempty('p'), p=12; end
if nargin<3 || isempty('Nb'), Nb=5; end
if nargin<4 || isempty('Nrs'), Nr=[0 0]; end   % no corner refinements
if nargin<5 || isempty('reffacs'), reffacs=[3 3]; end
if nargin<6, opts=[]; end
if ~isfield(opts,'qtype'), opts.qtype='g'; end

if opts.qtype=='g'        % the std panel quad on [-1,1]: nodes tp, wei wp
  [tp wp] = gauss(p);
elseif opts.qtype=='u'
  tp = -1 + 2*(0.5:p-0.5)/p; wp = (2/p)*ones(1,p);
end

Np = Nb + sum(Nrs);       % total # pans
L = 1.0/Nb;               % param length of each base panel
tpan = zeros(1,Np);         % ts(1)=0; fill other start points in order...
tpan(2:Nrs(1)+1) = L*reffacs(1).^-(Nrs(1):-1:1);   % first corner at 0
tpan(Nrs(1)+(2:Nb-1)) = L*(1:Nb-2);             % middle base panels, arithm
ic2 = Np+(-Nrs(2):0);     % indices for panels at 2nd corner
tpan(ic2) = -L*reffacs(2).^-(0:Nrs(2));   % 2nd corner; don't offset to 1 yet
tl = diff([tpan 0]);      % since all 0-offset, get rel acc
tpan(ic2) = 1.0 + tpan(ic2);  % fix the offset for 2nd corner stuff
tl(ic2(1)-1) = tpan(ic2(1))-tpan(ic2(1)-1);   % fix the one wrong length
s.t = zeros(p,Np); s.w = s.t;     % allocate output nodes, weights
for i=1:Np                % output nodes and weights
  s.t(:,i) = tpan(i) + (tp+1)/2*tl(i);
  s.w(:,i) = wp/2*tl(i);
  zpan(:,i) = [s.Z(tpan(i)); s.Z(tpan(i)+tl(i))];
end
s.t = s.t(:);             % one big col vec
s.x = s.Z(s.t);
if isfield(s,'Zp')
  nx = -1i * s.Zp(s.t);
  sp = abs(nx);
  if isfield(s,'Zpp')
    kap = -real(conj(nx).*s.Zpp(s.t))./sp.^3;
    s.cur = kap(:);
  end
  nx = nx./sp;   % x bdry pts, nx normals, sp speed, kap curvature
  s.w = sp(:).*s.w(:);
  s.nx = nx(:);
end

%%%%%
function test_panelarc
a = 5;
t.Z = @(t) exp(1i*a*(t-2)); t.Zp = @(t) 1i*a*exp(1i*a*(t-2));  % Z, Z' in [0,1]
[t tpan zpan] = panelarc(t,12,3,[5 10],[3 2]);
figure; plot(t.x,'.'); axis equal tight; hold on; plot(t.Z(tpan),'*');
plot([t.x t.x+0.1*t.nx].', 'c-'); plot(zpan,'r-'); title('test panelarc');
fprintf('rel err in total wei: %.3g\n',(sum(t.w)-a)/a)


