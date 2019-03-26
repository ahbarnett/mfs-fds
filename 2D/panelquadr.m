function [s ts pans] = panelquadr(s,p,Nb,Nr,r,egap,opts)
% PANELQUADR - set up toy 2D panel quadratures for smooth curve or one corner
%
% s = panelquadr(s,p,Nb) fills struct s with quadrature nodes (in complex
%  plane) and other geometric info.
% [s ts pan0 pan1] = panelquadr(s,p,Nb) also outputs:
%   ts: parameters values of panel breakpts.
%   pans: 2-by-Np complex coords of panel starts and ends
%
% Inputs:
% s = curve struct with s.Z param funcs on [-pi,pi].
% p = nodes per panel
% Nb = number of "base" panels before refinement
%
% s = panelquadr(s,p,Nb,Nr,r,egap) does refinement about parameter t=0 only:
% Nr = # refinement levels (=0 no refinement)
% r = geometric refinement factor (2=dyadic; but 3-5 is more efficient!)
% egap : 0=no gap, 1=2eps-sized gap in r-adic refinement
% If any inputs are absent or empty, default values are used.
%
% s = panelquadr(s,p,Nb,Nr,r,egap,opts) controls options such as:
%   opts.qtype = 'g' (Gauss-Legendre, default) vs 'u' (uniform)
%
% Barnett 10/30/14, simplified 3/26/19

if nargin<2 || isempty('p'), p=16; end
if nargin<3 || isempty('Nb'), Nb=20; end
if nargin<4 || isempty('Nr'), Nr=0; end
if nargin<5 || isempty('r'), r=3; end
if nargin<5 || isempty('egap'), egap=0; end
if nargin<6, opts=[]; end
if ~isfield(opts,'qtype'), opts.qtype='g'; end

if Nr>0, meth='d'; else meth='n'; end   % whether to refine

Nb = ceil(Nb/2)*2; % make even # base panels so a split occurs at t=0.
% Set up r-adic panels...
ts = [-Nb/2:Nb/2-1]/Nb*2*pi; % param starts of base panels, in param order
tl = 0*ts + 2*pi/Nb; % param lengths of base panels
if strcmp(meth,'d') % set up r-adic refinement via longer ts, tl lists...
  tr = r.^(0:-1:-Nr)*tl(1); % refined panel starts
  ts = [ts(1:Nb/2-1), -tr, 0, tr(end:-1:1), ts(Nb/2+3:end)];
  tl = diff([ts ts(1)+2*pi]);
  if egap, ts = [ts(1:Nb/2-1+Nr), ts(Nb/2+2+Nr:end)];
    tl = [tl(1:Nb/2-1+Nr), tl(Nb/2+2+Nr:end)]; end % remove inner panels
end

if opts.qtype=='g'
  [tp wp] = gauss(p);      % nodes per panel; set up panels from ts,tl
elseif opts.qtype=='u'
  tp = -1 + 2*(0.5:p-0.5)/p; wp = (2/p)*ones(1,p);
end
Np = numel(ts); t = zeros(p,Np); w = t; % Fill p*Np matrices: t param vals...
pans = nan(2,Np,1);    % panel starts and ends in complex plane
for i=1:Np
  t(:,i) = ts(i) + (tp+1)/2*tl(i); w(:,i) = wp/2*tl(i);   % param pts
  pans(:,i)=[s.Z(ts(i)); s.Z(ts(i)+tl(i))];
end
x = s.Z(t);
if isfield(s,'Zp')
  nx = -1i * s.Zp(t); sp = abs(nx);
  if isfield(s,'Zpp')
    kap = -real(conj(nx).*s.Zpp(t))./sp.^3;
    s.cur = kap(:);
  end
  nx = nx./sp;   % x bdry pts, xn normals, sp speed, kap curvature
  s.w = sp(:).*w(:); s.nx = nx(:); 
end
% load into seg struct
s.x = x(:); s.t = t(:);
