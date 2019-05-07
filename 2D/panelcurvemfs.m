function s = panelcurvemfs(p,t,tpan,opts)
% PANELCURVEMFS   MFS source points given a panelization of one closed curve
%
% s = panelcurvemfs(p,t,tpan,opts) generates MFS source points given a
%  boundary panelization of one closed curve.
%
% Inputs:
%  p = number of sources per panel
%  t = surf curve struct including t.Z the parameterization from [0,2pi) to
%      the complex plane.
%  tpan = parameter start points of surface panels
%  opts - optional struct with:
%         opts.reldist - relative distance for straight MFS segments
%                        (negative for interior)
%         opts.meth = 'd' disconnected straight MFS panels
%                     'c' connected straight MFS panels
% Outputs:
%  s = struct with fields: x giving complex locations of sources
%
% Called w/o args, does self-test.

% Barnett, Mar 2019; meths & new interface 4/22/19.

if nargin==0, test_panelcurvemfs; return; end
if nargin<4, opts = []; end
if ~isfield(opts,'reldist'), opts.reldist = 0.4; end
if ~isfield(opts,'meth'), opts.meth='c'; end
Np = numel(tpan);   % # panels
tpan = tpan(:).';   % make row vec

mfspans = nan(2,Np);   % build MFS source panel starts and ends as C #s...
zpan = t.Z(tpan);      % get panel start locs - only time t.Z used for now.
zpane = circshift(zpan,-1);   % end locs of each surf panel
v = zpane-zpan;        % surf panel start-to-end vectors, as row of C #s
n = -1i * (circshift(v,1)+v)/2;    % crude est outward normals
n = n./abs(n);                     % unit normals
if opts.meth=='d'
  mfspans = [zpan;zpane] + ([v;v]/1i) * opts.reldist;  % transl in normal dir
elseif opts.meth=='c'
  l = abs(v);    % surf panel lengths
  lpts = (circshift(l,1)+l)/2;     % mean lengths at the start pts of surf pans
  mfspans(1,:) = zpan + opts.reldist*n.*lpts;
  mfspans(2,:) = circshift(mfspans(1,:),-1);  % end pt is next start pt (closes)
else
  error(sprintf('unknown opts.meth: %s\n',opts.meth))
end
  
s.x = nan(p*Np,1);  % now build list of src pts, via unif linear panels...
for i=1:Np
  s.x((i-1)*p+(1:p)) = mfspans(1,i) + (mfspans(2,i)-mfspans(1,i))*(0.5:p-0.5)/p;
end

%%%%%%%
function test_panelcurvemfs   % just plot vis test fow now
a = 0.3; w = 3; b = 0.5; t.Z = @(t) (1+a*cos(w*t+b)).*exp(1i*t); % trefoil
Np = 20;
tpan = 2*pi*(1:Np)/Np;
p = 10;
o.meth = 'c';
o.reldist = -0.4;      % <0 for interior src
s = panelcurvemfs(p,t,tpan,o);
zpan = t.Z(tpan);
figure; n=1e3; plot(t.Z(2*pi*(1:n)/n),'-'); hold on; plot(zpan,'+');
plot(s.x,'r.'); axis equal tight; title('test panelcurvemfs');
