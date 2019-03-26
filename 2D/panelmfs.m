function s = panelmfs(p,t,ts,pans,opts)
% PANELMFS   choose source points in MFS given a panelization of bdry.
%
% s = panelmfs(p,pans,reldist)
% Inputs:
%  p = number of sources per panel
%  t = surface panel struct with t.Z analytic formula
%  pans = starts and ends of surface panels
%  opts - optional struct with:
%         opts.reldist - relative distance for straight MFS segments
%                        (negative for interior)

if nargin<5, opts = []; end
if ~isfield(opts,'reldist'), opts.reldist = 0.4; end
Np = numel(ts);   % # panels

s.x = nan(p*Np,1);
% crude linear source segments for now...
for i=1:Np
  v = pans(2,i)-pans(1,i);   % panel start-to-end vector
  mfspan = pans(:,i) + (v/1i) * opts.reldist;
  s.x((i-1)*p+(1:p)) = mfspan(1) + v*(0.5:p-0.5)/p;
end
