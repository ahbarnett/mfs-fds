function [reffac Nr reldist] = cornermfsrules(ang, eps, lpan, opts)
% CORNERMFSRULES.  recipe for MFS parameters to handle a corner to accuracy eps
%
% [reffac Nr reldist] = cornermfsrules(ang, eps, lpan)
%  returns refinement factor, number of panel refinements, and relative distance
%  parameter to handle a corner of angle ang via the MFS to accuracy eps. ang
%  is measured on the physical side, ie, ang>pi is reentrant (tough) corner.
%  lpan is the physical length of the base panel, used to decide when to clip Nr
%
% [reffac Nr reldist] = cornermfsrules(ang, eps, lpan, opts)
%  also controls opts such as:
%       opts.geomeps - smallest allowed panel size
%
%  Notes:
%  * ad-hoc for now; really such rules depend on singularity of PDE
%  * it's assumed p for MFS (and p for surface collocation) are set to match
%    eps.
%  * reldist>0; needs to be negated for exterior case
%
% See also: scat_panel_onecorner_conv_demos, where the recipes developed.

% Barnett 5/7/19

if nargin<3, opts = []; end
% prevent panels smaller than this size...
if ~isfield(opts,'geomeps'), opts.geomeps = 1e-14; end

reldist = 0.5;       % 0.5, defaults
reffac = 3.0;      % was 3.0, but needed 2.0 for skyline.
if ang>3*pi/2        % "bad" reentrant corner. As approach 2pi, go...
  reldist = 0.5 - 0.14*(ang-3*pi/2)^2;   % closer
end
if ang>pi            % "bad" reentrant corner. As approach 2pi, go...
  reffac = max(1.5, reffac - 0.2*(ang-pi)^2);     % slower refinement
end
Nr = ceil( (ang/pi) * log(1/eps)/log(reffac) );    % assuming Laplace BVP sing
Nr = min(Nr, floor(log(lpan/opts.geomeps)/log(reffac)) );  % prevent too close
