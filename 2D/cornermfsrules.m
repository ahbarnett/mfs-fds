function [reffac Nr reldist] = cornermfsrules(ang, eps, opts)
% CORNERMFSRULES.  recipe for MFS parameters to handle a corner to accuracy eps
%
% [reffac Nr reldist] = cornermfsrules(ang, eps)
%  returns refinement factor, number of panel refinements, and relative distance
%  parameter to handle a corner of angle ang via the MFS to accuracy eps. ang
%  is measured on the physical side, ie, ang>pi is reentrant (tough) corner.
%
% [reffac Nr reldist] = cornermfsrules(ang, eps, opts)
%  also controls opts such as:
%     opts.geomeps - smallest allowed panel size
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


reldist = 0.5;       % defaults
reffac = 3.0;
if ang>3*pi/2        % "bad" reentrant corner. As approach 2pi, go...
  reldist = 0.5 - 0.14*(ang-3*pi/2)^2;   % closer
  reffac = 3.0 - 0.4*(ang-3*pi/2)^2;     % slower refinement
end
Nr = ceil( (ang/pi) * log(1/eps)/log(reffac) );     % assuming Laplace BVP
Nr = min(Nr, floor(log(1/opts.geomeps)/log(reffac)) );
