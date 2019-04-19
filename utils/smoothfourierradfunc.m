function [R Rp Rpp] = smoothfourierradfunc(aj,bj)
% SMOOTHFOURIERRADFUNC.  general Fourier series smooth radial function
%
%  [R Rp Rpp] = SMOOTHFOURIERRADFUNC(aj, bj) returns function handles to
%   R(theta) polar function defining a smooth closed curve with Fourier
%   series coeffs aj (cosine) and bj (sine), for 1=1,..,nterms.
%   The zero cosine (DC) term is fixed at 1.  aj and bj may be row or col
%   vectors, but numel(bj) must be at least numel(aj).
%   Rp and Rpp are the 1st and 2nd analytic theta-derivative function handles.

% Barnett 4/18/19 based on smoothfourier from MPSpack.
n = numel(aj);
if numel(bj)~=n, error('bj must have n=numel(aj) elements!'); end
aj = aj(:).'; bj = bj(:).';  % make row vecs
j = (1:n).';                 % col vec of indices
if n==0    % circle case
  R = @(t) 1+0*t; Rp = @(t) 0*t; Rpp = @(t) 0*t;
else       % now use outer products inside the trig funcs...
  R = @(t) 1 + reshape(aj*cos(j*t(:).') + bj*sin(j*t(:).'),size(t));
  Rp = @(t) reshape(-(j'.*aj)*sin(j*t(:).') +(j'.*bj)*cos(j*t(:).'),size(t));
  Rpp = @(t) reshape(-(j'.^2.*aj)*cos(j*t(:).') - (j'.^2.*bj)*sin(j*t(:).') ,size(t));
end
