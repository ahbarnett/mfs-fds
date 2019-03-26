function i = insideradfunc(R,z)
% INSIDERADFUNC   boolean array of whether points inside a radial function
%
% i = insideradfunc(R,z)
%  each element of i is true if the corresponding element of z lies inside
%  the curve r = R(theta) in the complex plane.

% Barnett 3/26/19
i = abs(z) < R(angle(z));
