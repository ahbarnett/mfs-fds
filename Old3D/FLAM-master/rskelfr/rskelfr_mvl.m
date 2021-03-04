% RSKELFR_MVL  Multiply by L factor in rectangular recursive skeletonization
%              factorization F = L*D*U.
%
%    See also RSKELFR, RSKELFR_MV.

function Y = rskelfr_mvl(F,X,trans)

  % set default parameters
  if nargin < 3 || isempty(trans)
    trans = 'n';
  end

  % check inputs
  assert(strcmpi(trans,'n') || strcmpi(trans,'c'), ...
         'FLAM:rskelfr_mvl:invalidTrans', ...
         'Transpose parameter must be either ''N'' or ''C''.')

  % initialize
  n = F.lvp(end);
  Y = X;

  % no transpose
  if strcmpi(trans,'n')
    for i = n:-1:1
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      Y(sk,:) = Y(sk,:) + F.factors(i).E*Y(rd,:);
      Y(rd,:) = Y(rd,:) + F.factors(i).rT*Y(sk,:);
    end

  % conjugate transpose
  else
    for i = 1:n
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      Y(sk,:) = Y(sk,:) + F.factors(i).rT'*Y(rd,:);
      Y(rd,:) = Y(rd,:) + F.factors(i).E'*Y(sk,:);
    end
  end
end