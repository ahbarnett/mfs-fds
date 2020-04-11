% RSKELFM_SV  Solve using recursive skeletonization factorization for MFS.
%
%    For a factorization F = L*D*U, where L and U are rectangular and D is
%    square, this applies the pseudoinverse PINV(F) = U\INV(D)/L, provided that
%    FASTSV = 'N' in RSKELFM. If FASTSV ~= 'N', then an approximate
%    pseudoinverse is used instead.
%
%    Typical complexity: same as RSKELFM_MV.
%
%    Y = RSKELFM_SV(F,X) produces the matrix Y by applying the pseudoinverse of
%    the factored matrix F to the matrix X.
%
%    Y = RSKELFM_SV(F,X,TRANS) computes Y = F\X if TRANS = 'N' (default),
%    Y = F.'\X if TRANS = 'T', and Y = F'\X if TRANS = 'C'.
%
%    See also RSKELFM, RSKELFM_MV.

function Y = rskelfm_sv(F,X,trans)

  % set default parameters
  if nargin < 3 || isempty(trans), trans = 'n'; end

  % check inputs
  trans = chktrans(trans);

  % handle transpose by conjugation
  if trans == 't', Y = conj(rskelfr_sv(F,conj(X),'c')); return; end

  % initialize
  n = F.lvp(end);

  % no transpose
  if trans == 'n'

    % upward sweep in row space
    for i = 1:n
      if F.factors(i).pm == 'r' && ~isempty(F.factors(i).pL)
        psk = F.factors(i).psk;
        prd = F.factors(i).prd;
        tmp = X(psk,:) + F.factors(i).pT*X(prd,:);
        X(psk,:) = F.factors(i).pL'\(F.factors(i).pL\tmp);
      end
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      X(rd,:) = X(rd,:) - F.factors(i).rT'*X(sk,:);
      X(rd,:) = F.factors(i).L\X(rd(F.factors(i).p),:);
      X(sk,:) = X(sk,:) - F.factors(i).E*X(rd,:);
    end

    % transfer from row to column space
    Y = zeros(F.N,size(X,2));
    Y([F.factors.crd],:) = X([F.factors.rrd],:);

    % downward sweep in column space
    for i = n:-1:1
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      Y(rd,:) = Y(rd,:) - F.factors(i).F*Y(sk,:);
      Y(rd,:) = F.factors(i).U\Y(rd,:);
      Y(sk,:) = Y(sk,:) - F.factors(i).cT*Y(rd,:);
      if F.factors(i).pm == 'c' && ~isempty(F.factors(i).pL)
        psk = F.factors(i).psk;
        prd = F.factors(i).prd;
        Y(psk,:) = F.factors(i).pL'\(F.factors(i).pL\Y(psk,:));
        Y(prd,:) = F.factors(i).pT'*Y(psk,:);
      end
    end

  % conjugate transpose
  else

    % upward sweep in column space
    for i = 1:n
      if F.factors(i).pm == 'c' && ~isempty(F.factors(i).pL)
        psk = F.factors(i).psk;
        prd = F.factors(i).prd;
        tmp = X(psk,:) + F.factors(i).pT*X(prd,:);
        X(psk,:) = F.factors(i).pL'\(F.factors(i).pL\tmp);
      end
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      X(rd,:) = X(rd,:) - F.factors(i).cT'*X(sk,:);
      X(rd,:) = F.factors(i).U'\X(rd,:);
      X(sk,:) = X(sk,:) - F.factors(i).F'*X(rd,:);
    end

    % transfer from column to row space
    Y = zeros(F.M,size(X,2));
    Y([F.factors.rrd],:) = X([F.factors.crd],:);

    % downward sweep in row space
    for i = n:-1:1
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      Y(rd,:) = Y(rd,:) - F.factors(i).E'*Y(sk,:);
      Y(rd(F.factors(i).p),:) = F.factors(i).L'\Y(rd,:);
      Y(sk,:) = Y(sk,:) - F.factors(i).rT*Y(rd,:);
      if F.factors(i).pm == 'r' && ~isempty(F.factors(i).pL)
        psk = F.factors(i).psk;
        prd = F.factors(i).prd;
        Y(psk,:) = F.factors(i).pL'\(F.factors(i).pL\Y(psk,:));
        Y(prd,:) = F.factors(i).pT'*Y(psk,:);
      end
    end
  end
end
