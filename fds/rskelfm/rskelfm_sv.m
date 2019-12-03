% RSKELFM_SV  Solve using recursive skeletonization factorization for MFS.
%
%    For a factorization F = L*D*U, where L and U are square and D is
%    rectangular, this applies the "factored pseudoinverse"
%    INV(U)*PINV(D)*INV(L) consisting of the pseudoinverses of each factor.
%    While this is not a true pseudoinverse, we will continue to use the same
%    notation when referring to it for convenience.
%
%    Typical complexity: same as RSKELFM_MV.
%
%    Y = RSKELFM_SV(F,X) produces the matrix Y by applying the factored
%    pseudoinverse of the factored matrix F to the matrix X.
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
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      X(rd,:) = X(rd,:) - F.factors(i).rT*X(sk,:);
      if i < n
        X(rd,:) = F.factors(i).L\X(rd(F.factors(i).p),:);
        X(sk,:) = X(sk,:) - F.factors(i).E*X(rd,:);
      end
    end

    % transfer from row to column space
     Y = zeros(F.N,size(X,2));
    for i = 1:n
      rrd = F.factors(i).rrd;
      crd = F.factors(i).crd;
      if i < n, Y(crd,:) = X(rrd,:); continue; end
      if length(rrd) > length(crd)
        Y(crd,:) = F.factors(i).U \(F.factors(i).L'*X(rrd,:));
      else
        Y(crd,:) = F.factors(i).U'*(F.factors(i).L \X(rrd,:));
      end
    end

    % downward sweep in column space
    for i = n:-1:1
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      if i < n
        Y(rd,:) = Y(rd,:) - F.factors(i).F*Y(sk,:);
        Y(rd,:) = F.factors(i).U\Y(rd,:);
      end
      Y(sk,:) = Y(sk,:) - F.factors(i).cT*Y(rd,:);
    end

  % conjugate transpose
  else

    % upward sweep in column space
    for i = 1:n
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      X(rd,:) = X(rd,:) - F.factors(i).cT'*X(sk,:);
      if i < n
        X(rd,:) = F.factors(i).U'\X(rd,:);
        X(sk,:) = X(sk,:) - F.factors(i).F'*X(rd,:);
      end
    end

    % transfer from column to row space
    Y = zeros(F.M,size(X,2));
    for i = 1:n
      rrd = F.factors(i).rrd;
      crd = F.factors(i).crd;
      if i < n, Y(rrd,:) = X(crd,:); continue; end
      if length(rrd) > length(crd)
        Y(rrd,:) = F.factors(i).L *(F.factors(i).U'\X(crd,:));
      else
        Y(rrd,:) = F.factors(i).L'\(F.factors(i).U *X(crd,:));
      end
    end

    % downward sweep in row space
    for i = n:-1:1
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      if i < n
        Y(rd,:) = Y(rd,:) - F.factors(i).E'*Y(sk,:);
        Y(rd(F.factors(i).p),:) = F.factors(i).L'\Y(rd,:);
      end
      Y(sk,:) = Y(sk,:) - F.factors(i).rT'*Y(rd,:);
    end
  end
end
