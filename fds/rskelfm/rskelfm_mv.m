% RSKELFM_MV  Multiply using recursive skeletonization factorization for MFS.
%
%    Typical complexity: O(N + (M - N)*LOG^2(N)) in 1D and
%    O(N^(2*(1 - 1/D)) + (M - N)*N^(1 - 1/D)) in D dimensions.
%
%    Y = RSKELFM_MV(F,X) produces the matrix Y by applying the factored matrix F
%    to the matrix X.
%
%    Y = RSKELFM_MV(F,X,TRANS) computes Y = F*X if TRANS = 'N' (default),
%    Y = F.'*X if TRANS = 'T', and Y = F'*X if TRANS = 'C'.
%
%    See also RSKELFM, RSKELFM_SV.

function Y = rskelfm_mv(F,X,trans)

  % set default parameters
  if nargin < 3 || isempty(trans), trans = 'n'; end

  % check inputs
  trans = chktrans(trans);

  % handle transpose by conjugation
  if trans == 't', Y = conj(rskelfm_mv(F,conj(X),'c')); return; end

  % initialize
  n = F.lvp(end);

  % no transpose
  if trans == 'n'

    % upward sweep in column space
    for i = 1:n
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      X(sk,:) = X(sk,:) + F.factors(i).cT*X(rd,:);
      if i < n
        X(rd,:) = F.factors(i).U*X(rd,:);
        X(rd,:) = X(rd,:) + F.factors(i).F*X(sk,:);
      end
    end

    % transfer from column to row space
    Y = zeros(F.M,size(X,2));
    for i = 1:n
      rrd = F.factors(i).rrd;
      crd = F.factors(i).crd;
      if i < n, Y(rrd,:) = X(crd,:);
      else,     Y(rrd,:) = F.factors(i).L*(F.factors(i).U*X(crd,:));
      end
    end

    % downward sweep in row space
    for i = n:-1:1
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      if i < n
        Y(sk,:) = Y(sk,:) + F.factors(i).E*Y(rd,:);
        Y(rd(F.factors(i).p),:) = F.factors(i).L*Y(rd,:);
      end
      Y(rd,:) = Y(rd,:) + F.factors(i).rT*Y(sk,:);
    end

  % conjugate transpose
  else

    % upward sweep in row space
    for i = 1:n
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      X(sk,:) = X(sk,:) + F.factors(i).rT'*X(rd,:);
      if i < n
        X(rd,:) = F.factors(i).L'*X(rd(F.factors(i).p),:);
        X(rd,:) = X(rd,:) + F.factors(i).E'*X(sk,:);
      end
    end

    % transfer from row to column space
    Y = zeros(F.N,size(X,2));
    for i = 1:n
      rrd = F.factors(i).rrd;
      crd = F.factors(i).crd;
      if i < n, Y(crd,:) = X(rrd,:);
      else,     Y(crd,:) = F.factors(i).U'*(F.factors(i).L'*X(rrd,:));
      end
    end

    % downward sweep in column space
    for i = n:-1:1
      sk = F.factors(i).csk;
      rd = F.factors(i).crd;
      if i < n
        Y(sk,:) = Y(sk,:) + F.factors(i).F'*Y(rd,:);
        Y(rd,:) = F.factors(i).U'*Y(rd,:);
      end
      Y(rd,:) = Y(rd,:) + F.factors(i).cT'*Y(sk,:);
    end
  end
end
