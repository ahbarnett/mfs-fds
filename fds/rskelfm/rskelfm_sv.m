% RSKELFM_SV  Solve using recursive skeletonization factorization for MFS.
%
%    For a factorization F = L1*L2*...*D*...*U2*U1, where the L and U factors
%    are rectangular and D is square, this applies the "factored pseudoinverse"
%    U1\U2\...\INV(D)/.../L2/L1, provided that FASTSV = 'N' in RSKELFM. If
%    FASTSV ~= 'N', then a further approximation is used instead.
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
      if F.factors(i).pm == 'r' && ~isempty(F.factors(i).pL)
        psk = F.factors(i).psk;
        prd = F.factors(i).prd;
        tmp = X(psk,:) + F.factors(i).pT*X(prd,:);
        X(psk,:) = F.factors(i).pL'\(F.factors(i).pL\tmp);
      end
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      X(rd,:) = X(rd,:) - F.factors(i).rT'*X(sk,:);
      if i < n
        X(rd,:) = F.factors(i).L\X(rd(F.factors(i).p),:);
        X(sk,:) = X(sk,:) - F.factors(i).E*X(rd,:);
      end
    end

    % transfer from row to column space
    Y = zeros(F.N,size(X,2));
    Y([F.factors(1:n-1).crd],:) = X([F.factors(1:n-1).rrd],:);
    rrd = F.factors(n).rrd;
    crd = F.factors(n).crd;
    nrrd = length(rrd);
    ncrd = length(crd);
    if nrrd == ncrd     % LU
      Y(crd,:) = F.factors(n).U\(F.factors(n).L\X(rrd(F.factors(n).p)));
    elseif nrrd > ncrd  % QR
      Y(crd,:) = F.factors(n).U \(F.factors(n).L'*X(rrd,:));
    else                % LQ
      Y(crd,:) = F.factors(n).U'*(F.factors(n).L \X(rrd,:));
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
      if i < n
        X(rd,:) = F.factors(i).U'\X(rd,:);
        X(sk,:) = X(sk,:) - F.factors(i).F'*X(rd,:);
      end
    end

    % transfer from column to row space
    Y = zeros(F.M,size(X,2));
    Y([F.factors(1:n-1).rrd],:) = X([F.factors(1:n-1).crd],:);
    rrd = F.factors(n).rrd;
    crd = F.factors(n).crd;
    nrrd = length(rrd);
    ncrd = length(crd);
    if nrrd == ncrd     % LU
      Y(rrd(F.factors(n).p),:) = F.factors(n).L'\(F.factors(n).U'\X(crd,:));
    elseif nrrd > ncrd  % QR
      Y(rrd,:) = F.factors(i).L *(F.factors(i).U'\X(crd,:));
    else                % LQ
      Y(rrd,:) = F.factors(i).L'\(F.factors(i).U *X(crd,:));
    end

    % downward sweep in row space
    for i = n:-1:1
      sk = F.factors(i).rsk;
      rd = F.factors(i).rrd;
      if i < n
        Y(rd,:) = Y(rd,:) - F.factors(i).E'*Y(sk,:);
        Y(rd(F.factors(i).p),:) = F.factors(i).L'\Y(rd,:);
      end
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
