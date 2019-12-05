% RDPIVOT  Pivoting for redundant point subselection.
%
%    For rectangular matrices, a given block will generally have different
%    numbers of redundant rows and columns (in the sense of the ID). This is a
%    problem for methods based on skeletonization factorization, which require
%    an invertible redundant-redundant submatrix to eliminate. One way to
%    reconcile these is to subselect redundant indices from the larger of the
%    row/column sets until the dimensions square up; this is equivalent to
%    augmenting the corresponding skeletons. We refer to the manner in which
%    this subselection is done as "pivoting" since we further want the resultant
%    redundant-redundant submatrix to be well-conditioned.
%
%    [SK,RD,T] = RDPIVOT(A,SK,RD,T,MODE) takes as input the current redundant-
%    redundant submatrix A with corresponding row ID skeleton and redundant
%    indices SK and RD, respectively (i.e., SIZE(A,1) = LENGTH(RD)), and
%    interpolation matrix T, and finds a new subset of redundant indices such
%    that the new redundant-redundant submatrix is square; the other ID outputs
%    are updated accordingly. The pivoting method is specified by the parameter
%    MODE:
%
%        MODE = 'L': LU pivoting
%        MODE = 'Q': QR pivoting (requires transpose)
%        MODE = 'R': random pivoting
%
%    Note that no error checking is done for the input MODE; this is delegated
%    to higher-level callers. This function is deliberately written in a row-
%    oriented manner in order to make LU pivoting (the tentative default) fast.
%
%    See also ID.

function [sk,rd,T] = rdpivot(A,sk,rd,T,mode)
  [m,n] = size(A);
  if m <= n || m ~= length(rd)
    warning('FLAM:rdpivot:badRdDim','Unexpected redundant submatrix size.');
    return  % warn and exit if short-and-fat or dimensions don't match
  end
  if n == 0  % special case (Octave LU misbehaves)
    sk = [sk rd];
    rd = [];
    T = zeros(length(sk),0);
    return
  end
  if     mode == 'l', [~,~,p] = lu(A ,'vector');  % find "best" redundant ...
  elseif mode == 'q', [~,~,p] = qr(A','vector');  % ... points to eliminate
  elseif mode == 'r', p = randperm(m);
  end
  sk = [sk rd(p(n+1:end))];  % augment skeletons with remainder
  idx = p(1:n);
  rd = rd(idx);
  T = [T(:,idx); zeros(length(sk)-size(T,1),n)];
end
