% RSKELFM  Recursive skeletonization factorization for MFS.
%
%    This is a generalization of RSKELF in FLAM for the MFS setting, i.e.,
%    consistent linear least squares with moderate solution norm. Specifically,
%    we make two important modifications during the processing of each block:
%
%      - If a block is too rectangular -- namely, to the point that the entire
%        block row or column, including the diagonal, is low-rank, given the
%        assumed structure -- then it is explicitly rank-reduced first as a
%        preprocessing step.
%
%      - Then row/column skeletons are computed as usual but augmented (i.e.,
%        redundant indices are subselected) as necessary to form a square
%        redundant-redundant submatrix so that it can be eliminated within the
%        existing framework.
%
%    The combination of these two steps leads to an algorithm that works only
%    with "mostly square" blocks. The result is again basically a generalized
%    LDU decomposition A = L*D*U, where L and U can now be rectangular, as they
%    may include rank reductions, but D remains square. Least squares systems
%    A*X ~ B can then be solved approximately as X = U\D\L\B, which should
%    achieve a small residual provided that B lies in the range of A.
%
%    Typical complexity for [M,N] = SIZE(A) with M >= N without loss of
%    generality: O(M + N) in 1D and O(M + N^(3*(1 - 1/D))) in D dimensions.
%
%    F = RSKELFM(A,RX,CX,OCC,RANK_OR_TOL) produces a factorization F of the
%    matrix A acting on the row and column points RX and CX, respectively, using
%    tree occupancy parameter OCC and local precision parameter RANK_OR_TOL. See
%    HYPOCT and ID for details. Note that both row/column points are sorted
%    together in the same tree, so OCC should be set roughly twice the desired
%    leaf size. Since no proxy function is supplied, this simply performs a
%    naive compression of all off-diagonal blocks.
%
%    F = RSKELFM(A,RX,CX,OCC,RANK_OR_TOL,PXYFUN) accelerates the compression
%    using the proxy function PXYFUN to capture the far field. This is a
%    function of the form
%
%      [KPXY,NBR] = PXYFUN(RC,RX,CX,SLF,NBR,L,CTR)
%
%    that is called for every block, where
%
%      - KPXY: interaction matrix against artificial proxy points
%      - NBR:  block neighbor point indices (can be modified)
%      - RC:   flag to specify row or column compression ('R' or 'C')
%      - RX:   input row points
%      - CX:   input column points
%      - SLF:  block point indices
%      - L:    block node size
%      - CTR:  block node center
%
%    The relevant arguments will be passed in by the algorithm; the user is
%    responsible for handling them. See the examples for further details. If
%    PXYFUN is not provided or empty (default), then the code uses the naive
%    global compression scheme.
%
%    F = RSKELFM(A,RX,CX,OCC,RANK_OR_TOL,PXYFUN,OPTS) also passes various
%    options to the algorithm. Valid options include:
%
%      - LVLMAX: maximum tree depth (default: LVLMAX = Inf). See HYPOCT.
%
%      - EXT: set the root node extent to [EXT(D,1) EXT(D,2)] along dimension D.
%             If EXT is empty (default), then the root extent is calculated from
%             the data. See HYPOCT.
%
%      - RDPIV: pivoting strategy for redundant point subselection. LU pivoting
%               is used if RDPIV = 'L', QR pivoting if RDPIV = 'Q', and random
%               if RDPIV = 'R' (default: RDPIV = 'L').
%
%      - FASTSV: fast solve mode by skipping rank reduction pseudoinverses (and
%                generation of their factors). Skip none if FASTSV = 'N', skip
%                for row reductions if FASTSV = 'R', skip for column reductions
%                if FASTSV = 'C', and skip for both if FASTSV = 'B' (default:
%                FASTSV = 'N').
%
%      - VERB: display status info if VERB = 1 (default: VERB = 0). This prints
%              to screen a table tracking row/column compression statistics
%              through level. Special levels: 'T', tree sorting.
%
%    See also HYPOCT, ID, RSKELFM_MV, RSKELFM_SV.

function F = rskelfm(A,rx,cx,occ,rank_or_tol,pxyfun,opts)

  % set default parameters
  if nargin < 6, pxyfun = []; end
  if nargin < 7, opts = []; end
  if ~isfield(opts,'lvlmax'), opts.lvlmax = Inf; end
  if ~isfield(opts,'ext'), opts.ext = []; end
  if ~isfield(opts,'rdpiv'), opts.rdpiv = 'l'; end
  if ~isfield(opts,'fastsv'), opts.fastsv = 'n'; end
  if ~isfield(opts,'verb'), opts.verb = 0; end

  % check inputs
  opts.rdpiv = chkrdpiv(opts.rdpiv);

  % print header
  if opts.verb
    fprintf([repmat('-',1,69) '\n'])
    fprintf('%3s | %6s | %19s | %19s | %10s\n', ...
            'lvl','nblk','start/end npts','start/end npts/blk','time (s)')
    fprintf([repmat('-',1,69) '\n'])
  end

  % build tree
  M = size(rx,2);
  N = size(cx,2);
  ts = tic;
  t = hypoct([rx cx],occ,opts.lvlmax,opts.ext);  % bundle row/col points
  te = toc(ts);
  for i = 1:t.lvp(t.nlvl+1)
    xi = t.nodes(i).xi;
    idx = xi <= M;
    t.nodes(i).rxi = xi( idx);      % extract row indices
    t.nodes(i).cxi = xi(~idx) - M;  % extract col indices
    t.nodes(i).xi = [];
  end
  if opts.verb, fprintf('%3s | %63.2e\n','t',te); end

  % count nonempty boxes at each level
  pblk = zeros(t.nlvl+1,1);
  for lvl = 1:t.nlvl
    pblk(lvl+1) = pblk(lvl);
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      if isempty([t.nodes(i).rxi t.nodes(i).cxi]), continue; end
      pblk(lvl+1) = pblk(lvl+1) + 1;
    end
  end

  % initialize
  nbox = t.lvp(end);
  e = cell(nbox,1);
  F = struct('pm',zeros(nbox,1),'psk',e,'prd',e,'pT',e,'pL',e,'rsk',e,'rrd', ...
             e,'csk',e,'crd',e,'rT',e,'cT',e,'L',e,'U',e,'p',e,'E',e,'F',e);
  F = struct('M',M,'N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F);
  nlvl = 0;
  n = 0;
  rrem = true(M,1); crem = true(N,1);  % which row/cols remain?
  S = cell(nbox,1);                    % storage for modified diagonal blocks
  rI = zeros(M,1); cI = zeros(N,1);    % for indexing

  % loop over tree levels
  for lvl = t.nlvl:-1:1
    ts = tic;
    nlvl = nlvl + 1;
    nrrem1 = nnz(rrem); ncrem1 = nnz(crem);  % remaining row/cols at start
    l = t.l(:,lvl);

    % pull up skeletons from children
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      t.nodes(i).rxi = [t.nodes(i).rxi [t.nodes(t.nodes(i).chld).rxi]];
      t.nodes(i).cxi = [t.nodes(i).cxi [t.nodes(t.nodes(i).chld).cxi]];
    end

    % loop over nodes
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      rslf = t.nodes(i).rxi;
      cslf = t.nodes(i).cxi;
      rnbr = [t.nodes(t.nodes(i).nbor).rxi];
      cnbr = [t.nodes(t.nodes(i).nbor).cxi];
      nrslf = length(rslf);
      ncslf = length(cslf);
      ctr = t.nodes(i).ctr;

      % generate modified diagonal block
      S{i} = full(A(rslf,cslf));
      if lvl < t.nlvl
        rI(rslf) = 1:nrslf;
        cI(cslf) = 1:ncslf;
        for j = t.nodes(i).chld
          ri = t.nodes(j).rxi;
          ci = t.nodes(j).cxi;
          S{i}(rI(ri),cI(ci)) = S{j};  % overwrite subblock from child
          S{j} = [];                   % clear child storage
        end
      end

      % preprocess if too rectangular
      pm = 'n'; psk = []; prd = []; pT = []; pL = [];
      rKpxy = []; cKpxy = [];
      if nrslf > ncslf
        [rKpxy,cnbr] = genpxy('r',pxyfun,rslf,cnbr,find(crem),rx,cx,lvl,l,ctr);
        if nrslf > ncslf + length(cnbr) + size(rKpxy,2)  % too many rows
          pm = 'r';
          K = [S{i} full(A(rslf,cnbr)) rKpxy]';
          [psk,prd,pT] = id(K,0);                        % square up
          if opts.fastsv ~= 'r' && opts.fastsv ~= 'b'
            pL = chol(eye(length(psk)) + pT*pT','lower');
          end
          rKpxy = rKpxy(psk,:);
          S{i} = S{i}(psk,:);
          psk = rslf(psk); prd = rslf(prd); rslf = psk;
          rrem(prd) = 0;
        end
      elseif ncslf > nrslf
        [cKpxy,rnbr] = genpxy('c',pxyfun,cslf,rnbr,find(rrem),rx,cx,lvl,l,ctr);
        if ncslf > nrslf + length(rnbr) + size(cKpxy,1);  % too many cols
          pm = 'c';
          K = [S{i}; full(A(rnbr,cslf)); cKpxy];
          [psk,prd,pT] = id(K,0);                         % square up
          if opts.fastsv ~= 'c' && opts.fastsv ~= 'b'
            pL = chol(eye(length(psk)) + pT*pT','lower');
          end
          cKpxy = cKpxy(:,psk);
          S{i} = S{i}(:,psk);
          psk = cslf(psk); prd = cslf(prd); cslf = psk;
          crem(prd) = 0;
        end
      end

      % compress row space
      if isempty(rKpxy)
        [rKpxy,cnbr] = genpxy('r',pxyfun,rslf,cnbr,find(crem),rx,cx,lvl,l,ctr);
      end
      K = [full(A(rslf,cnbr)) rKpxy]';
      [rsk,rrd,rT] = id(K,rank_or_tol);

      % compress column space
      if isempty(cKpxy)
        [cKpxy,rnbr] = genpxy('c',pxyfun,cslf,rnbr,find(rrem),rx,cx,lvl,l,ctr);
      end
      K = [full(A(rnbr,cslf)); cKpxy];
      [csk,crd,cT] = id(K,rank_or_tol);

      % find good redundant pivots
      K = S{i};
      nrrd = length(rrd);
      ncrd = length(crd);
      if nrrd > ncrd
        [rsk,rrd,rT] = rdpivot(K(rrd,crd) ,rsk,rrd,rT,opts.rdpiv);
      elseif nrrd < ncrd
        [csk,crd,cT] = rdpivot(K(rrd,crd)',csk,crd,cT,opts.rdpiv);
      end

      % move on if no compression
      if isempty(prd) && isempty(rrd) && isempty(crd), continue; end

      % compute factors
      K(rrd,:) = K(rrd,:) - rT'*K(rsk,:);
      K(:,crd) = K(:,crd) - K(:,csk)*cT;
      [L,U,p] = lu(K(rrd,crd),'vector');
      E = K(rsk,crd)/U;
      G = L\K(rrd(p),csk);
      S{i} = S{i}(rsk,csk) - E*G;  % update self-interaction

      % store matrix factors
      n = n + 1;
      F.factors(n).pm = pm;
      F.factors(n).psk = psk;
      F.factors(n).prd = prd;
      F.factors(n).pT = pT;
      F.factors(n).pL = pL;
      F.factors(n).rsk = rslf(rsk);
      F.factors(n).rrd = rslf(rrd);
      F.factors(n).csk = cslf(csk);
      F.factors(n).crd = cslf(crd);
      F.factors(n).rT = rT;
      F.factors(n).cT = cT;
      F.factors(n).L = L;
      F.factors(n).U = U;
      F.factors(n).p = p;
      F.factors(n).E = E;
      F.factors(n).F = G;

      % restrict to skeletons for next level
      t.nodes(i).rxi = rslf(rsk);
      t.nodes(i).cxi = cslf(csk);
      rrem(rslf(rrd)) = 0;
      crem(cslf(crd)) = 0;
    end
    F.lvp(nlvl+1) = n;
    te = toc(ts);

    % print summary
    if opts.verb
      nrrem2 = nnz(rrem); ncrem2 = nnz(crem);
      nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);  % nonempty up to this level
      fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e\n', ...
              lvl,nblk,nrrem1,nrrem2,nrrem1/nblk,nrrem2/nblk,te)
      fprintf('%3s | %6d | %8d | %8d | %8.2f | %8.2f | %10s\n', ...
              ' ',nblk,ncrem1,ncrem2,ncrem1/nblk,ncrem2/nblk,'')
    end
  end

  % finish
  F.factors = F.factors(1:n);
  if opts.verb, fprintf([repmat('-',1,69) '\n']), end
end

% wrapper for proxy interaction matrix generation
function [K,nbr] = genpxy(rc,pxyfun,slf,nbr,rem,rx,cx,lvl,l,ctr)
  nslf = length(slf);
  if rc ==  'r'
    K = zeros(nslf,0);
    if lvl > 2
      if isempty(pxyfun), nbr = setdiff(rem,slf);  % no proxy needed
      else, [K,nbr] = pxyfun('r',rx,cx,slf,nbr,l,ctr);
      end
    end
  else
    K = zeros(0,nslf);
    if lvl > 2
      if isempty(pxyfun), nbr = setdiff(rem,slf);  % no proxy needed
      else, [K,nbr] = pxyfun('c',rx,cx,slf,nbr,l,ctr);
      end
    end
  end
end
