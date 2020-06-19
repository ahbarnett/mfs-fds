% RSKELFM  Recursive skeletonization factorization for MFS.
%
%    This is a generalization of RSKELF in FLAM for the MFS setting, i.e.,
%    consistent linear least squares with moderate solution norm. Specifically,
%    we make two important modifications during the processing of each block:
%
%      - If an entire block row/column, including the diagonal, is detected as
%        low-rank, then it is explicitly rank-reduced first as a preprocessing
%        step.
%
%      - Row/column skeletons are then computed as usual but augmented (i.e.,
%        redundant indices are subselected) as necessary to form a square
%        redundant-redundant submatrix so that it can be eliminated within the
%        existing framework.
%
%    The combination of these two steps leads to an algorithm that works only
%    with "mostly square" blocks. The result is again basically a multilevel
%    generalized LDU decomposition A = L1*L2*...*D*...*U2*U1, where the L and U
%    factors can now be rectangular, as they may include rank reductions, but D
%    remains square. Least squares systems A*X ~ B can then be solved
%    approximately as X = U1\U2\...\D\...\L2\L1\B, which should achieve a small
%    residual provided that B lies in the range of A.
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
%      - LVLMAX: maximum tree depth (default: LVLMAX = INF). See HYPOCT.
%
%      - EXT: set the root node extent to [EXT(D,1) EXT(D,2)] along dimension D.
%             If EXT is empty (default), then the root extent is calculated from
%             the data. See HYPOCT.
%
%      - TMAX: ID interpolation matrix entry bound (default: TMAX = 2). See ID.
%
%      - RRQR_ITER: maximum number of RRQR refinement iterations in ID (default:
%                   RRQR_ITER = INF). See ID.
%
%      - RRATIO: aspect ratio for rank preprocessing. A given block is row-
%                compressed in its entirety if NRSLF > RRATIO*(NCSLF + NRSK),
%                where NRSLF and NCSLF are the number of row and column points,
%                respectively, in the block, and NRSK is the number of row
%                skeletons from compression of the corresponding off-diagonal
%                block row. An analogous criterion holds for column compression.
%                The default RRATIO = 0 corresponds to preprocessing all blocks
%                and is suitable for rank-deficient problems; otherwise,
%                RRATIO = 1 may be more appropriate.
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
  if ~isfield(opts,'Tmax'), opts.Tmax = 2; end
  if ~isfield(opts,'rrqr_iter'), opts.rrqr_iter = Inf; end
  if ~isfield(opts,'rratio'), opts.rratio = 0; end
  if ~isfield(opts,'rdpiv'), opts.rdpiv = 'l'; end
  if ~isfield(opts,'fastsv'), opts.fastsv = 'n'; end
  if ~isfield(opts,'verb'), opts.verb = 0; end

  % check inputs
  assert(opts.rratio >= 0, 'mfs-fds:rskelfm:invalidRratio', ...
         'Rank preprocessing ratio must be nonnegative.')
  opts.rdpiv = chkrdpiv(opts.rdpiv);
  opts.fastsv = chkfastsv(opts.fastsv);

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
  F = struct('rsk',e,'rrd',e,'rT',e,'csk',e,'crd',e,'cT',e, ...
             'psk',e,'prd',e,'pT',e,'pL',e,'qsk',e,'qrd',e,'qT',e,'qL',e, ...
             'L',e,'U',e,'p',e,'E',e,'F',e);
  F = struct('M',M,'N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F);
  nlvl = 0;
  n = 0;
  rrem = true(M,1); crem = true(N,1);  % which row/cols remain?
  S = cell(nbox,1);                    % storage for modified diagonal blocks
  rI = zeros(M,1); cI = zeros(N,1);    % for indexing
  tol = rem(rank_or_tol,1);            % relative tolerance for off-diagonal

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

      % compress off-diagonal block row
      Kpxy = zeros(nrslf,0);
      if lvl > 2
        if isempty(pxyfun), cnbr = setdiff(find(crem),cslf);
        else, [Kpxy,cnbr] = pxyfun('r',rx,cx,rslf,cnbr,l,t.nodes(i).ctr);
        end
      end
      rK = [full(A(rslf,cnbr)) Kpxy];
      [rsk,rrd,rT] = id(rK',rank_or_tol,opts.Tmax,opts.rrqr_iter);

      % compress off-diagonal block column
      Kpxy = zeros(0,ncslf);
      if lvl > 2
        if isempty(pxyfun), rnbr = setdiff(find(rrem),rslf);
        else, [Kpxy,rnbr] = pxyfun('c',rx,cx,cslf,rnbr,l,t.nodes(i).ctr);
        end
      end
      cK = [full(A(rnbr,cslf)); Kpxy];
      [csk,crd,cT] = id(cK,rank_or_tol,opts.Tmax,opts.rrqr_iter);

      % compress entire block row
      psk = []; prd = []; pT = []; pL = [];
      nskmax = ncslf + length(rsk);  % maximum rank
      if nrslf > opts.rratio*nskmax
        rK = [S{i} rK];
        ptol = nskmax + tol;         % up to tolerance or maximum rank
        [psk,prd,pT] = id(rK',ptol,opts.Tmax,opts.rrqr_iter,rsk);
        if opts.fastsv ~= 'r' && opts.fastsv ~= 'b'
          pL = chol(eye(length(psk)) + pT*pT','lower');
        end
        idx = ismember(rrd,psk); rrd = rrd(idx); rT = rT(:,idx);
        psk = rslf(psk); prd = rslf(prd);
        rrem(prd) = false;
      end

      % compress entire block column
      qsk = []; qrd = []; qT = []; qL = [];
      nskmax = nrslf + length(csk);  % maximum rank
      if ncslf > opts.rratio*nskmax
        cK = [S{i}; cK];
        qtol = nskmax + tol;         % up to tolerance or maximum rank
        [qsk,qrd,qT] = id(cK,qtol,opts.Tmax,opts.rrqr_iter,csk);
        if opts.fastsv ~= 'c' && opts.fastsv ~= 'b'
          qL = chol(eye(length(qsk)) + qT*qT','lower');
        end
        idx = ismember(crd,qsk); crd = crd(idx); cT = cT(:,idx);
        qsk = cslf(qsk); qrd = cslf(qrd);
        crem(qrd) = false;
      end

      % find good redundant pivots
      K = S{i};
      nrrd = length(rrd);
      ncrd = length(crd);
      if lvl > 1
        if nrrd > ncrd
          [rsk,rrd,rT] = rdpivot(K(rrd,crd) ,rsk,rrd,rT,opts.rdpiv);
          nrrd = length(rrd);
        elseif nrrd < ncrd
          [csk,crd,cT] = rdpivot(K(rrd,crd)',csk,crd,cT,opts.rdpiv);
          ncrd = length(crd);
        end
      end

      % move on if no compression
      if isempty(rrd) && isempty(crd) && isempty(prd) && isempty(qrd)
        continue
      end

      % compute factors
      rrs = [rrd rsk];
      crs = [crd csk];
      K(rrd,crs) = K(rrd,crs) - rT'*K(rsk,crs);
      K(rrs,crd) = K(rrs,crd) - K(rrs,csk)*cT;
      if nrrd == ncrd  % for all non-root
        [L,U,p] = lu(K(rrd,crd),'vector');
        E = K(rsk,crd)/U;
        G = L\K(rrd(p),csk);
        S{i} = S{i}(rsk,csk) - E*G;  % update self-interaction
      else             % can only happen at root
        if nrrd > ncrd, [L,U] = qr(K(rrd,crd) ,0);
        else,           [Q,R] = qr(K(rrd,crd)',0); L = R'; U = Q';
       end
        p = [];
        E = zeros(0,nrrd);
        G = zeros(ncrd,0);
      end

      % store matrix factors
      n = n + 1;
      F.factors(n).rsk = rslf(rsk);
      F.factors(n).rrd = rslf(rrd);
      F.factors(n).rT = rT;
      F.factors(n).csk = cslf(csk);
      F.factors(n).crd = cslf(crd);
      F.factors(n).cT = cT;
      F.factors(n).psk = psk;
      F.factors(n).prd = prd;
      F.factors(n).pT = pT;
      F.factors(n).pL = pL;
      F.factors(n).qsk = qsk;
      F.factors(n).qrd = qrd;
      F.factors(n).qT = qT;
      F.factors(n).qL = qL;
      F.factors(n).L = L;
      F.factors(n).U = U;
      F.factors(n).p = p;
      F.factors(n).E = E;
      F.factors(n).F = G;

      % restrict to skeletons for next level
      t.nodes(i).rxi = rslf(rsk);
      t.nodes(i).cxi = cslf(csk);
      rrem(rslf(rrd)) = false;
      crem(cslf(crd)) = false;
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
