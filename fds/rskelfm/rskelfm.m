% RSKELFM  Recursive skeletonization factorization for MFS.
%
%    This is a generalization of RSKELF in FLAM for the MFS setting, i.e.,
%    consistent linear least squares. Specifically, row/column ID skeletons are
%    computed in the same way, but only a square subset is eliminated at each
%    step in order to fit into the existing framework. The resulting generalized
%    LDU decomposition now has rectangular D, where all blocks except the last
%    one are square and the final block has the same "excess rectangularity" as
%    the original matrix. Least squares systems A*X ~ B can then be solved
%    approximately as X = U\D\L\B, which should achieve small residuals provided
%    that B lies in RANGE(A).
%
%    Typical complexity for [M,N] = SIZE(A) with M >= N without loss of
%    generality: O(N + (M - N)*LOG^3(N)) in 1D and
%    O(N^(3*(1 - 1/D)) + (M - N)*N^(2*(1 - 1/D))) in D dimensions.
%
%    F = RSKELFM(A,RX,CX,OCC,RANK_OR_TOL) produces a factorization F of the
%    matrix A acting on the row and column points RX and CX, respectively, using
%    tree occupancy parameter OCC and local precision parameter RANK_OR_TOL. See
%    HYPOCT and ID for details. Since no proxy function is supplied, this simply
%    performs a naive compression of all off-diagonal blocks.
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
%      - VERB: display status info if VERB = 1 (default: VERB = 0). This prints
%              to screen a table tracking row/column compression statistics
%              through level. Two sets of data are shown for each type in terms
%              of the actual number of points remaining and the number of
%              compression "candidates" remaining, which constitute a subset of
%              the former and should scale as in the standard RSKELF. Special
%              levels: 'T', tree sorting.
%
%    See also HYPOCT, ID, RSKELFM_MV, RSKELFM_SV.

function F = rskelfm(A,rx,cx,occ,rank_or_tol,pxyfun,opts)

  % set default parameters
  if nargin < 6, pxyfun = []; end
  if nargin < 7, opts = []; end
  if ~isfield(opts,'lvlmax'), opts.lvlmax = Inf; end
  if ~isfield(opts,'ext'), opts.ext = []; end
  if ~isfield(opts,'rdpiv'), opts.rdpiv = 'l'; end
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
  F = struct('rsk',e,'rrd',e,'csk',e,'crd',e,'rT',e,'cT',e,'L',e,'U',e,'p', ...
             e,'E',e,'F',e);
  F = struct('M',M,'N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F);
  nlvl = 0;
  n = 0;
  rrem  = true(M,1); crem  = true(N,1);  % which row/cols remain?
  rcrem = true(M,1); ccrem = true(N,1);  % which row/cols remain for ID?
  S = cell(nbox,1);                      % storage for modified diagonal blocks
  rI = zeros(M,1); cI = zeros(N,1);      % auxiliary arrays for indexing

  % loop over tree levels
  for lvl = t.nlvl:-1:1
    ts = tic;
    nlvl = nlvl + 1;
    nrrem1  = sum(rrem ); ncrem1  = sum(crem );  % remaining row/cols at start
    nrcrem1 = sum(rcrem); nccrem1 = sum(ccrem);  % remaining candidates at start
    l = t.lrt/2^(lvl - 1);

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

      % find restriction to compression candidates
      rcslf = find(rcrem(rslf))';
      ccslf = find(ccrem(cslf))';
      rnbr = rnbr(find(rcrem(rnbr)));
      cnbr = cnbr(find(ccrem(cnbr)));

      % compress row space
      Kpxy = zeros(nrslf,0);
      if lvl > 2
        if isempty(pxyfun), cnbr = setdiff(find(crem),cslf);
        else, [Kpxy,cnbr] = pxyfun('r',rx,cx,rslf,cnbr,l,t.nodes(i).ctr);
        end
      end
      K = [full(A(rslf,cnbr)) Kpxy]';          % full matrix for block
      [rsk,~,~] = id(K(:,rcslf),rank_or_tol);  % compress among candidates only
      rsk = rcslf(rsk);
      rrd = setdiff(1:nrslf,rsk);
      rT = K(:,rsk)\K(:,rrd);                  % full interpolation operator

      % compress column space
      Kpxy = zeros(0,ncslf);
      if lvl > 2
        if isempty(pxyfun), rnbr = setdiff(find(rrem),rslf);
        else, [Kpxy,rnbr] = pxyfun('c',rx,cx,cslf,rnbr,l,t.nodes(i).ctr);
        end
      end
      K = [full(A(rnbr,cslf)); Kpxy];          % full matrix for block
      [csk,~,~] = id(K(:,ccslf),rank_or_tol);  % compress among candidates only
      csk = ccslf(csk);
      crd = setdiff(1:ncslf,csk);
      cT = K(:,csk)\K(:,crd);                  % full interpolation operator

      % update candidates
      rcrem(rslf(rrd)) = 0;
      ccrem(cslf(crd)) = 0;

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
      if isempty(rrd) && isempty(crd), continue; end

      % compute factors
      K(rrd,:) = K(rrd,:) - rT'*K(rsk,:);
      K(:,crd) = K(:,crd) - K(:,csk)*cT;
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
      F.factors(n).csk = cslf(csk);
      F.factors(n).crd = cslf(crd);
      F.factors(n).rT = rT';
      F.factors(n).cT = cT ;
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
      nrrem2  = sum(rrem ); ncrem2  = sum(crem );  % remaining row/cols at end
      nrcrem2 = sum(rcrem); nccrem2 = sum(ccrem);  % remaining candidates at end
      nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);  % nonempty up to this level
      fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e\n', ...
              lvl,nblk,nrrem1,nrrem2,nrrem1/nblk,nrrem2/nblk,te)
      fprintf('%3s | %6s | %8d | %8d | %8.2f | %8.2f | %10s\n', ...
              ' ',' ',nrcrem1,nrcrem2,nrcrem1/nblk,nrcrem2/nblk,'')
      fprintf('%3s | %6d | %8d | %8d | %8.2f | %8.2f | %10s\n', ...
              ' ',nblk,ncrem1,ncrem2,ncrem1/nblk,ncrem2/nblk,'')
      fprintf('%3s | %6s | %8d | %8d | %8.2f | %8.2f | %10s\n', ...
              ' ',' ',nccrem1,nccrem2,nccrem1/nblk,nccrem2/nblk,'')
    end
  end

  % finish
  F.factors = F.factors(1:n);
  if opts.verb, fprintf([repmat('-',1,69) '\n']), end
end
