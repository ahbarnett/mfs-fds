% CHKRDPIV  Check redundant pivoting parameter.
%
%    RDPIV = CHKRDPIV(RDPIV) validates the redundant pivoting parameter RDPIV
%    and converts it to lowercase.

function rdpiv = chkrdpiv(rdpiv)
  rdpiv = lower(rdpiv);
  assert(strcmp(rdpiv,'l') || strcmpi(rdpiv,'q') || strcmp(rdpiv,'r'), ...
         'FLAM:chkrdpiv:invalidRdpiv', ...
         'Redundant pivoting parameter must be one of ''L'', ''Q'', or ''R''.')
end
