% CHKFASTSV  Check fast solve mode parameter.
%
%    FASTSV = CHKFASTSV(FASTSV) validates the fast solve mode parameter FASTSV
%    and converts it to lowercase.

function fastsv = chkfastsv(fastsv)
  fastsv = lower(fastsv);
  assert(strcmp(fastsv,'n') || strcmpi(fastsv,'r') || strcmp(fastsv,'c') || ...
         strcmp(fastsv,'b'),'mfs-fds:chkfastsv:invalidFastsv', ...
         ['Fast solve mode parameter must be one of ''N'', ''R'', ''C'', or' ...
          ' ''B''.'])
end
