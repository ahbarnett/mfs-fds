function [Kpxy,nbr] = pxyfun(rc,rx,cx,slf,nbr,l,ctr)
global proxy N
pxy = bsxfun(@plus,proxy*l,ctr');

%
%
% if strcmpi(rc,'r')
%     %      Kpxy = Kfun(rx(:,slf),pxy,'s')*(2*pi/N); % monopoles for proxy
%     % note: (2*pi/N) is a quadrature weight for unit circle geometry
%     Kpxy = Kfun(rx(:,slf),pxy)*(2*pi/N);
%     dx = cx(1,nbr) - ctr(1);
%     dy = cx(2,nbr) - ctr(2);
% elseif strcmpi(rc,'c')
%     %      Kpxy = Kfun(pxy,cx(:,slf),'d')*(2*pi/N); % why flips to dipoles here?
%     Kpxy = Kfun(pxy,cx(:,slf))*(2*pi/N);
%     dx = rx(1,nbr) - ctr(1);            % (makes no difference)
%     dy = rx(2,nbr) - ctr(2);
% end




if strcmpi(rc,'r')
%      Kpxy = Kfun(rx(:,slf),pxy,'s')*(2*pi/N);
    Kpxy = Kfun(rx(:,slf),pxy,'s');  % no need to scale for straight kernel interactions
    dx = cx(1,nbr) - ctr(1);
    dy = cx(2,nbr) - ctr(2);
    dz = cx(3,nbr) - ctr(3);
elseif strcmpi(rc,'c')
%      Kpxy = Kfun(pxy, cx(:,slf),'s')*(2*pi/N);
    Kpxy = Kfun(pxy, cx(:,slf),'s');

    dx = rx(1,nbr) - ctr(1);
    dy = rx(2,nbr) - ctr(2);
    dz = rx(3,nbr) - ctr(3);
end
dist = sqrt(dx.^2 + dy.^2 + dz.^2);
nbr = nbr(dist/l < 1.5);
end