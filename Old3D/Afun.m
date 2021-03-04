function A = Afun(i,j)
global rx cx 
% if isempty(i) || isempty(j)
%     A = zeros(length(i),length(j));
%     return
% end
% [I,J] = ndgrid(i,j);
% A = bsxfun(@times,Kfun(rx(:,i),cx(:,j),'s',nu(:,j)),area(j));
% M = spget(i,j);
% idx = M ~= 0;
% A(idx) = M(idx);
% A(I == J) = -0.5;


A = Kfun(rx(:,i),cx(:,j), 's');


end