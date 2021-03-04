function [Y,cres,niter] = ls(X)
global A nC tau N
n = size(X,2);
[Y,cres,niter] = lsedc(@lsfun,A(nC+1:end,:),zeros(N,n),A(1:nC,:),X,tau);
Y = Y(1:N,:);
end