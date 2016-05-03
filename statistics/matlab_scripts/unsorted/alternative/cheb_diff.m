function [D] = cheb_diff(y)

% Computes the differentiation matrix D

n = length(y);

c   = [2 ; ones(n-2,1) ; 2].*(-1).^(1:n)';
y1  = repmat(y,1,n);
dy1 = y1 - y1';

D   = (c*(1./c)')./(dy1+(eye(n))); % off-diagonal entries
D   = D - diag(sum(D'));           % diagonal entries







