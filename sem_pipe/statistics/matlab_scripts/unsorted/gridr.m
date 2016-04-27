function [r, drdeta, d2rdeta, deta] = gridr(L, C, N)

rj = 0;
%
eta0 = asinh(-rj*sinh(C)/L)/C;
deta = (1+abs(eta0))/(N-1);
%
eta = eta0:deta:1;
%
drdeta = L*C*cosh(C*eta)/sinh(C);
d2rdeta = L*C^2*sinh(C*eta)/sinh(C);
%
r = rj + L*sinh(C*eta)/sinh(C);

%plot(r)
