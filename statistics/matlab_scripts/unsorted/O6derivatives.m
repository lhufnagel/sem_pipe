function [D1, D2] = O6derivatives(dh, d2h, deta, n)

%     Lele's Pade scheme 3/4/6 (6th order in interior)

M = zeros(n,n);
N = zeros(n,n);

a = 3.;
b = (2. + 4.*a)/3.;
c = (4. - a)/3.;
ca1 = 16.*(2.*a + 1.0)/(40. - a);
cb1 = (4.0*ca1 + 2.0)/3.0;
cc1 = (4.0 - ca1)/3.0;

a2 = 5.5;
b2 = (4.*a2 - 4.)/3.;
c2 = (10. - a2)/3.;
ca2 = a2;
cb2 = b2;
cc2 = c2;

for i = 4:n-3
    
    M(i, i-1) = 1;
    M(i, i) = a;
    M(i, i+1) = 1;
    
    N(i, i-2) = -c/(4*deta);
    N(i, i-1) = -b/(2*deta);
    N(i, i+1) = b/(2*deta);
    N(i, i+2) = c/(4*deta);
    
    M2(i, i+1) = 1;
    M2(i, i) = a2;
    M2(i, i-1) = 1;
    
    N2(i, i-2) = c2/(4*deta^2);
    N2(i, i-1) = b2/(deta^2);
    N2(i, i) = -2*((b2/(deta^2))+(c2/(4*deta^2)));
    N2(i, i+1) = b2/(deta^2);
    N2(i, i+2) = c2/(4*deta^2);
end

%Boundary Conditions

M(1, 1) = 2.0;
M(1, 2) = 4.0;  
M(2, 1) = 1.0;  
M(2, 2) = 4.0;
M(2, 3) = 1.0;
M(3, 2) = 1.0;
M(3, 3) = ca1;
M(3, 4) = 1.0;
M(n-2, n-3) = 1.0;
M(n-2, n-2) = ca1;
M(n-2, n-1) = 1.0;
M(n-1, n-2) = 1.0;
M(n-1, n-1) = 4.0;
M(n-1, n) = 1.0;
M(n, n-1) = 4.0;
M(n, n) = 2.0;

N(1, 1) = -5/deta;
N(1, 2) = 4/deta;
N(1, 3) = 1.0/deta;
N(2, 1) = -3/deta;
N(2, 3) = 3/deta;
N(3, 1) = -cc1/(4*deta);
N(3, 2) = -cb1/(2*deta);
N(3, 4) = cb1/(2*deta);
N(3, 5) = cc1/(4*deta);
N(n-2, n-4) = -cc1/(4*deta);
N(n-2, n-3) = -cb1/(2*deta);
N(n-2, n-1) = cb1/(2*deta);
N(n-2, n) = cc1/(4*deta);
N(n-1, n-2) = -3/deta;
N(n-1, n) = 3/deta;
N(n, n-2) = -1.0/deta;
N(n, n-1) = -4/deta;
N(n, n) = 5/deta;


M2(1, 1) = 1.0;
M2(1, 2) = 11.0;
M2(2, 1) = 1.0;
M2(2, 2) = 10.0;
M2(2, 3) = 1.0;
M2(3, 2) = 1.0;
M2(3, 3) = ca2;
M2(3, 4) = 1.0;
M2(n-2, n-3) = 1.0;
M2(n-2, n-2) = ca2;
M2(n-2, n-1) = 1.0;
M2(n-1, n-2) = 1.0;
M2(n-1, n-1) = 10.0;
M2(n-1, n) = 1.0;
M2(n, n-1) = 11.0;
M2(n, n) = 1.0;

N2(1, 1) = 13/(deta^2);
N2(1, 2) = -27/(deta^2);
N2(1, 3) = 15/(deta^2);
N2(1, 4) = -1/(deta^2);
N2(2, 1) = 12/(deta^2);
N2(2, 2) = -24/(deta^2);
N2(2, 3) = 12/(deta^2);
N2(3, 1) = cc2/(4*deta^2);
N2(3, 2) = cb2/(deta^2);
N2(3, 3) = (-2*cb2/(deta^2))-(2*cc2/(4*deta^2));
N2(3, 4) = cb2/(deta^2);
N2(3, 5) = cc2/(4*deta^2);
N2(n-2, n-4) = cc2/(4*deta^2);
N2(n-2, n-3) = cb2/(deta^2);
N2(n-2, n-2) = (-2*cb2/(deta^2))-(2*cc2/(4*deta^2));
N2(n-2, n-1) = cb2/(deta^2);
N2(n-2, n) = cc2/(4*deta^2);
N2(n-1, n-2) = 12/(deta^2);
N2(n-1, n-1) = -24/(deta^2);
N2(n-1, n) = 12/(deta^2);
N2(n, n-3) = -1/(deta^2);
N2(n, n-2) = 15/(deta^2);
N2(n, n-1) = -27/(deta^2);
N2(n, n) = 13/(deta^2);

d1 = inv(M)*N;
d2 = inv(M2)*N2;

for i = 1:n
    for j = 1:n
        D1(i,j) = d1(i,j)/dh(i);
        D2(i,j) = (d2(i,j)/dh(i)^2) - (d2h(i)*d1(i,j)/dh(i)^3);
    end
end