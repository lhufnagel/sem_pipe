
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Compute and plot the  statistics for %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   turbulent pipe flow simulation     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    G. K. El Khoury & A. Noorani      %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
%clear all
%clc

format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% The record in the binary file: stat000N is saved as follows %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  <U>                              % Q1
% 2.  <V>                              % Q2
% 3.  <W>                              % Q3
% 4.  <P>                              % Q4

% 5.  <uu> (u: instantaneous)          % Q5
% 6.  <vv> (v: instantaneous)          % Q6
% 7.  <ww> (w: instantaneous)          % Q7
% 8.  <pp> (p: instantaneous)          % Q8

% 9.  <uv> (u, v: instantaneous)       % Q9
% 10. <vw> (v, w: instantaneous)       % Q10
% 11. <uw> (u, w: instantaneous)       % Q11

% 12. <pu> (p, u: instantaneous)       % Q12
% 13. <pv> (p, v: instantaneous)       % Q13
% 14. <pw> (p, w: instantaneous)       % Q14

% 15. <pdudx> (p, dudx: instantaneous) % Q15
% 16. <pdudy> (p, dudy: instantaneous) % Q16
% 17. <pdudz> (p, dudz: instantaneous) % Q17

% 18. <pdvdx> (p, dvdx: instantaneous) % Q18
% 19. <pdvdy> (p, dvdy: instantaneous) % Q19
% 20. <pdvdz> (p, dvdz: instantaneous) % Q20

% 21. <pdwdx> (p, dwdx: instantaneous) % Q21
% 22. <pdwdy> (p, dwdy: instantaneous) % Q22
% 23. <pdwdz> (p, dwdz: instantaneous) % Q23

% 24. <uuu> (u: instantaneous)         % Q24
% 25. <vvv> (v: instantaneous)         % Q25
% 26. <www> (w: instantaneous)         % Q26
% 27. <ppp> (p: instantaneous)         % Q27

% 28. <uuv> (u, v: instantaneous)      % Q28
% 29. <uuw> (u, w: instantaneous)      % Q29
% 30. <vvu> (v, u: instantaneous)      % Q30
% 31. <vvw> (v, w: instantaneous) 	   % Q31   
% 32. <wwu> (w, u: instantaneous)      % Q32
% 33. <wwv> (w, v: instantaneous)      % Q33
% 34. <uvw> (u, v, w: instantaneous)   % Q34

% 35. <uuuu> (u, v: instantaneous)     % Q35
% 36. <vvvv> (u, w: instantaneous)     % Q36
% 37. <wwww> (v, u: instantaneous)     % Q37
% 38. <pppp> (v, w: instantaneous)     % Q38

% 39. <uuuv> (u: instantaneous)        % Q39
% 40. <uuvv> (v: instantaneous)        % Q40
% 41. <uvvv> (w: instantaneous)	       % Q41 

% 42. e11: <(du/dx.du/dx + du/dy.du/dy + du/dz.du/dz)> (u: instantaneous)    % Q42 
% 43. e22: <(dv/dx.dv/dx + dv/dy.dv/dy + dv/dz.dv/dz)> (v: instantaneous)    % Q43 
% 44. e33: <(dw/dx.dw/dx + dw/dy.dw/dy + dw/dz.dw/dz)> (w: instantaneous)    % Q44 
% 45. e12: <(du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz)> (u, v: instantaneous) % Q45  
% 46. e13: <(du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz)> (u, w: instantaneous) % Q46 
% 47. e23: <(dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz)> (v, w: instantaneous) % Q47 

% 48. <du/dx>       (u  : instantaneous) % Q48
% 49. <du/dy>       (u  : instantaneous) % Q49 
% 50. <du/dz>       (u  : instantaneous) % Q50 

% 51. <dv/dx>       (v  : instantaneous) % Q51
% 52. <dv/dy>       (v  : instantaneous) % Q52 
% 53. <dv/dz>       (v  : instantaneous) % Q53 

% 54. <dw/dx>       (w  : instantaneous) % Q54
% 55. <dw/dy>       (w  : instantaneous) % Q55 
% 56. <dw/dz>       (w  : instantaneous) % Q56 

% 57. <omegar>          (omegar  : instantaneous) % Q57
% 58. <omegat>          (omegat  : instantaneous) % Q58 
% 59. <omegaz>          (omegaZ  : instantaneous) % Q59 

% 60. <omegar*omegar>   (omegar  : instantaneous) % Q60
% 61. <omegat*omegat>   (omegat  : instantaneous) % Q61 
% 62. <omegaz*omegaz>   (omegaz  : instantaneous) % Q62 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Read the interpolated field %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fname')
  fname = '../recordings/polar_z_0.0';
end
[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
CH            = fread(fid,hdr,'*char')    ;
dum1          = fread(fid,1,'*float64')   ;
Rer           = fread(fid,1,'*float64')   ;
Domain        = fread(fid,3,'*float64')   ;
nel           = fread(fid,1,'int32')      ;
Poly          = fread(fid,3,'int32')      ;
nstat         = fread(fid,1,'int32')      ;
times         = fread(fid,1,'*float64')   ;
timee         = fread(fid,1,'*float64')   ;
atime         = fread(fid,1,'*float64')   ;
DT            = fread(fid,1,'*float64')   ;
nrec          = fread(fid,1,'int32')      ;
ng_r          = fread(fid,1,'int32')      ;
ng_theta      = fread(fid,1,'int32')      ;
CR            = fread(fid,1,'int32')      ;
dum2          = fread(fid,1,'*float64')   ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ri = 1;
% r1 = floor(ng_r/2    ); % index of middle grid - 1
% rm = floor(ng_r/2 + 1); % index of middle grid
% r2 = floor(ng_r/2 + 1); % index of middle grid + 1 
% re = ng_r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Genvarname is a huge piece of crap
for i=1:nstat
  clear(['Q' num2str(i)]);
end

for i =1:nstat

if i == 1
   fseek(fid,0,'cof');
else
   fseek(fid,8,'cof');
end
% dum4 = fread(fid,1   ,'int32')               ;
 MatR = fread(fid,[ng_r, ng_theta],'*float64'); 
 Mat  = flipdim(MatR,1)                       ;
 eval(['Q' num2str(i) '= Mat;'])
% dum4 = fread(fid,1   ,'int32')               ;
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reb = U(bulk)*D/nu 
% Ub: bulk velocity. 
% D : pipe diameter.
% nu: kinematic viscosity
D   = round(Domain(1));
Reb = Rer*D;
Ub  = 1;
R   = D/2;
rho = 1;
nu  = Ub*D/Reb; 
my  = nu*rho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Interpolated matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C(1,1)  corresponds to r = +1, and to theta = 0
% 
%          theta=0                              theta=2*pi 
% r = +1 [   *      *       *       *       *       *    ]
%        [   *      *       *       *       *       *    ]
%        [   *      *       *       *       *       *    ]
%        [   *      *       *       *       *       *    ]
% r = 0  [   *      *       *       *       *       *    ]  

% Later in the code the matrices are arranged in the following way:
%          theta=0                              theta=2*pi 
% r =  0 [   *      *       *       *       *       *    ]
%        [   *      *       *       *       *       *    ]
% r = +1 [   *      *       *       *       *       *    ]  

fname_r           = '../polar_mesh.bin';
[fid_r,message_r] = fopen(fname_r,'r','ieee-le');
hdr_r             = fread(fid_r,1,'int32');
dum3              = fread(fid_r,3,'int32');
dum4              = fread(fid_r,1,'*float64');
r                 = fread(fid_r,[ng_r, 1],'*float64');

fname_t           = '../theta_mesh.bin';
[fid_t,message_t] = fopen(fname_t,'r','ieee-le');
hdr_t             = fread(fid_t,1,'int32');
dum5              = fread(fid_t,3,'int32');
dum6              = fread(fid_t,1,'*float64');
theta             = fread(fid_t,[ng_theta, 1],'*float64');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r01 = (r(end:-1:1));  %%%%% grid in radial direction from 0 (center) to 1 (wall)
r10  = R - r01     ;  %%%%% grid in radial direction from 1 (center) to 0 (wall)
rr = length(r10)   ;  %%%%% 

atheta   = theta*180/pi;
n_theta  = ng_theta;  
n_r      = ng_r;

[th,r_] = meshgrid(theta,r01);
[X,Y] = pol2cart(th,r_);

ssp = rr:-1:1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% Use the 6th order modified Pade Scheme to calculate the derivatives
%%%%% The method is 3rd order at the boundaries
%%%%% Note that the derivatives matrices DM1 and DM2 need to be generated
%%%%% with the same constants L and C used to generate the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[rdum, dh, d2h, deta] = gridr(R, CR, n_r);
[DM1, DM2] = O6derivatives(dh, d2h, deta, n_r);

[theta_fr, Dfr] = fourdif(n_theta,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Mean velocities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  <U>  % Q1
% 2.  <V>  % Q2
% 3.  <W>  % Q3

RR1_car(1:3,1:n_r,1:n_theta) = 0;

for i = 1:n_r
   for j=1:n_theta
       RR1_car(:,i,j) = [Q1(i,j); Q2(i,j); Q3(i,j)];
    end
end

% squeeze(RR1_car(1,:,:)) = U(:,:)
% squeeze(RR1_car(2,:,:)) = V(:,:)
% squeeze(RR1_car(3,:,:)) = W(:,:) 

RR1 = cyl_Vel_R(RR1_car,theta,n_r,n_theta); 

% squeeze(RR1(1,:,:)) = Ur(:,:)
% squeeze(RR1(2,:,:)) = Ut(:,:)
% squeeze(RR1(3,:,:)) = Uz(:,:)

% R1(:,1) = Ur
% R1(:,2) = Ut    
% R1(:,3) = Uz

R1(1:n_r,1:3)     = 0;
R1_car(1:n_r,1:3) = 0;

for i=1:3
     R1(:,i)     = sum(squeeze(RR1(i,:,:)),2)./n_theta;   
     R1_car(:,i) = sum(squeeze(RR1_car(i,:,:)),2)./n_theta;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dUidxj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 48. <du/dx>       (u  : instantaneous) % Q48
% 49. <du/dy>       (u  : instantaneous) % Q49 
% 50. <du/dz>       (u  : instantaneous) % Q50 

% 51. <dv/dx>       (v  : instantaneous) % Q51
% 52. <dv/dy>       (v  : instantaneous) % Q52 
% 53. <dv/dz>       (v  : instantaneous) % Q53 

% 54. <dw/dx>       (w  : instantaneous) % Q54
% 55. <dw/dy>       (w  : instantaneous) % Q55 
% 56. <dw/dz>       (w  : instantaneous) % Q56 

dUidxj(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
         dUidxj(:,:,i,j) = [Q48(i,j) Q49(i,j) Q50(i,j);
                            Q51(i,j) Q52(i,j) Q53(i,j); 
                            Q54(i,j) Q55(i,j) Q56(i,j)];
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% Calculate the derivatives by the following method
%%%%%
%%%%% [ d(Ui)/dr ]   [  cos(theta)     sin(theta)    ][ d(Ui)/dx ]
%%%%% [          ]   [                               ][          ]
%%%%% [          ] = [                               ][          ]
%%%%% [          ]   [                               ][          ]
%%%%% [ d(Ui)/dt ]   [  -r*sin(theta)   r*cos(theta) ][ d(Ui)/dy ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
DRR1(1:n_r,1:n_theta,1:3) = 0;
DR1(1:n_r,1:3)            = 0;

for i = 1:3
for j = 1:n_theta
    
    DRR1(:,j,i) = squeeze(dUidxj(i,1,:,j)).*cos(theta(j)) + ...
                  squeeze(dUidxj(i,2,:,j)).*sin(theta(j));
        
end
end

for i = 1:3
   
    DR1(:,i) = sum(squeeze(DRR1(:,:,i)),2)/n_theta;
    
end
% % % dUdz1 = flipdim(-DM1*R1(ssp,3),1);
% % % dUdz  = DR1(:,3);
% % % DM3 = flipdim(DM1,2);
% % % dUdzt = DM3*R1(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Wall friction & Friction velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tau_wall = my*-DR1(rr,3);
u_tau    = sqrt(tau_wall/rho);
  
lstar = nu/u_tau;
tstar = nu/u_tau^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Present Reynolds number based on wall-friction velocity %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Re_tau = u_tau*R/nu;
Re_cl  = R1(1,3)*R/nu;
%DR1(:,3) = DR1(:,3).*(lstar/u_tau); % lstar/u_tau = my/tau_wall                                                                                                                                                             
                                        
r10_plus = r10*u_tau/nu;
r01_plus = r01*u_tau/nu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Velocity defect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Uzd   = R1(1,3) - R1(:,3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Shaper factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% The shape factor in boundary layer flows is defined as:
%%% H12 = δ1/δ2
%%%
%%% δ1: is the displacement thickness  
%%% δ1 = int(0,inf)[1 - u(y)/uo]dy 
%%% uo is the velocity in the free-stream outside the boundary layer
%%%
%%% δ2: is the momentum thickness 
%%% δ2 = int(0,inf)[u(y)/uo*{1 - u(y)/uo}]dy
%%%
%%% The shape factor in pipe flows is defined as:
%%% H12 = δ1/δ2
%%%
%%% δ1: is the displacement thickness  
%%% (δ1)(2R-δ1) = 2*int(0,R)[(r)*(1 - Uz(r)/Uz_cl)]dr 
%%% uo is the velocity in the free-stream outside the boundary layer
%%%
%%% δ2: is the momentum thickness 
%%% (δ2)(2R-δ2) = 2*int(0,R)[(r)*(Uz(r)/Uz_cl)*(1 - Uz(r)/Uz_cl)]dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cte_del1 = trapz(rdel,udel1)
%%% cte_del2 = trapz(rdel,udel2)

udel1 = 2*(r01).*(1 - R1(:,3)/R1(1,3));
udel2 = 2*(r01).*(1 - R1(:,3)/R1(1,3)).*(R1(:,3)/R1(1,3));
rdel  = r01;

Iudel1 = flipdim(-DM1\udel1(ssp),1);
Iudel2 = flipdim(-DM1\udel2(ssp),1);

cte_del1 = Iudel1(end) - Iudel1(1);
cte_del2 = Iudel2(end) - Iudel2(1);

a1 = 1; b1 = -2*R; c1 = cte_del1;
del1_1 = (-b1  - sqrt(b1^2 - 4*a1*c1))/(2*a1);
del1_2 = (-b1  + sqrt(b1^2 - 4*a1*c1))/(2*a1);

a2 = 1; b2 = -2*R; c2 = cte_del2;
del2_1 = (-b2  - sqrt(b2^2 - 4*a2*c2))/(2*a2);
del2_2 = (-b2  + sqrt(b2^2 - 4*a2*c2))/(2*a2);

H12 = del1_1/del2_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Reynolds stress tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  <U>                              % Q1
% 2.  <V>                              % Q2
% 3.  <W>                              % Q3

% 5.  <uu> (u: instantaneous)          % Q5
% 6.  <vv> (v: instantaneous)          % Q6
% 7.  <ww> (w: instantaneous)          % Q7
% 9.  <uv> (u, v: instantaneous)       % Q9
% 10. <vw> (v, w: instantaneous)       % Q10
% 11. <uw> (u, w: instantaneous)       % Q11

RR2_car(1:3,1:3,1:n_r,1:n_theta) = 0;
 
for j = 1:n_theta
   for i = 1:n_r
        RR2_car(:,:,i,j) = [(Q5(i,j) -Q1(i,j).*Q1(i,j)) (Q9(i,j) -Q1(i,j).*Q2(i,j))  (Q11(i,j) -Q1(i,j).*Q3(i,j));
                            (Q9(i,j) -Q1(i,j).*Q2(i,j)) (Q6(i,j) -Q2(i,j).*Q2(i,j))  (Q10(i,j) -Q2(i,j).*Q3(i,j)); 
                            (Q11(i,j)-Q1(i,j).*Q3(i,j)) (Q10(i,j)-Q2(i,j).*Q3(i,j))  (Q7(i,j)  -Q3(i,j).*Q3(i,j))];                        
    end
end    
   
RR2 = cyl_Rey_ten_R(RR2_car,theta,n_r,n_theta); 

R2(1:3,1:3,1:n_r)     = 0; 
R2_car(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      R2(:,:,i)     = R2(:,:,i)     + RR2(:,:,i,j)./n_theta;
      R2_car(:,:,i) = R2_car(:,:,i) + RR2_car(:,:,i,j)./n_theta;
    end
end

k=0.5*(RR2(1,1,:,:)+RR2(2,2,:,:)+RR2(3,3,:,:));
% mdf_k=squeeze(k)./(u_tau^2);
% >> mdf_k(:,end+1)=mdf_k(:,1);
% >> pcolor(XX,YY,mdf_k);
% >> shading flat;axis equal tight;colorbar;rita;axis off;
% >> shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Total shear stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tau_rt  = squeeze(R2(3,1,:))./tau_wall - DR1(:,3).*(lstar/u_tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Residual Rsd = tau_rt - r/R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Rsd  = r01/R - tau_rt;
SRsd = Rsd.^2;
IRsd = flipdim(-DM1\SRsd(ssp),1);

%%% The error in the discrete L2-norm.
error_L2  = sqrt(IRsd(end)-IRsd(1));

%%% dUdz     = flipdim(-DM1*R1(ssp,3),1) ;
%%% Uz       = flipdim(-DM1\dUdz1(ssp),1);
%%% error_L2_1 = sqrt(trapz(r01,Rsd.^2)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Mean, r.m.s, Skewness and Flatness of pressure %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4 . <P>                              % Q4
% 8 . <pp>   (p: instantaneous)        % Q8
% 27. <ppp>  (p: instantaneous)        % Q27
% 38. <pppp> (v, w: instantaneous)     % Q38

P_mean    = Q4;
pp_mean   = Q8;
ppp_mean  = Q27;
pppp_mean = Q38;

pfpf_mean = pp_mean - (P_mean).*(P_mean);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Skewness  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
pfpfpf_mean = ppp_mean - P_mean.*P_mean.*P_mean - 3*P_mean.*pfpf_mean;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flatness  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
pfpfpfpf_mean = pppp_mean -4*pfpfpf_mean.*P_mean - ...
   6*pfpf_mean.*P_mean.*P_mean - P_mean.*P_mean.*P_mean.*P_mean;

P    = sum(P_mean,2       )./n_theta;
pp   = sum(pfpf_mean,2    )./n_theta; 
ppp  = sum(pfpfpf_mean,2  )./n_theta; 
pppp = sum(pfpfpfpf_mean,2)./n_theta; 

prms = sqrt(pp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
skew_p = ppp./(pp).^(3/2);
flat_p = pppp./(pp).^(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Skewness and Flatness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 24. <uuu> (u: instantaneous)         % Q24
% 25. <vvv> (v: instantaneous)         % Q25
% 26. <www> (w: instantaneous)         % Q26
% 28. <uuv> (u, v: instantaneous)      % Q28
% 29. <uuw> (u, w: instantaneous)      % Q29
% 30. <vvu> (v, u: instantaneous)      % Q30
% 31. <vvw> (v, w: instantaneous) 	   % Q31   
% 32. <wwu> (w, u: instantaneous)      % Q32
% 33. <wwv> (w, v: instantaneous)      % Q33
% 34. <uvw> (u, v, w: instantaneous)   % Q34

% The tensor has the following form:
%
% [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ] 
% [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw wwv ]
% [ uuw uvw uww ] [ uvw vvw vww ] [ uww wwv www ]
%

tensRR3_car(1:3,1:3,1:3,1:n_r,1:n_theta)  = 0;
tensRR3f_car(1:3,1:3,1:3,1:n_r,1:n_theta) = 0;
tensRR3(1:3,1:3,1:3,1:n_r,1:n_theta)      = 0;
tensR3(1:3,1:3,1:3,1:n_r)                 = 0;  

for j=1:n_theta
    for i=1:n_r        

tensRR3_car(:,:,1,i,j) = [Q24(i,j) Q28(i,j) Q29(i,j); 
                          Q28(i,j) Q30(i,j) Q34(i,j); 
                          Q29(i,j) Q34(i,j) Q32(i,j)];
             
tensRR3_car(:,:,2,i,j) = [Q28(i,j) Q30(i,j) Q34(i,j); 
                          Q30(i,j) Q25(i,j) Q31(i,j); 
                          Q34(i,j) Q31(i,j) Q33(i,j)];
             
tensRR3_car(:,:,3,i,j) = [Q29(i,j) Q34(i,j) Q32(i,j); 
                          Q34(i,j) Q31(i,j) Q33(i,j); 
                          Q32(i,j) Q33(i,j) Q26(i,j)];   
        
    end
end

for it=1:3
 for jt=1:3
  for kt=1:3
      for j = 1:n_theta
        
   tensRR3f_car(it,jt,kt,:,j) = squeeze(tensRR3_car(it,jt,kt,:,j)) - ...
            squeeze(RR2_car(it,jt,:,j)).*squeeze(RR1_car(kt,:,j))' - ...
            squeeze(RR2_car(it,kt,:,j)).*squeeze(RR1_car(jt,:,j))' - ...
            squeeze(RR2_car(jt,kt,:,j)).*squeeze(RR1_car(it,:,j))' - ...
            squeeze(RR1_car(it,:,j))'.*squeeze(RR1_car(jt,:,j))'.*squeeze(RR1_car(kt,:,j))';
 
      end
   end
 end
end

for j = 1:n_theta
   
 Rot_m   = [ cos(theta(j)) sin(theta(j)) 0;
            -sin(theta(j)) cos(theta(j)) 0;
                  0           0          1];    
for m=1:3
for n = 1:3
for p = 1:3
%------------------------
   for rx = 1:3
   for sx = 1:3
   for tx = 1:3
    tensRR3(m,n,p,:,j) = tensRR3(m,n,p,:,j) + ...
               Rot_m(m,rx).*Rot_m(n,sx).*Rot_m(p,tx).*tensRR3f_car(rx,sx,tx,:,j);
   end
   end
   end
%-------------------------        
end
end
end 
            
clear Rot_m 
 
end

for i = 1:n_r
 for j = 1:n_theta 
   tensR3(:,:,:,i)  = tensR3(:,:,:,i) + tensRR3(:,:,:,i,j)./n_theta;
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Skewness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
skew_ur = squeeze(tensR3(1,1,1,:))./(squeeze(R2(1,1,:))).^(3/2);
skew_ut = squeeze(tensR3(2,2,2,:))./(squeeze(R2(2,2,:))).^(3/2);
skew_uz = squeeze(tensR3(3,3,3,:))./(squeeze(R2(3,3,:))).^(3/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Flatness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 35. <uuuu> (u, v: instantaneous)     % Q35
% 36. <vvvv> (u, w: instantaneous)     % Q36
% 37. <wwww> (v, u: instantaneous)     % Q37
% 38. <pppp> (v, w: instantaneous)     % Q38

urururur_inst(1:n_r,1:n_theta) = 0;
utututut_inst(1:n_r,1:n_theta) = 0;

for j = 1:n_theta 

urururur_inst(:,j) = Q35(:,j).*cos(theta(j)).^4 + ...
                   6*Q40(:,j).*cos(theta(j)).^2.*sin(theta(j)).^2 + ...
                   4*Q39(:,j).*cos(theta(j)).^3.*sin(theta(j))    + ...
                   4*Q41(:,j).*cos(theta(j)).*sin(theta(j)).^3    + ...
                     Q36(:,j).*sin(theta(j)).^4;
                 
utututut_inst(:,j) = Q35(:,j).*sin(theta(j)).^4 + ...
                   6*Q40(:,j).*cos(theta(j)).^2.*sin(theta(j)).^2 - ...
                   4*Q39(:,j).*cos(theta(j)).*sin(theta(j)).^3    - ...
                   4*Q41(:,j).*cos(theta(j)).^3.*sin(theta(j))    + ...
                     Q36(:,j).*cos(theta(j)).^4;            
                                    
end

urururur_2d = urururur_inst - 4*squeeze(tensRR3(1,1,1,:,:)).*squeeze(RR1(1,:,:)) - ...
     6*squeeze(RR2(1,1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)) - ... 
       squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)); 

utututut_2d = utututut_inst - 4*squeeze(tensRR3(2,2,2,:,:)).*squeeze(RR1(2,:,:)) - ...
              6*squeeze(RR2(2,2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)) - ...
        squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:));           
          
uzuzuzuz_2d = Q37 - 4*squeeze(tensRR3(3,3,3,:,:)).*squeeze(RR1(3,:,:)) - ...
      6*squeeze(RR2(3,3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)) - ...
        squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)); 
          
urururur = sum(urururur_2d,2)./n_theta;   
utututut = sum(utututut_2d,2)./n_theta;  
uzuzuzuz = sum(uzuzuzuz_2d,2)./n_theta;        

flat_ur = urururur./squeeze(R2(1,1,:)).^(2);
flat_ut = utututut./squeeze(R2(2,2,:)).^(2);
flat_uz = uzuzuzuz./squeeze(R2(3,3,:)).^(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Convection tensor Cij, as a function of r %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Crr = 0;
Ctt = 0;
Czz = 0;
Crt = 0;
Crz = 0;
Ctz = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Production tensor Pij, as a function of r %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PIJ = -(<uiuk>*dUi/dxk + <ujuk>*dUj/dxk)
% P11 = 2<uu>dU/dx + 2<uv>dU/dy 
% P12 = <uu>dV/dx  + <uv>dV/dy  + <uv>dU/dx + <vv>dU/dy
% P13 = <uu>dW/dx  + <uv>dW/dy  + <uw>dU/dx + <vw>dU/dy
% P22 = 2<uv>dV/dx + 2<vv>dV/dy 
% P23 = <uv>dW/dx  + <vv>dW/dy  + <uw>dV/dx + <vw>dV/dy
% P33 = 2<uw>dW/dx + 2<vw>dW/dy 

PR_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for it = 1:3
 for jt = 1:3

  for j = 1:n_theta
   for dummy = 1:3
         
    PR_car(it,jt,:,j) = squeeze(PR_car(it,jt,:,j)) - ( ...
                        squeeze(RR2_car(it,dummy,:,j)).*squeeze(dUidxj(jt,dummy,:,j)) + ...
                        squeeze(RR2_car(jt,dummy,:,j)).*squeeze(dUidxj(it,dummy,:,j)));
                 
   end
  end
        
 end
end

PR_cyl = cyl_Rey_ten_R(PR_car,theta,n_r,n_theta); 

PRij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      PRij(:,:,i) = PRij(:,:,i) + PR_cyl(:,:,i,j)./n_theta;
    end
end

Prr = zeros(n_r,1); 
Ptt = zeros(n_r,1);
Pzz = -2*squeeze(R2(3,1,:)).*(DR1(:,3)); 
Prt = zeros(n_r,1); 
Prz = -squeeze(R2(1,1,:)).*(DR1(:,3));
Ptz = zeros(n_r,1); 

% Prr = squeeze(PRij(1,1,:));
% Ptt = squeeze(PRij(2,2,:));
% Pzz = squeeze(PRij(3,3,:));
% Prt = squeeze(PRij(1,2,:));
% Prz = squeeze(PRij(1,3,:));
% Ptz = squeeze(PRij(2,3,:));
% ----------------------------------------------- 
% Prr = 0
% Ptt = 0
% Pzz = -2*<uzur>*dUz/dr
% Prt = 0
% Prz = -<urur>*dUz/dr
% Ptz = 0
% Prr1 = zeros(n_r,1); 
% Ptt1 = zeros(n_r,1);
% Pzz1 = -2*squeeze(R2(3,1,:)).*(DR1(:,3)); 
% Prt1 = zeros(n_r,1); 
% Prz1 = -squeeze(R2(1,1,:)).*(DR1(:,3));
% Ptz1 = zeros(n_r,1); 
% ----------------------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Turbulent diffusion tensor TDij, as a function of r %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Tij = -d<uiujuk>/dxk
% T11 = - (d<uuuk>)/dxk, K=1,3 = -[d<uuu>/dx + d<uuv>/dy + d<uuw>/dz]
% T22 = - (d<vvuk>)/dxk, K=1,3 = -[d<vvu>/dx + d<vvv>/dy + d<vvw>/dz]
% T33 = - (d<wwuk>)/dxk, K=1,3 = -[d<wwu>/dx + d<wwv>/dy + d<www>/dz]
% T12 = - (d<uvuk>)/dxk, K=1,3 = -[d<uvu>/dx + d<uvv>/dy + d<uvw>/dz]
% T13 = - (d<uwuk>)/dxk, K=1,3 = -[d<uwu>/dx + d<uwv>/dy + d<uww>/dz]
% T23 = - (d<vwuk>)/dxk, K=1,3 = -[d<vwu>/dx + d<vwv>/dy + d<vww>/dz]

% TDrr = -1/r*(d(r<ururur>)/dr) + 2/r*<urutut> = -[(1/r)*<ururur> + d(ururur)/dr] + (2/r)*<urutut> 
% TDtt = -1/r*(d(r<urutut>)/dr) - 2/r*<urutut> = -[(1/r)*<urutut> + d(urutut)/dr] - (2/r)*<urutut> 
% TDzz = -1/r*(d(r<uruzuz>)/dr)                = -[(1/r)*<uruzuz> + d(uruzuz)/dr]
% TDrt = 0
% TDrz = -1/r*(d(r<ururuz>)/dr) + 1/r*<ututuz> = -[(1/r)*<ururuz> + d(ururuz)/dr] + (1/r)*<ututuz>
% TDtz = 0

TD_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for it = 1:3
for jt = 1:3

   [dQdxj] = dQ3dxj_t(it,jt,n_r,n_theta,r01,theta,ssp,DM1,Dfr,tensRR3f_car); 
  
   for dummy = 1:3
                                               
    TD_car(it,jt,:,:) = squeeze(TD_car(it,jt,:,:)) - squeeze(dQdxj(dummy,:,:));    

  end
        
end
end

TD_cyl = cyl_Rey_ten_R(TD_car,theta,n_r,n_theta); 

TDij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      TDij(:,:,i) = TDij(:,:,i) + TD_cyl(:,:,i,j)./n_theta;
    end
end

TDrr = squeeze(TDij(1,1,:));
TDtt = squeeze(TDij(2,2,:));
TDzz = squeeze(TDij(3,3,:));
TDrt = squeeze(TDij(1,2,:));
TDrz = squeeze(TDij(1,3,:));
TDtz = squeeze(TDij(2,3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%% Calculated without tensor rotation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% TDrr_1 = (1./r01).*squeeze(tensR3(1,1,1,:));
%%%% TDrr_2 = flipdim(-DM1*squeeze(tensR3(1,1,1,ssp)),1);
%%%% TDrr_3 = (2./r01).*squeeze(tensR3(2,2,1,:));
%%%% 
%%%% TDtt_1 = (1./r01).*squeeze(tensR3(2,2,1,:));
%%%% TDtt_2 = flipdim(-DM1*squeeze(tensR3(2,2,1,ssp)),1);
%%%% TDtt_3 = (2./r01).*squeeze(tensR3(2,2,1,:));
%%%% 
%%%% 
%%%% TDzz_1 = (1./r01).*squeeze(tensR3(3,3,1,:)); 
%%%% TDzz_2 = flipdim(-DM1*squeeze(tensR3(3,3,1,ssp)),1);
 
%%%% TDrz_1 = (1./r01).*squeeze(tensR3(1,3,1,:));
%%%% TDrz_2 = flipdim(-DM1*squeeze(tensR3(1,3,1,ssp)),1);
%%%% TDrz_3 = (1./r01).*squeeze(tensR3(2,2,3,:));
 
%%%% TDrr = -(TDrr_1 + TDrr_2) + TDrr_3;
%%%% TDtt = -(TDtt_1 + TDtt_2) - TDtt_3;
%%%% TDzz = -(TDzz_1 + TDzz_2);
%%%% TDrt = zeros(n_r,1);
%%%% TDrz = -(TDrz_1 + TDrz_2) + TDrz_3;
%%%% TDtz = zeros(n_r,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Pressure diffusion tensor PDij, as a function of r %%%%%%%%%%%%%%%%
%%%%%%%%%%% Divergence of pressure-velocity correlation %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. <pu> (p, u: instantaneous)       % Q12
% 13. <pv> (p, v: instantaneous)       % Q13
% 14. <pw> (p, w: instantaneous)       % Q14

% Gij = -(1/rho)*(d<puj>/dxi + d<pui>/dxj)
% G11 = -(1/rho)*(d<pu>/dx + d<pu>/dx)
% G12 = -(1/rho)*(d<pv>/dx + d<pu>/dy) 
% G13 = -(1/rho)*(d<pw>/dx + d<pu>/dz) 
% G22 = -(1/rho)*(d<pv>/dy + d<pv>/dy)
% G23 = -(1/rho)*(d<pw>/dy + d<pv>/dz) 
% G33 = -(1/rho)*(d<pw>/dz + d<pw>/dz)

% <pu> = <pu> - PU;
% <pv> = <pv> - PV;
% <pw> = <pw> - PW;

PU_car(1:3,1:n_r,1:n_theta) = 0;

for i = 1:n_r
   for j=1:n_theta
       PU_car(:,i,j) = [Q12(i,j) - Q4(i,j).*Q1(i,j)  ;...
                        Q13(i,j) - Q4(i,j).*Q2(i,j) ;...
                        Q14(i,j) - Q4(i,j).*Q3(i,j)];
    end
end

[dpuidxj] = dUidxj_t(n_r,n_theta,r01,theta,ssp,DM1,Dfr,PU_car); 

G11 = (-2/rho)*(squeeze(dpuidxj(1,1,:,:)));
G22 = (-2/rho)*(squeeze(dpuidxj(2,2,:,:)));
G33 = (-2/rho)*(squeeze(dpuidxj(3,3,:,:)));
G12 = (-1/rho)*(squeeze(dpuidxj(2,1,:,:)) + squeeze(dpuidxj(1,2,:,:)));
G13 = (-1/rho)*(squeeze(dpuidxj(3,1,:,:)) + squeeze(dpuidxj(1,3,:,:)));
G23 = (-1/rho)*(squeeze(dpuidxj(3,2,:,:)) + squeeze(dpuidxj(2,3,:,:)));

PD_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
        PD_car(:,:,i,j) = [G11(i,j) G12(i,j) G13(i,j);
                           G12(i,j) G22(i,j) G23(i,j); 
                           G13(i,j) G23(i,j) G33(i,j)];
   end

end    
   
PD_cyl = cyl_Rey_ten_R(PD_car,theta,n_r,n_theta); 

PDij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      PDij(:,:,i) = PDij(:,:,i) + PD_cyl(:,:,i,j)./n_theta;
    end
end

PDrr = squeeze(PDij(1,1,:));
PDtt = squeeze(PDij(2,2,:));
PDzz = squeeze(PDij(3,3,:));
PDrt = squeeze(PDij(1,2,:));
PDrz = squeeze(PDij(1,3,:));
PDtz = squeeze(PDij(2,3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Pressure strain tensor PSij, as a function of r %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 15. <pdudx> (p, dudx: instantaneous) % Q15
% 16. <pdudy> (p, dudy: instantaneous) % Q16
% 17. <pdudz> (p, dudz: instantaneous) % Q17

% 18. <pdvdx> (p, dvdx: instantaneous) % Q18
% 19. <pdvdy> (p, dvdy: instantaneous) % Q19
% 20. <pdvdz> (p, dvdz: instantaneous) % Q20

% 21. <pdwdx> (p, dwdx: instantaneous) % Q21
% 22. <pdwdy> (p, dwdy: instantaneous) % Q22
% 23. <pdwdz> (p, dwdz: instantaneous) % Q23

% PIij = (1/rho)*(<pdui/dxj + pduj/dxi>)
% Pi11 = (1/rho)*(<pdu/dx> + <pdu/dx>)
% Pi12 = (1/rho)*(<pdu/dy> + <pdv/dx>)
% Pi13 = (1/rho)*(<pdu/dz> + <pdw/dx>)
% Pi22 = (1/rho)*(<pdv/dy> + <pdv/dy>)
% Pi23 = (1/rho)*(<pdv/dz> + <pdw/dy>)
% Pi33 = (1/rho)*(<pdw/dz> + <pdv/dy>)

% Pi11 = (2/rho)*(<pdu/dx> - <P>*dU/dx)
% Pi12 = (1/rho)*(<pdu/dy> + <pdv/dx>  - <P>dU/dy - <P>dV/dx)
% Pi13 = (1/rho)*(<pdu/dz> + <pdw/dx>  - <P>dU/dz - <P>dW/dx)
% Pi22 = (2/rho)*(<pdv/dy> - <P>*dV/dy)
% Pi23 = (1/rho)*(<pdv/dz> + <pdw/dy>  - <P>dV/dz - <P>dW/dy)
% Pi33 = (2/rho)*(<pdw/dz> - <P>dW/dz)

Pi11(1:n_r,1:n_theta) = 0; 
Pi22(1:n_r,1:n_theta) = 0;
Pi33(1:n_r,1:n_theta) = 0;
Pi12(1:n_r,1:n_theta) = 0;
Pi13(1:n_r,1:n_theta) = 0;
Pi23(1:n_r,1:n_theta) = 0;

for j = 1:n_theta
        
    Pi11(:,j) = (2/rho)*(Q15(:,j) - Q4(:,j).*(squeeze(dUidxj(1,1,:,j)))); 
    Pi22(:,j) = (2/rho)*(Q19(:,j) - Q4(:,j).*(squeeze(dUidxj(2,2,:,j)))); 
    Pi33(:,j) = (2/rho)*(Q23(:,j) - Q4(:,j).*(squeeze(dUidxj(3,3,:,j)))); 
    
    Pi12(:,j) = (1/rho)*(Q16(:,j) + Q18(:,j)   - ...
        Q4(:,j).*(squeeze(dUidxj(1,2,:,j)))    - ...
        Q4(:,j).*(squeeze(dUidxj(2,1,:,j))) ); 
    Pi13(:,j) = (1/rho)*(Q17(:,j) + Q21(:,j)   - ...
        Q4(:,j).*(squeeze(dUidxj(1,3,:,j)))    - ...
        Q4(:,j).*(squeeze(dUidxj(3,1,:,j))) ); 
    Pi23(:,j) = (1/rho)*(Q20(:,j) + Q22(:,j)   - ...
        Q4(:,j).*(squeeze(dUidxj(2,3,:,j)))    - ...
        Q4(:,j).*(squeeze(dUidxj(3,2,:,j))) ); 
    
end

PS_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
        PS_car(:,:,i,j) = [Pi11(i,j) Pi12(i,j) Pi13(i,j);
                           Pi12(i,j) Pi22(i,j) Pi23(i,j); 
                           Pi13(i,j) Pi23(i,j) Pi33(i,j)];
   end

end    
   
PS_cyl = cyl_Rey_ten_R(PS_car,theta,n_r,n_theta); 

PSij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      PSij(:,:,i) = PSij(:,:,i) + PS_cyl(:,:,i,j)./n_theta;
    end
end

PSrr = squeeze(PSij(1,1,:));
PStt = squeeze(PSij(2,2,:));
PSzz = squeeze(PSij(3,3,:));
PSrt = squeeze(PSij(1,2,:));
PSrz = squeeze(PSij(1,3,:));
PStz = squeeze(PSij(2,3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Viscous diffusion tensor VDij, as a function of r %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% VDrr = 1/r*(d[r*nu*d<urur>/dr]/dr) + 2*nu/r2*(<utut> - <urur>) =  nu*{[d2(<urur>)/dr2 + (1/r)*d<urur>/dr] + (2/r2)*(<utut> - <urur>)]} 
% VDtt = 1/r*(d[r*nu*d<utut>/dr]/dr) - 2*nu/r2*(<utut> - <urur>) =  nu*{[d2(<utut>)/dr2 + (1/r)*d<utut>/dr] - (2/r2)*(<utut> - <urur>)]} 
% VDzz = 1/r*(d[r*nu*d<uzuz>/dr]/dr)                             =  nu*{[d2(<uzuz>)/dr2 + (1/r)*d<uzuz>/dr]} 
% VDrt = 0
% VDrz = 1/r*(d[r*nu*d<uruz>/dr]/dr) - 1*nu/r2*<uruz>            =  nu*{[d2(<uruz>)/dr2 +(1/r)*d<uruz>/dr] - (1/r2)*<uruz>} 
% VDtz = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Apply Chain rule for Higher derivatives of multivariable functions  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

VD_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for it = 1:3
for jt = 1:3

[d2QdT2] = d2Chain(n_r,n_theta,r01,ssp,DM1,DM2,Dfr,squeeze(RR2_car(it,jt,:,:)));

VD_car(it,jt,:,:) = d2QdT2;

clear dQ2dT2
     
end
end

VD_cyl = cyl_Rey_ten_R(VD_car,theta,n_r,n_theta); 

VDij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      VDij(:,:,i) = VDij(:,:,i) + VD_cyl(:,:,i,j)./n_theta;
    end
end

VDrr = nu*squeeze(VDij(1,1,:));
VDtt = nu*squeeze(VDij(2,2,:));
VDzz = nu*squeeze(VDij(3,3,:));
VDrt = nu*squeeze(VDij(1,2,:));
VDrz = nu*squeeze(VDij(1,3,:));
VDtz = nu*squeeze(VDij(2,3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Viscousdissipation tensor Dij, as a function of r %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 42. e11: <(du/dx.du/dx + du/dy.du/dy + du/dz.du/dz)> (u: instantaneous)    % Q42 
% 43. e22: <(dv/dx.dv/dx + dv/dy.dv/dy + dv/dz.dv/dz)> (v: instantaneous)    % Q43 
% 44. e33: <(dw/dx.dw/dx + dw/dy.dw/dy + dw/dz.dw/dz)> (w: instantaneous)    % Q44 
% 45. e12: <(du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz)> (u, v: instantaneous) % Q45  
% 46. e13: <(du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz)> (u, w: instantaneous) % Q46 
% 47. e23: <(dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz)> (v, w: instantaneous) % Q47 

% e11 = e11 - (dU/dx)^2 - (dU/dy)^2 
% e22 = e22 - (dV/dx)^2 - (dV/dy)^2 
% e33 = e33 - (dW/dx)^2 - (dW/dy)^2 
% e12 = e12 - (dU/dx*dV/x) - (dU/dx*dV/dy)
% e13 = e13 - (dU/dx*dW/x) - (dU/dx*dW/dy)
% e23 = e23 - (dV/dx*dW/x) - (dV/dx*dW/dy)

e11(1:n_r,1:n_theta) = 0; 
e22(1:n_r,1:n_theta) = 0;
e33(1:n_r,1:n_theta) = 0;
e12(1:n_r,1:n_theta) = 0;
e13(1:n_r,1:n_theta) = 0;
e23(1:n_r,1:n_theta) = 0;

for j = 1:n_theta
        
    e11(:,j) = Q42(:,j) - ...
        (squeeze(dUidxj(1,1,:,j))).^2  - (squeeze(dUidxj(1,2,:,j))).^2; 
    e22(:,j) = Q43(:,j) - ...
        (squeeze(dUidxj(2,1,:,j))).^2  - (squeeze(dUidxj(2,2,:,j))).^2; 
    e33(:,j) = Q44(:,j) - ...
        (squeeze(dUidxj(3,1,:,j))).^2  - (squeeze(dUidxj(3,2,:,j))).^2; 
    e12(:,j) = Q45(:,j) - ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(2,1,:,j)) - ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(2,2,:,j));
    e13(:,j) = Q46(:,j) - ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(3,1,:,j)) - ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(3,2,:,j));
    e23(:,j) = Q47(:,j) - ...
        squeeze(dUidxj(2,1,:,j)).*squeeze(dUidxj(3,1,:,j)) - ...
        squeeze(dUidxj(2,2,:,j)).*squeeze(dUidxj(3,2,:,j));    
      
end

DS_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
        DS_car(:,:,i,j) = [e11(i,j) e12(i,j) e13(i,j);
                           e12(i,j) e22(i,j) e23(i,j); 
                           e13(i,j) e23(i,j) e33(i,j)];
   end

end    
   
DS_cyl = cyl_Rey_ten_R(DS_car,theta,n_r,n_theta); 

DSij(1:3,1:3,1:n_r) = 0; 

for i = 1:n_r
    for j = 1:n_theta 
      DSij(:,:,i) = DSij(:,:,i) + DS_cyl(:,:,i,j)./n_theta;
    end
end

Drr = -2*nu*squeeze(DSij(1,1,:));
Dtt = -2*nu*squeeze(DSij(2,2,:));
Dzz = -2*nu*squeeze(DSij(3,3,:));
Drt = -2*nu*squeeze(DSij(1,2,:));
Drz = -2*nu*squeeze(DSij(1,3,:));
Dtz = -2*nu*squeeze(DSij(2,3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%% Sum of budget terms for each Reynolds stress tensor %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SUM_rr = Prr + TDrr + PDrr + PSrr + VDrr + Drr;
SUM_tt = Ptt + TDtt + PDtt + PStt + VDtt + Dtt;
SUM_zz = Pzz + TDzz + PDzz + PSzz + VDzz + Dzz;
SUM_rt = Prt + TDrt + PDrt + PSrt + VDrt + Drt;
SUM_rz = Prz + TDrz + PDrz + PSrz + VDrz + Drz;
SUM_tz = Ptz + TDtz + PDtz + PStz + VDtz + Dtz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%% (omega)rms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 57. <omegar>          (omegar  : instantaneous) % Q57
% 58. <omegat>          (omegat  : instantaneous) % Q58 
% 59. <omegaz>          (omegaZ  : instantaneous) % Q59 

% 60. <omegar*omegar>   (omegar  : instantaneous) % Q60
% 61. <omegat*omegat>   (omegat  : instantaneous) % Q61 
% 62. <omegaz*omegaz>   (omegaz  : instantaneous) % Q62 

omr_mean    = sum(Q57,2)/n_theta;
omt_mean    = sum(Q58,2)/n_theta;
omz_mean    = sum(Q59,2)/n_theta;
omromr_mean = sum(Q60,2)/n_theta;
omtomt_mean = sum(Q61,2)/n_theta;
omzomz_mean = sum(Q62,2)/n_theta;

omrfomrf_mean = omromr_mean - omr_mean.*omr_mean;
omtfomtf_mean = omtomt_mean - omt_mean.*omt_mean;
omzfomzf_mean = omzomz_mean - omz_mean.*omz_mean;

omr_rms = nu*sqrt(omrfomrf_mean)./u_tau^2;
omt_rms = nu*sqrt(omtfomtf_mean)./u_tau^2;
omz_rms = nu*sqrt(omzfomzf_mean)./u_tau^2;

tauw_rms = omt_rms(end);
taut_rms = omz_rms(end);

ZB_W = sqrt(squeeze(R2(3,3,rr-1)))./R1(rr-1,3);
ZB_T = sqrt(squeeze(R2(2,2,rr-1)))./R1(rr-1,3);

duzrmsdz = flipdim(-DM1*sqrt(squeeze(R2(3,3,ssp))),1);
dutrmsdz = flipdim(-DM1*sqrt(squeeze(R2(2,2,ssp))),1);
 
ZB_WHOP = duzrmsdz(2:end)./DR1(2:end,3);
ZB_THOP = dutrmsdz(2:end)./DR1(2:end,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

