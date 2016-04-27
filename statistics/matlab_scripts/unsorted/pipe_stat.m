%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Compute and plot the  statistics for %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   turbulent pipe flow simulation     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; 
%clear all; 
%clc
double precision;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% The record in the binary file: stat00000N is saved as follows %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  <U>                              % Q
% 2.  <V>                              % Q1
% 3.  <W>                              % Q2
% 4.  <P>                              % Q3

% 5.  <uu> (u: instantaneous)          % Q4
% 6.  <vv> (v: instantaneous)          % Q5
% 7.  <ww> (w: instantaneous)          % Q6
% 8.  <pp> (p: instantaneous)          % Q7

% 9.  <uv> (u, v: instantaneous)       % Q8
% 10. <vw> (v, w: instantaneous)       % Q9
% 11. <uw> (u, w: instantaneous)       % Q10

% 12. <pu> (p, u: instantaneous)       % Q11
% 13. <pv> (p, v: instantaneous)       % Q12
% 14. <pw> (p, w: instantaneous)       % Q13

% 15. <pdudx> (p, dudx: instantaneous) % Q14
% 16. <pdudy> (p, dudy: instantaneous) % Q15
% 17. <pdudz> (p, dudz: instantaneous) % Q16

% 18. <pdvdx> (p, dvdx: instantaneous) % Q17
% 19. <pdvdy> (p, dvdy: instantaneous) % Q18
% 20. <pdvdz> (p, dvdz: instantaneous) % Q19

% 21. <pdwdx> (p, dwdx: instantaneous) % Q20
% 22. <pdwdy> (p, dwdy: instantaneous) % Q21
% 23. <pdwdz> (p, dwdz: instantaneous) % Q22

% 24. <uuu> (u: instantaneous)         % Q23
% 25. <vvv> (v: instantaneous)         % Q24
% 26. <www> (w: instantaneous)         % Q25
% 27. <ppp> (p: instantaneous)         % Q26

% 28. <uuv> (u, v: instantaneous)      % Q27
% 29. <uuw> (u, w: instantaneous)      % Q28
% 30. <vvu> (v, u: instantaneous)      % Q29
% 31. <vvw> (v, w: instantaneous) 	   % Q30   
% 32. <wwu> (w, u: instantaneous)      % Q31
% 33. <wwv> (w, v: instantaneous)      % Q32
% 34. <uvw> (u, v, w: instantaneous)   % Q33

% 35. <uuuu> (u, v: instantaneous)     % Q34
% 36. <vvvv> (u, w: instantaneous)     % Q35
% 37. <wwww> (v, u: instantaneous)     % Q36
% 38. <pppp> (v, w: instantaneous)     % Q37

% 39. <uuuv> (u: instantaneous)        % Q38
% 40. <uuvv> (v: instantaneous)        % Q39
% 41. <uvvv> (w: instantaneous)	       % Q40 
   
% 42. e11: <(d2u2/dx2 + d2u2/dy2 + d2u2/dz2)> (u: instantaneous)             % Q41 
% 43. e22: <(d2v2/dx2 + d2v2/dy2 + d2v2/dz2)> (v: instantaneous)             % Q42 
% 44. e33: <(d2w2/dx2 + d2w2/dy2 + d2w2/dz2)> (w: instantaneous)             % Q43 
% 45. e12: <(du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz)> (u, v: instantaneous) % Q44  
% 46. e13: <(du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz)> (u, w: instantaneous) % Q45 
% 47. e23: <(dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz)> (v, w: instantaneous) % Q46 

% 48. <omz>     (omz: instantaneous)     % Q47
% 49. <omz*omz> (omz: instantaneous)     % Q48

% 50. <dw/dx*dw/dx> (w: instantaneous)	 % Q49 
% 51. <dw/dy*dw/dy> (w: instantaneous)	 % Q50 
% 52. <dw/dx*dw/dy> (w: instantaneous)	 % Q51 

% 53. <du/dx*du/dx> (u: instantaneous)	 % Q52
% 54. <du/dy*du/dy> (u: instantaneous)	 % Q53 
% 55. <du/dx*du/dy> (u: instantaneous)	 % Q54 

% 56. <dv/dx*dv/dx> (v: instantaneous)	 % Q55
% 57. <dv/dy*dv/dy> (v: instantaneous)	 % Q56 
% 58. <dv/dx*dv/dy> (v: instantaneous)	 % Q57 

% 59. <du/dx*dv/dx> (u,v: instantaneous)	 % Q58
% 60. <du/dy*dv/dy> (u,v: instantaneous)	 % Q59 
% 61. <du/dx*dv/dy> (u,v: instantaneous)	 % Q60 
% 62. <du/dy*dv/dx> (u,v: instantaneous)	 % Q61

% 63. <du/dx>       (u  : instantaneous)	 % Q62
% 64. <du/dy>       (u  : instantaneous)	 % Q63 
% 65. <du/dz>       (u  : instantaneous)	 % Q64 

% 66. <dv/dx>       (v  : instantaneous)	 % Q65
% 67. <dv/dy>       (v  : instantaneous)	 % Q66 
% 68. <dv/dz>       (v  : instantaneous)	 % Q67 

% 69. <dw/dx>       (w  : instantaneous)	 % Q68
% 70. <dw/dy>       (w  : instantaneous)	 % Q69 
% 71. <dw/dz>       (w  : instantaneous)	 % Q70 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Read the interpolated field %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fname')
  fname = '../recordings/polar_x_0.0';
end
[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
CH            = fread(fid,hdr,'*char')    ;
dum1          = fread(fid,1,'*float64')   ;
Rer           = fread(fid,1,'*float64')   ;
Domain        = fread(fid,3,'*float64')   ;
nel           = fread(fid,3,'int32')      ;
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
clear('Q');
for i=1:70
  clear(['Q' num2str(i)]);
end


for i =1:nstat
    
if i == 1
   fseek(fid,0,'cof');
else
   fseek(fid,8,'cof');
end

 MatR = fread(fid,[ng_r, ng_theta],'*float64'); 
 Mat  = flipdim(MatR,1);
 
 v  = genvarname('Q', who);
 evalc([v '=INOUT(Mat)']);   
    
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reb = U(bulk)*D/nu 
% Ub: bulk velocity. 
% D : pipe diameter.
% nu: kinematic viscosity

Ub  = 1;
%D   = round(Domain(1));
D = 1;
Reb = Rer*D;
R   = D/2;
rho = 1;
nu  = 1/Rer; 
mu  = nu*rho;

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
r10  = R - r01     ; %%%%% grid in radial direction from 1 (center) to 0 (wall)
rr = length(r10);      %%%%% 

atheta   = theta*180/pi;
n_theta  = ng_theta;  
n_r      = ng_r;

ssp = rr:-1:1; 

[th,r_] = meshgrid(theta,r01);
[X,Y] = pol2cart(th,r_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Straight pipe data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid = fopen('ZZ_Budget_180.dat')               ;
% FX  = textscan(fid, '%f %f %f %f %f %f %f %f','headerlines',16);
% fclose(fid)                                              ;
% 
% r10_r10      = FX{1};
% prms_1000       = FX{2};
% omr_rms_1000    = FX{3};
% omt_rms_1000    = FX{4};
% omz_rms_1000    = FX{5};
% clear FX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Bendiks Jan Boersma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nr_bjb      = 200             ;
% nt_bjb      = 360             ;
% nz_bjb      = 640             ;
% Re_tau_bjb  = 1100.00000000000;     
% Re_bulk_bjb = 18931.0235426785;     
% Re_cen_bjb  = 23846.1293443586; 
% 
% fid = fopen('DATA_BJB/RET_550/mean_01100');
% FX  = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',28);
% fclose(fid);
% 
% rp05_bjb        = FX{1} ; 
% rp_bjb          = FX{2} ; 
% U_bjb           = FX{3} ; 
% V_bjb           = FX{4} ; 
% W_bjb           = FX{5} ;
% P_bjb           = FX{6} ; 
% uu_mean_bjb     = FX{7} ; 
% vv_mean_bjb     = FX{8} ;
% ww_mean_bjb     = FX{9} ;
% pp_mean_bjb     = FX{10};
% uw_mean_bjb     = FX{11};
% tau_visc_bjb    = FX{12}; % mu*dw/dr,
% tau_t_bjb       = FX{13};
% tau_viscRr_bjb  = FX{14}; %-(R-r)*dw/dr
% uuu_mean_bjb    = FX{15};
% vvv_mean_bjb    = FX{16};
% www_mean_bjb    = FX{17};
% uuuu_mean_bjb   = FX{18};
% vvvv_mean_bjb   = FX{19};
% wwww_mean_bjb   = FX{20};    
% 
% clear FX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Wu & Moin-data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ub_wm       =  1.00000000E+00;
% tau_wall_wm =  3.10475214E-03;
% Re_tau_wm   =  6.84802521E+02;
% delta_xp    =  1.00312869E+01;
% delta_yp    =  1.43759831E-01;
% delta_zp    =  4.20101354E+00;
% 
% u_tau_wm   = sqrt(tau_wall_wm);
% %------------------------
% fid = fopen('DATA_MOIN/RE_24k/1_re24k_stat.txt');
% FX  = textscan(fid, '%f %f %f %f %f %f %f','headerlines',3);
% fclose(fid);
% 
% ym1_wm = FX{1}; 
% U_wm   = FX{2}; 
% uf_wm  = FX{3}; 
% W_wm   = FX{3};
% wf_wm  = FX{5}; 
% P_wm   = FX{6}; 
% pf_wm  = FX{7};
% 
% ym1_wm1 = 1 - ym1_wm;
% clear FX
% %------------------------      
% fid = fopen('DATA_MOIN/RE_24k/2_re24k_stat.txt');
% FX  = textscan(fid, '%f %f %f %f %f %f','headerlines',3);  
% fclose(fid);
% 
% y2_wm       = FX{1}; 
% V_wm        = FX{2}; 
% vf_wm       = FX{3}; 
% uv_wm       = FX{4};
% tau_sgs_wm  = FX{5}; 
% tau_visc_wm = FX{6};
%       
% y2_wm1 = 1 - y2_wm;
% clear FX
% %------------------------   
% fid = fopen('DATA_MOIN/RE_24k/3_re24k_stat.txt');
% FX  = textscan(fid, '%f %f %f %f %f %f %f %f','headerlines',3);        
% fclose(fid);
% 
% ymp3_wm = FX{1}; 
% Up_wm   = FX{2}; 
% ufp_wm  = FX{3}; 
% Wp_wm   = FX{4};
% wfp_wm  = FX{5}; 
% Pp_wm   = FX{6}; 
% pfp_wm  = FX{7}; 
% loglaw3 = FX{8};       
% 
% clear FX
% %------------------------  
% fid = fopen('DATA_MOIN/RE_24k/4_re24k_stat.txt');
% FX  = textscan(fid,'%f %f %f %f %f %f %f','headerlines',3);        
% fclose(fid);
% 
% ymp4_wm      = FX{1}; 
% Vp_wm        = FX{2}; 
% vfp_wm       = FX{3}; 
% uvp_wm       = FX{4}; 
% taup_sgs_wm  = FX{5}; 
% taup_visc_wm = FX{6}; 
% convgp       = FX{7}; 
% 
% clear FX
% %------------------------ 
% fid = fopen('DATA_MOIN/RE_24k/5_re24k_stat.txt');
% FX  = textscan(fid,'%f %f %f %f %f','headerlines',3);        
% fclose(fid);
% 
% y5_wm       = FX{1};
% ym5_wm      = FX{2}; 
% nut_wm      = FX{3}; 
% theta_uv_wm = FX{4}; 
% residual_wm = FX{5}; 
% 
% clear FX

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


RR1_car(1:3,1:n_r,1:n_theta) = 0;

for i = 1:n_r
   for j=1:n_theta
       RR1_car(:,i,j) = [Q(i,j); Q1(i,j); Q2(i,j)];
    end
end


RR1 = cyl_Vel_R(RR1_car,theta,n_r,n_theta); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dUidxj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

dUidxj(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
         dUidxj(:,:,i,j) = [Q62(i,j) Q63(i,j) Q64(i,j);
                            Q65(i,j) Q66(i,j) Q67(i,j); 
                            Q68(i,j) Q69(i,j) Q70(i,j)];
    end
end   


dUidxj_cyl = cyl_Rey_ten_R(dUidxj,theta,n_r,n_theta);
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
% DR1(1:n_r,1:3)            = 0;

for i = 1:3
for j = 1:n_theta
    
    DRR1(:,j,i) = squeeze(dUidxj(i,1,:,j)).*cos(theta(j)) + ...
                  squeeze(dUidxj(i,2,:,j)).*sin(theta(j));
        
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Velocity defect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Uzd   = R1(1,3) - R1(:,3); 
%U_wmd = U_wm(1) - U_wm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Shaper factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% The shape factor in boundary layer flows is defined as:
%%% H12 = δ1/δ2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% δ1: is the displacement thickness  
%%% δ1 = int(0,inf)[1 - u(y)/uo]dy 
%%% uo is the velocity in the free-stream outside the boundary layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% δ2: is the momentum thickness 
%%% δ2 = int(0,inf)[u(y)/uo*{1 - u(y)/uo}]dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% for i =1:rr
%     if R1(i,3) < 0.1*R1(1,3)/Ub
%        Uo = R1(i,3);
%        break   
%     end
%     
% end
% 
% udel1 = 1 - R1(end:-1:1,3)/R1(1,3);
% rdel1 = 1 - r01(end:-1:1);
% 
% udel2 = (R1(end:-1:1,3)/R1(1,3)).*udel1;
% 
% del1 = trapz(rdel1,udel1);
% del2 = trapz(rdel1,udel2);
% 
% H12 = del1/del2;

%%%%% plot(rdel1,udel1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Reynolds stress tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

RR2_car(1:3,1:3,1:n_r,1:n_theta) = 0;
 
for j = 1:n_theta
   for i = 1:n_r
        RR2_car(:,:,i,j) = [(Q4(i,j) -Q(i,j).*Q(i,j) ) (Q8(i,j)-Q(i,j).*Q1(i,j) )  (Q10(i,j)-Q(i,j).*Q2(i,j) );
                            (Q8(i,j) -Q(i,j).*Q1(i,j)) (Q5(i,j)-Q1(i,j).*Q1(i,j))  (Q9(i,j) -Q1(i,j).*Q2(i,j)); 
                            (Q10(i,j)-Q(i,j).*Q2(i,j)) (Q9(i,j)-Q1(i,j).*Q2(i,j))  (Q6(i,j) -Q2(i,j).*Q2(i,j))];
    end
end    

RR2 = cyl_Rey_ten_R(RR2_car,theta,n_r,n_theta); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Total shear stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%tau_rt  = squeeze(R2(3,1,:))./tau_wall - DR1(:,3).*(l_visc/u_tau);
%taup_wm = uvp_wm + taup_visc_wm;
%tau_wm  = uv_wm  + tau_visc_wm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Residual Rsd = tau_rt - r/R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%Rsd    = tau_rt - r01/R;
%Rsd_wm = taup_wm(end:-1:1) - y5_wm(1:end-1)/R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Mean, r.m.s, Skewness and Flatness of pressure %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

P_mean    = Q3;
pp_mean   = Q7;
ppp_mean  = Q26;
pppp_mean = Q37;

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

% P    = sum(P_mean,2       )./n_theta;
% pp   = sum(pfpf_mean,2    )./n_theta; 
% ppp  = sum(pfpfpf_mean,2  )./n_theta; 
% pppp = sum(pfpfpfpf_mean,2)./n_theta; 

P    = (P_mean);
pp   = (pfpf_mean); 
ppp  = (pfpfpf_mean); 
pppp = (pfpfpfpf_mean); 



prms = sqrt(pp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

skew_p = ppp./(pp).^(3/2);
flat_p = pppp./(pp).^(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Skewness and Flatness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
% Q23 = uuu_mean;
% Q24 = vvv_mean;
% Q25 = www_mean;
% 
% Q27 = uuv_mean;
% Q28 = uuw_mean;
% Q29 = vvu_mean;
% Q30 = vvw_mean;
% Q31 = wwu_mean;
% Q32 = wwv_mean;
% Q33 = uvw_mean;

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

tensRR3_car(:,:,1,i,j) = [Q23(i,j) Q27(i,j) Q28(i,j); 
                          Q27(i,j) Q29(i,j) Q33(i,j); 
                          Q28(i,j) Q33(i,j) Q31(i,j)];
             
tensRR3_car(:,:,2,i,j) = [Q27(i,j) Q29(i,j) Q33(i,j); 
                          Q29(i,j) Q24(i,j) Q30(i,j); 
                          Q33(i,j) Q30(i,j) Q32(i,j)];
             
tensRR3_car(:,:,3,i,j) = [Q28(i,j) Q33(i,j) Q31(i,j); 
                          Q33(i,j) Q30(i,j) Q32(i,j); 
                          Q31(i,j) Q32(i,j) Q25(i,j)];   
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
 
% for i = 1:n_r
%  for j = 1:n_theta 
%    tensR3(:,:,:,i)  = tensR3(:,:,:,i) + tensRR3(:,:,:,i,j)./n_theta;
%  end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Skewness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
skew_tol = 0;
skew_ur  = squeeze(tensRR3(1,1,1,1:end-skew_tol,:))./((squeeze(RR2(1,1,1:end-skew_tol,:))).^(3/2));
skew_uth = squeeze(tensRR3(2,2,2,1:end-skew_tol,:))./((squeeze(RR2(2,2,1:end-skew_tol,:))).^(3/2));
skew_us  = squeeze(tensRR3(3,3,3,1:end-skew_tol,:))./((squeeze(RR2(3,3,1:end-skew_tol,:))+2e-20).^(3/2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Flatness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Q34 = uuuu_mean;
% Q35 = vvvv_mean;
% Q36 = wwww_mean;
% 
% Q38 = uuuv_mean;
% Q39 = uuvv_mean;
% Q40 = uvvv_mean;

urururur_inst(1:n_r,1:n_theta) = 0;
utututut_inst(1:n_r,1:n_theta) = 0;

for j = 1:n_theta 

urururur_inst(:,j) = Q34(:,j).*cos(theta(j)).^4 + ...
                   6*Q39(:,j).*cos(theta(j)).^2.*sin(theta(j)).^2 + ...
                   4*Q38(:,j).*cos(theta(j)).^3.*sin(theta(j))    + ...
                   4*Q40(:,j).*cos(theta(j)).*sin(theta(j)).^3    + ...
                     Q35(:,j).*sin(theta(j)).^4;

utututut_inst(:,j) = Q34(:,j).*sin(theta(j)).^4 + ...
                   6*Q39(:,j).*cos(theta(j)).^2.*sin(theta(j)).^2 - ...
                   4*Q38(:,j).*cos(theta(j)).*sin(theta(j)).^3    - ...
                   4*Q40(:,j).*cos(theta(j)).^3.*sin(theta(j))    + ...
                     Q35(:,j).*cos(theta(j)).^4;            
            
end

urururur_2d = urururur_inst - 4*squeeze(tensRR3(1,1,1,:,:)).*squeeze(RR1(1,:,:)) - ...
     6*squeeze(RR2(1,1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)) - ... 
       squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)).*squeeze(RR1(1,:,:)); 

utututut_2d = utututut_inst - 4*squeeze(tensRR3(2,2,2,:,:)).*squeeze(RR1(2,:,:)) - ...
              6*squeeze(RR2(2,2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)) - ...
        squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:)).*squeeze(RR1(2,:,:));           
          
uzuzuzuz_2d = Q36 - 4*squeeze(tensRR3(3,3,3,:,:)).*squeeze(RR1(3,:,:)) - ...
      6*squeeze(RR2(3,3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)) - ...
        squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)).*squeeze(RR1(3,:,:)); 
          
% urururur = sum(urururur_2d,2)./n_theta;   
% utututut = sum(utututut_2d,2)./n_theta;  
% uzuzuzuz = sum(uzuzuzuz_2d,2)./n_theta;        

flat_tol = 0;

flat_ur  = urururur_2d(1:end-flat_tol,:)./(squeeze(RR2(1,1,1:end-flat_tol,:)).^(2));
flat_uth = utututut_2d(1:end-flat_tol,:)./(squeeze(RR2(2,2,1:end-flat_tol,:)).^(2));
flat_us  = uzuzuzuz_2d(1:end-flat_tol,:)./(squeeze(RR2(3,3,1:end-flat_tol,:)).^(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Convection tensor Cij, as a function of r %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Cij = Uk*d<uiuj>/dxk
% C11 = Uk*(d<uuk>)/dxk, K=1,3 = [U*d<uu>/dx + V*d<uu>/dy + W*d<uu>/dz]
% C22 = Uk*(d<vvk>)/dxk, K=1,3 = [U*d<vv>/dx + V*d<vv>/dy + W*d<vv>/dz]
% C33 = Uk*(d<wwk>)/dxk, K=1,3 = [U*d<ww>/dx + V*d<ww>/dy + W*d<ww>/dz]
% C12 = Uk*(d<uvk>)/dxk, K=1,3 = [U*d<uv>/dx + V*d<uv>/dy + W*d<uv>/dz]
% C13 = Uk*(d<uwk>)/dxk, K=1,3 = [U*d<uw>/dx + V*d<uw>/dy + W*d<uw>/dz]
% C23 = Uk*(d<vwk>)/dxk, K=1,3 = [U*d<vw>/dx + V*d<vw>/dy + W*d<vw>/dz]

% Crr = ... 
% Ctt = ...
% Czz = ...
% Crt = ...
% Crz = ...
% Ctz = ...

C_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for it = 1:3
    for jt = 1:3
        
        [dRijdxj] = dQ2dxj_t(it,jt,n_r,n_theta,r01,theta,ssp,DM1,Dfr,RR2_car);
        
        for dummy = 1:3
            
            C_car(it,jt,:,:) = squeeze(C_car(it,jt,:,:)) + squeeze(dRijdxj(dummy,:,:)) ...
                .* squeeze(RR1_car(dummy,:,:));
        end
        
    end
end

C_cyl = cyl_Rey_ten_R(C_car,theta,n_r,n_theta); 

% Cij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       Cij(:,:,i) = Cij(:,:,i) + C_cyl(:,:,i,j)./n_theta;
%     end
% end


Crr   = squeeze(C_cyl(1,1,:,:));%./normal;
Cthth = squeeze(C_cyl(2,2,:,:));%./normal;
Css   = squeeze(C_cyl(3,3,:,:));%./normal;
Crth  = squeeze(C_cyl(1,2,:,:));%./normal;
Crs   = squeeze(C_cyl(1,3,:,:));%./normal;
Csth  = squeeze(C_cyl(2,3,:,:));%./normal;

Ctke = 0.5*(Crr+Cthth+Css);

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

% PRij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       PRij(:,:,i) = PRij(:,:,i) + PR_cyl(:,:,i,j)./n_theta;
%     end
% end

Prr   = squeeze(PR_cyl(1,1,:,:));%./normal;
Pthth = squeeze(PR_cyl(2,2,:,:));%./normal;
Pss   = squeeze(PR_cyl(3,3,:,:));%./normal;
Prth  = squeeze(PR_cyl(1,2,:,:));%./normal;
Prs   = squeeze(PR_cyl(1,3,:,:));%./normal;
Psth  = squeeze(PR_cyl(2,3,:,:));%./normal;

Ptke = 0.5*(Prr+Pthth+Pss);
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

% TDij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       TDij(:,:,i) = TDij(:,:,i) + TD_cyl(:,:,i,j)./n_theta;
%     end
% end

TDrr   = squeeze(TD_cyl(1,1,:,:));%./normal;
TDthth = squeeze(TD_cyl(2,2,:,:));%./normal;
TDss   = squeeze(TD_cyl(3,3,:,:));%./normal;
TDrth  = squeeze(TD_cyl(1,2,:,:));%./normal;
TDrs   = squeeze(TD_cyl(1,3,:,:));%./normal;
TDsth  = squeeze(TD_cyl(2,3,:,:));%./normal;

TDtke = 0.5*(TDrr+TDthth+TDss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Pressure diffusion tensor PDij, as a function of r %%%%%%%%%%%%%%%%
%%%%%%%%%%% Divergence of pressure-velocity correlation %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 12. pu_mean (p, u: instantaneous)       % Q11
% 13. pv_mean (p, v: instantaneous)       % Q12
% 14. pw_mean (p, w: instantaneous)       % Q13

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
       PU_car(:,i,j) = [Q11(i,j) - Q3(i,j).*Q(i,j)  ;...
                        Q12(i,j) - Q3(i,j).*Q1(i,j) ;...
                        Q13(i,j) - Q3(i,j).*Q2(i,j)];
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

% PDij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       PDij(:,:,i) = PDij(:,:,i) + PD_cyl(:,:,i,j)./n_theta;
%     end
% end

PDrr   = squeeze(PD_cyl(1,1,:,:));%./normal;
PDthth = squeeze(PD_cyl(2,2,:,:));%./normal;
PDss   = squeeze(PD_cyl(3,3,:,:));%./normal;
PDrth  = squeeze(PD_cyl(1,2,:,:));%./normal;
PDrs   = squeeze(PD_cyl(1,3,:,:));%./normal;
PDsth  = squeeze(PD_cyl(2,3,:,:));%./normal;

PDtke = 0.5*(PDrr+PDthth+PDss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Pressure strain tensor PSij, as a function of r %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 15. pdudx_mean (p, dudx: instantaneous) % Q14
% 16. pdudy_mean (p, dudy: instantaneous) % Q15
% 17. pdudz_mean (p, dudz: instantaneous) % Q16

% 18. pdvdx_mean (p, dvdx: instantaneous) % Q17
% 19. pdvdy_mean (p, dvdy: instantaneous) % Q18
% 20. pdvdz_mean (p, dvdz: instantaneous) % Q19

% 21. pdwdx_mean (p, dwdx: instantaneous) % Q20
% 22. pdwdy_mean (p, dwdy: instantaneous) % Q21
% 23. pdwdz_mean (p, dwdz: instantaneous) % Q22

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
        
    Pi11(:,j) = (2/rho)*(Q14(:,j) - Q3(:,j).*(squeeze(dUidxj(1,1,:,j)))); 
    Pi22(:,j) = (2/rho)*(Q18(:,j) - Q3(:,j).*(squeeze(dUidxj(2,2,:,j)))); 
    Pi33(:,j) = (2/rho)*(Q22(:,j) - Q3(:,j).*(squeeze(dUidxj(3,3,:,j)))); 
    
    Pi12(:,j) = (1/rho)*(Q15(:,j) + Q17(:,j)   - ...
        Q3(:,j).*(squeeze(dUidxj(1,2,:,j)))    - ...
        Q3(:,j).*(squeeze(dUidxj(2,1,:,j))) ); 
    Pi13(:,j) = (1/rho)*(Q16(:,j) + Q20(:,j)   - ...
        Q3(:,j).*(squeeze(dUidxj(1,3,:,j)))    - ...
        Q3(:,j).*(squeeze(dUidxj(3,1,:,j))) ); 
    Pi23(:,j) = (1/rho)*(Q19(:,j) + Q21(:,j)   - ...
        Q3(:,j).*(squeeze(dUidxj(2,3,:,j)))    - ...
        Q3(:,j).*(squeeze(dUidxj(3,2,:,j))) ); 
    
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

% PSij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       PSij(:,:,i) = PSij(:,:,i) + PS_cyl(:,:,i,j)./n_theta;
%     end
% end

PSrr   = squeeze(PS_cyl(1,1,:,:));%./normal;
PSthth = squeeze(PS_cyl(2,2,:,:));%./normal;
PSss   = squeeze(PS_cyl(3,3,:,:));%./normal;
PSrth  = squeeze(PS_cyl(1,2,:,:));%./normal;
PSrs   = squeeze(PS_cyl(1,3,:,:));%./normal;
PSsth  = squeeze(PS_cyl(2,3,:,:));%./normal;

PStke = 0.5*(PSrr+PSthth+PSss);
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
%%%
%%% d(W)/dx2 = dW/dr*(sin(theta)^2)/r                 + 
%%%            cos(theta)^2*d2(W)/dr2                 - 
%%%            cos(theta)*sin(theta)/r                -
%%%            cos(theta)*sin(theta)*d2(W)/drd(theta) - 
%%%            2*sin(theta)*cos(theta)*dW/d(theta)    +
%%%            cos(theta)*d2(W)/drd(theta)            -
%%%            sin(theta)/r*d2W/d(theta)2

VD_car(1:3,1:3,1:n_r,1:n_theta) = 0;

for it = 1:3
for jt = 1:3

    
 [d2Qdx2,d2Qdy2] = d2Chain(n_r,n_theta,r01,theta,ssp,DM1,DM2,Dfr,...
     squeeze(RR2_car(it,jt,:,:)));
% %%%%%%%%%%%%%%%
% [d2Qdx2,d2Qdy2] = d2Chain(n_r,n_theta,r01,theta,ssp,DM1,DM2,Dfr,...
%     squeeze(RR2_car(3,3,:,:)));
% %%%%%%%%%%%%%%%

VD_car(it,jt,:,:) = d2Qdx2 + d2Qdy2;

clear d2Qdx2 d2Qdy2 
 
end
end




VD_cyl = cyl_Rey_ten_R(VD_car,theta,n_r,n_theta); 

% VDij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       VDij(:,:,i) = VDij(:,:,i) + VD_cyl(:,:,i,j)./n_theta;
%     end
% end

VDrr   = nu*squeeze(VD_cyl(1,1,:,:));%./normal;
VDthth = nu*squeeze(VD_cyl(2,2,:,:));%./normal;
VDss   = nu*squeeze(VD_cyl(3,3,:,:));%./normal;
VDrth  = nu*squeeze(VD_cyl(1,2,:,:));%./normal;
VDrs   = nu*squeeze(VD_cyl(1,3,:,:));%./normal;
VDsth  = nu*squeeze(VD_cyl(2,3,:,:));%./normal;

VDtke = 0.5*(VDrr+VDthth+VDss);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Viscous dissipation tensor Dij, as a function of r %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 42. e11: (d2u2/dx2 + d2u2/dy2 + d2u2/dz2)_mean (u: instantaneous)             % Q41 
% 43. e22: (d2v2/dx2 + d2v2/dy2 + d2v2/dz2)_mean (v: instantaneous)             % Q42 
% 44. e33: (d2w2/dx2 + d2w2/dy2 + d2w2/dz2)_mean (w: instantaneous)             % Q43 
% 45. e12: (du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz)_mean (u, v: instantaneous) % Q44  
% 46. e13: (du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz)_mean (u, w: instantaneous) % Q45 
% 47. e23: (dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz)_mean (v, w: instantaneous) % Q46 

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
        
    e11(:,j) = Q41(:,j) - ...
        (squeeze(dUidxj(1,1,:,j))).^2  - (squeeze(dUidxj(1,2,:,j))).^2; 
    e22(:,j) = Q42(:,j) - ...
        (squeeze(dUidxj(2,1,:,j))).^2  - (squeeze(dUidxj(2,2,:,j))).^2; 
    e33(:,j) = Q43(:,j) - ...
        (squeeze(dUidxj(3,1,:,j))).^2  - (squeeze(dUidxj(3,2,:,j))).^2; 
    e12(:,j) = Q44(:,j) - ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(2,1,:,j)) - ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(2,2,:,j));
    e13(:,j) = Q45(:,j) - ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(3,1,:,j)) - ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(3,2,:,j));
    e23(:,j) = Q46(:,j) - ...
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

% DSij(1:3,1:3,1:n_r) = 0; 
% 
% for i = 1:n_r
%     for j = 1:n_theta 
%       DSij(:,:,i) = DSij(:,:,i) + DS_cyl(:,:,i,j)./n_theta;
%     end
% end

Drr   = -2*nu*squeeze(DS_cyl(1,1,:,:));%./normal;
Dthth = -2*nu*squeeze(DS_cyl(2,2,:,:));%./normal;
Dss   = -2*nu*squeeze(DS_cyl(3,3,:,:));%./normal;
Drth  = -2*nu*squeeze(DS_cyl(1,2,:,:));%./normal;
Drs   = -2*nu*squeeze(DS_cyl(1,3,:,:));%./normal;
Dsth  = -2*nu*squeeze(DS_cyl(2,3,:,:));%./normal;

Dtke = 0.5*(Drr+Dthth+Dss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%% Sum of budget terms for each Reynolds stress tensor %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

SUMrr   = -Crr  + Prr  + TDrr  + PDrr  + PSrr  + VDrr  + Drr;
SUMthth = -Cthth+Pthth + TDthth+PDthth + PSthth+ VDthth+ Dthth;
SUMss   = -Css  + Pss  + TDss  + PDss  + PSss  + VDss  + Dss;
SUMrth  = -(-Crth + Prth + TDrth + PDrth + PSrth + VDrth + Drth);
SUMrs   = -(-Crs  + Prs  + TDrs  + PDrs  + PSrs  + VDrs  + Drs);
SUMsth  = -(-Csth + Psth + TDsth + PDsth + PSsth + VDsth + Dsth);

SUMtke  = -Ctke + Ptke + TDtke + PDtke + PStke + VDtke + Dtke;


SUMij(1:3,1:3,1:n_r,1:n_theta) = 0; 
for j = 1:n_theta
   for i = 1:n_r
        SUMij(:,:,i,j) =         [ SUMrr(i,j)  SUMrth(i,j)  SUMrs(i,j) ;
                                   SUMrth(i,j) SUMthth(i,j) SUMsth(i,j); 
                                   SUMrs(i,j)  SUMsth(i,j)  SUMss(i,j)];
   end

end    

sum_all(1:3,1:3,1:n_r) = 0;
for i = 1:n_r
    for j = 1:n_theta 
      sum_all(:,:,i) = sum_all(:,:,i) + SUMij(:,:,i,j)./n_theta;
    end
end


sum_tke(1:n_r) = 0;
for i = 1:n_r
    for j = 1:n_theta 
      sum_tke(i) = sum_tke(i) + SUMtke(i,j)./n_theta;
    end
end

