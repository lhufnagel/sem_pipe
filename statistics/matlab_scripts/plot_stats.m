%clear all;
addpath([pwd '/unsorted'])

% :-X hardcoded definitions..
Re_tau = 180;
%Re_tau = 1.8105409983122749E+02; % George's Reference data at this Re_tau
pipe_radius = 0.5;
nu = 1/5300;

u_plot_scaling = 0.08;
rms_plot_scaling = 2;
k_plot_scaling = 100.0;

%yplus/y =Re_tau/pipe_radius
delta_tau = pipe_radius/Re_tau;
%u_tau = nu/delta_tau;

u_tau = 6.8322301823104711E-02;

%reference = importdata('/scratch/hufnagel/MSc/ElKhouryData/180_Re_1.dat'); % El Khoury data
reference = importdata('../../../ElKhouryData/180_Re_1.dat'); % El Khoury data
ref_rr_budg = importdata('../../../ElKhouryData/180_RR_Budget.dat'); % El Khoury data
ref_tt_budg = importdata('../../../ElKhouryData/180_TT_Budget.dat'); % El Khoury data
ref_zz_budg = importdata('../../../ElKhouryData/180_ZZ_Budget.dat'); % El Khoury data

%ref_eps_rr = u_tau^4/nu*ref_rr_budg.data(:,8);
ref_eps_rr = u_tau^4/nu*ref_rr_budg.data(:,8);
ref_eps_tt = u_tau^4/nu*ref_tt_budg.data(:,8);
ref_eps_zz = u_tau^4/nu*ref_zz_budg.data(:,8);

ref_r = reference.data(:,2);
ref_u = reference.data(:,3);
ref_ur_rms = reference.data(:,5);
ref_ut_rms = reference.data(:,6);
ref_uz_rms = reference.data(:,7);
ref_k = u_tau^2*.5*(ref_ur_rms.^2+ref_ut_rms.^2+ref_uz_rms.^2);

%Glob files
file_names = dir('../recordings/polar_x_*');
x_vals = zeros(size(file_names));
if (length(x_vals) == 0)
    error('no files found')
end
for i =1:length(file_names)
  x_vals(i)=str2num(file_names(i).name(9:end));
end

%my rotation_tens = @(theta) [1, 0 0;0 cos(theta) sin(theta);0,-sin(theta), cos(theta)];

% Small function for relative error that avoids div-by-zero
%calcRelErr = @(ref_x, actual_x) norm((ref_x(ref_x>1e-10) - actual_x(ref_x>1e-10))./ref_x(ref_x>1e-10),1)/length(ref_x(ref_x>1e-10));
calcRelErr = @(ref_x, actual_x) norm((ref_x(ref_x>1e-10) - actual_x(ref_x>1e-10)))/norm(ref_x(ref_x>1e-10));

figure(1);
clf;
figure(2);
clf;
figure(3);
clf;
figure(4);
clf;
figure(5);
clf;
figure(6);
clf;
figure(7);
clf;
figure(8);
clf;


for file = 1:length(file_names)

  fname  = ['../recordings/' file_names(file).name];
  pipe_stat;

  figure(1);
  u_mean = mean(Q,2);

  %visc = ref_r;
  %visc = visc(visc < 30);
  %loglaw = ref_r;
  %loglaw = loglaw(loglaw > 5);

  semilogy(x_vals(file) + u_plot_scaling*u_mean/u_tau, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) + u_plot_scaling*ref_u, ref_r, 'k--');

  % Plot viscous sublayer
  %semilogy(x_vals(file) +  u_plot_scaling*visc,visc,'b--');

  %% Plot log-law
  %semilogy(x_vals(file)*ones(size(loglaw)) + u_plot_scaling*(1/0.41*log(loglaw)+5.0),loglaw,'r--');

  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, ref_u, r10/delta_tau), u_mean/u_tau);
  text(x_vals(file)+.5*max(u_mean),max(ref_r)*.01+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);


  figure(2);

  RR2_mean = mean(RR2,4);

  %take the abs() because unfortunately some values are <0
  ux_rms = sqrt(abs(squeeze(RR2_mean(1,1,:))));
  ut_rms = sqrt(abs(squeeze(RR2_mean(2,2,:))));
  ur_rms = sqrt(abs(squeeze(RR2_mean(3,3,:))));
  k = .5*(ux_rms.^2+ur_rms.^2+ut_rms.^2);


  semilogy(x_vals(file) +  k_plot_scaling*k, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) + k_plot_scaling*ref_k, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);
  err_rel = calcRelErr(interp1(ref_r, ref_k, r10/delta_tau), k);
  text(x_vals(file)+.5*max(ref_k),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(3);

  semilogy(x_vals(file) +  ux_rms/u_tau, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) + ref_uz_rms, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, ref_uz_rms, r10/delta_tau), ux_rms/u_tau);
  text(x_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);


  figure(4);

  semilogy(x_vals(file) +  rms_plot_scaling*ut_rms/u_tau, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) + rms_plot_scaling*ref_ut_rms, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, ref_ut_rms, r10/delta_tau), ut_rms/u_tau);
  text(x_vals(file)+.5*max(ref_ut_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(5);

  semilogy(x_vals(file) +  rms_plot_scaling*ur_rms/u_tau, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) + rms_plot_scaling*ref_ur_rms, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, ref_ur_rms, r10/delta_tau), ur_rms/u_tau);
  text(x_vals(file)+.5*max(ref_ur_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(6);



 %Drr   = -2*nu*squeeze(DS_cyl(1,1,:,:));%./normal;
 %Dthth = -2*nu*squeeze(DS_cyl(2,2,:,:));%./normal;
 %Dss   = -2*nu*squeeze(DS_cyl(3,3,:,:));%./normal;
  eps_rr = -mean(Dss,2);

  semilogy(x_vals(file) +  k_plot_scaling*eps_rr, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) - k_plot_scaling*ref_eps_rr, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_rr, r10/delta_tau), eps_rr);
  text(x_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(7);

  eps_tt = -mean(Dthth,2);

  semilogy(x_vals(file) +  k_plot_scaling*eps_tt, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) - k_plot_scaling*ref_eps_tt, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_tt, r10/delta_tau), eps_tt);
  text(x_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(8);

  eps_xx = -mean(Drr,2);

  semilogy(x_vals(file) +  20*eps_xx, r10/delta_tau, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(x_vals(file) - 20*ref_eps_zz, ref_r, 'k--');
  % labels
  text(x_vals(file)+.1,1.5, ['x = ' num2str(x_vals(file))]);
  line([x_vals(file) x_vals(file)],[.1 max(r10/delta_tau)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_zz, r10/delta_tau), eps_xx);
  text(x_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

end

hold off;
figure(1);
title('Bulk-Velocity from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(u_plot_scaling)]); 
ylabel('r^+');
legend('u_{mean}^+', 'El Khoury');%, 'visc. sublayer', 'loglaw');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);


figure(2);
title('Turbulence kinetic energy from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Energy in kg*(m/s)^2 is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('k', 'El Khoury');%, 'visc. sublayer', 'loglaw');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(3);
title('u_{x,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s)^+ is not scaled']); 
ylabel('r^+');
legend('u_{\theta,rms}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(4);
title('u_{\theta,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(rms_plot_scaling)]); 
ylabel('r^+');
legend('u_{\theta,rms}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(5);
title('u_{r,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(rms_plot_scaling)]); 
ylabel('r^+');
legend('u_{r,rms}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(6);
title('(minus) \epsilon_{rr} from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s) is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('-\epsilon_{rr}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(7);
title('(minus) \epsilon_{\theta\theta} from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s) is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('-\epsilon_{\theta\theta}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);

figure(8);
title('(minus) \epsilon_{xx} from SEM along developing straight pipe, where L = 24*R')
xlabel(['x Downstream the pipe; Velocity in (m/s) is scaled by *' num2str(20)]); 
ylabel('r^+');
legend('-\epsilon_{xx}', 'El Khoury');
axis([min(x_vals) 1.25*max(x_vals) .5 2*Re_tau]);
