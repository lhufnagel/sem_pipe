%clear all;
addpath([pwd '/unsorted'])

% :-X hardcoded definitions..
Re_tau = 180;
%Re_tau = 1.8105409983122749E+02; % George's Reference data at this Re_tau
ref_u_tau = 6.8322301823104711E-02;
nu = 1/5300;

u_plot_scaling = 0.08;
rms_plot_scaling = 2;
k_plot_scaling = 100.0;


%reference = importdata('/scratch/hufnagel/MSc/ElKhouryData/180_Re_1.dat'); % El Khoury data
reference = importdata('../../../ElKhouryData/180_Re_1.dat'); % El Khoury data
ref_rr_budg = importdata('../../../ElKhouryData/180_RR_Budget.dat'); % El Khoury data
ref_tt_budg = importdata('../../../ElKhouryData/180_TT_Budget.dat'); % El Khoury data
ref_zz_budg = importdata('../../../ElKhouryData/180_ZZ_Budget.dat'); % El Khoury data

ref_eps_rr = ref_u_tau^4/nu*ref_rr_budg.data(:,8);
ref_eps_tt = ref_u_tau^4/nu*ref_tt_budg.data(:,8);
ref_eps_zz = ref_u_tau^4/nu*ref_zz_budg.data(:,8);

ref_r = reference.data(:,2);
ref_u = reference.data(:,3);
ref_ur_rms = reference.data(:,5);
ref_ut_rms = reference.data(:,6);
ref_uz_rms = reference.data(:,7);
ref_uzur = reference.data(:,8);
ref_k = ref_u_tau^2*.5*(ref_ur_rms.^2+ref_ut_rms.^2+ref_uz_rms.^2);

%Glob files
file_names = dir('../recordings/polar_z_*');
z_vals = zeros(size(file_names));
if (length(z_vals) == 0)
    error('no files found')
end
for i =1:length(file_names)
  z_vals(i)=str2num(file_names(i).name(9:end));
end

% Small function for relative error that avoids div-by-zero
%calcRelErr = @(ref_x, actual_x) norm((ref_x(ref_x>1e-10) - actual_x(ref_x>1e-10))./ref_x(ref_x>1e-10),1)/length(ref_x(ref_x>1e-10));
calcRelErr = @(ref_x, actual_x) norm(ref_x - actual_x)/norm(ref_x);

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
figure(9);
clf;
figure(10);
clf;
utaus = [];
tauws = [];
h12s = [];


for file = 1:length(file_names)

  fname  = ['../recordings/' file_names(file).name];
  pipe_stat;

  figure(1);
  u_mean = R1(:,3);

  %visc = ref_r;
  %visc = visc(visc < 30);
  %loglaw = ref_r;
  %loglaw = loglaw(loglaw > 5);

  semilogy(z_vals(file) + u_plot_scaling*u_mean/u_tau, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + u_plot_scaling*ref_u, ref_r, 'k--');

  % Plot viscous sublayer
  %semilogy(z_vals(file) +  u_plot_scaling*visc,visc,'b--');

  %% Plot log-law
  %semilogy(z_vals(file)*ones(size(loglaw)) + u_plot_scaling*(1/0.41*log(loglaw)+5.0),loglaw,'r--');

  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, ref_u, r10/lstar,'pchip'), u_mean/u_tau);
  text(z_vals(file)+.5*max(u_mean),max(ref_r)*.01+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);


  figure(2);

  %take the abs() because unfortunately some values are <0
  ur_rms = sqrt(abs(squeeze(R2(1,1,:))));
  ut_rms = sqrt(abs(squeeze(R2(2,2,:))));
  uz_rms = sqrt(abs(squeeze(R2(3,3,:))));

  uz_ur =  abs(squeeze(R2(1,3,:)));

  k = .5*(uz_rms.^2+ur_rms.^2+ut_rms.^2);

  semilogy(z_vals(file) +  k_plot_scaling*k, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + k_plot_scaling*ref_k, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);
  err_rel = calcRelErr(interp1(ref_r, ref_k, r10/lstar,'pchip'), k);
  text(z_vals(file)+.5*max(ref_k),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(3);

  semilogy(z_vals(file) +  uz_rms/u_tau, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + ref_uz_rms, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, ref_uz_rms, r10/lstar,'pchip'), uz_rms/u_tau);
  text(z_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(4);

  semilogy(z_vals(file) +  uz_ur/u_tau^2, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + ref_uzur, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, ref_uzur, r10/lstar,'pchip'), uz_ur/u_tau^2);
  text(z_vals(file)+.5*max(ref_uzur),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);


  figure(5);

  semilogy(z_vals(file) +  rms_plot_scaling*ut_rms/u_tau, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + rms_plot_scaling*ref_ut_rms, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, ref_ut_rms, r10/lstar,'pchip'), ut_rms/u_tau);
  text(z_vals(file)+.5*max(ref_ut_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(6);

  semilogy(z_vals(file) +  rms_plot_scaling*ur_rms/u_tau, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) + rms_plot_scaling*ref_ur_rms, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, ref_ur_rms, r10/lstar,'pchip'), ur_rms/u_tau);
  text(z_vals(file)+.5*max(ref_ur_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(7);

  eps_rr = -mean(Drr,2);

  semilogy(z_vals(file) +  k_plot_scaling*eps_rr, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) - k_plot_scaling*ref_eps_rr, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_rr, r10/lstar,'pchip'), eps_rr);
  text(z_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(8);

  eps_tt = -mean(Dtt,2);

  semilogy(z_vals(file) +  k_plot_scaling*eps_tt, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) - k_plot_scaling*ref_eps_tt, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_tt, r10/lstar,'pchip'), eps_tt);
  text(z_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  figure(9);

  eps_zz = -mean(Dzz,2);

  semilogy(z_vals(file) +  20*eps_zz, r10/lstar, 'x-');
  hold on; % Have to put hold after first semilog call - matlab bug

  % Plot El-Khoury reference
  semilogy(z_vals(file) - 20*ref_eps_zz, ref_r, 'k--');
  % labels
  text(z_vals(file)+.1,1.5, ['x = ' num2str(z_vals(file))]);
  line([z_vals(file) z_vals(file)],[.1 max(r10/lstar)]);

  err_rel = calcRelErr(interp1(ref_r, -ref_eps_zz, r10/lstar,'pchip'), eps_zz);
  text(z_vals(file)+.5*max(ref_uz_rms),4+mod(file,2) , ['||e_{rel}|| = ' num2str(err_rel)],'Color','red','FontSize',12);

  tauws = [tauws tauw_rms];
  h12s = [h12s H12];
  utaus = [utaus u_tau];
end

hold off;
figure(1);
title('Bulk-Velocity from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(u_plot_scaling)]); 
ylabel('r^+');
legend('u_{mean}^+', 'El Khoury');%, 'visc. sublayer', 'loglaw');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);


figure(2);
title('Turbulence kinetic energy from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Energy in kg*(m/s)^2 is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('k', 'El Khoury');%, 'visc. sublayer', 'loglaw');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(3);
title('u_{z,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s)^+ is not scaled']); 
ylabel('r^+');
legend('u_{\theta,rms}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(4);
title('<u_zu_r>^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s)^+ is not scaled']); 
ylabel('r^+');
legend('<u_zu_r>^+', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(5);
title('u_{\theta,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(rms_plot_scaling)]); 
ylabel('r^+');
legend('u_{\theta,rms}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(6);
title('u_{r,rms}^+ from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s)^+ is scaled by *' num2str(rms_plot_scaling)]); 
ylabel('r^+');
legend('u_{r,rms}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(7);
title('(minus) \epsilon_{rr} from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s) is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('-\epsilon_{rr}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(8);
title('(minus) \epsilon_{\theta\theta} from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s) is scaled by *' num2str(k_plot_scaling)]); 
ylabel('r^+');
legend('-\epsilon_{\theta\theta}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(9);
title('(minus) \epsilon_{zz} from SEM along developing straight pipe, where L = 24*R')
xlabel(['z in diameters downstream the pipe; Velocity in (m/s) is scaled by *' num2str(20)]); 
ylabel('r^+');
legend('-\epsilon_{zz}', 'El Khoury');
axis([min(z_vals) 1.25*max(z_vals) .5 2*Re_tau]);

figure(10);
subplot(3,1,1)
plot(z_vals,tauws, 'x-');
grid on;
legend('\tau_{w,rms}');
subplot(3,1,2)
plot(z_vals,h12s, 'x-');
grid on;
legend('H_{12}');
hold on;
subplot(3,1,3)
plot(z_vals,utaus/ref_u_tau, 'x-');
grid on;
legend('u_{\tau}/u_{\tau,ref}');
xlabel('z in diameters downstream the pipe'); 
