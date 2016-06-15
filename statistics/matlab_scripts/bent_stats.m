%clear all;
addpath([pwd '/unsorted'])

% :-X hardcoded definitions..
Re_tau = 180; nu = 1/5300; ref_u_tau = 0.0683; ref_h12 = 1.85; ref_tauws = 0.3503;
%Re_tau = 360; nu = 1/11700; ref_u_tau = 0.0617; ref_h12 = 1.62;ref_tauws = 0.3836;

u_plot_scaling = 0.08;
rms_plot_scaling = 2;
k_plot_scaling = 100.0;

reference   = importdata(['../../../ElKhouryData/' num2str(Re_tau) '_Re_1.dat']); % El Khoury data
ref_rr_budg = importdata(['../../../ElKhouryData/' num2str(Re_tau) '_RR_Budget.dat']); % El Khoury data
ref_tt_budg = importdata(['../../../ElKhouryData/' num2str(Re_tau) '_TT_Budget.dat']); % El Khoury data
ref_zz_budg = importdata(['../../../ElKhouryData/' num2str(Re_tau) '_ZZ_Budget.dat']); % El Khoury data

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
[z_vals, ind] = sort(z_vals);
file_names = file_names(ind);

% Small function for relative error that avoids div-by-zero
%calcRelErr = @(ref_x, actual_x) norm((ref_x(ref_x>1e-10) - actual_x(ref_x>1e-10))./ref_x(ref_x>1e-10),1)/length(ref_x(ref_x>1e-10));
calcRelErr = @(ref_x, actual_x) norm(ref_x - actual_x)/norm(ref_x);

figure(1);
clf;
hold on;
figure(2);
clf;
hold on;
figure(3);
clf;
hold on;

for file = 1:length(file_names)

  fname  = ['../recordings/' file_names(file).name];
  pipe_stat;

  X(:,end+1) = X(:,1);
  Y(:,end+1) = Y(:,1);


%  contour(X,Y,Q59,'LineWidth',2)
%  contour(X,Y,sqrt(Q1.^2+Q2.^2),'LineWidth',2)

  var1 = Q57;
  var2 = Q58;
  var3 = Q59;

  var1(:,end+1) = var1(:,1);
  var2(:,end+1) = var2(:,1);
  var3(:,end+1) = var3(:,1);
  
  figure(1)

  subplot(4,3,file)
  pcolor(X,Y,var1);
  shading interp;
%  contour(X,Y,var1,'LineWidth',2)

% caxis([0, 1.3]);
% axis equal;
  title(['z = ' num2str(z_vals(file)) ', [' num2str(min(var1(:))) ', ' num2str(max(var1(:))) ']']);
  
  figure(2)

% subplot(ceil(sqrt(length(file_names))),ceil(sqrt(length(file_names))),file)
  subplot(4,3,file)
  pcolor(X,Y,var2);
  shading interp;
% caxis([0, 1.7]);
% axis equal;
  title(['z = ' num2str(z_vals(file)) ', [' num2str(min(var2(:))) ', ' num2str(max(var2(:))) ']']);
  
  figure(3)

  subplot(4,3,file)
% subplot(ceil(sqrt(length(file_names))),ceil(sqrt(length(file_names))),file)
  pcolor(X,Y,var3);
% caxis([0, 1.7]);
  shading interp;
% axis equal;
  title(['z = ' num2str(z_vals(file)) ', [' num2str(min(var3(:))) ', ' num2str(max(var3(:))) ']']);
  
  
end

