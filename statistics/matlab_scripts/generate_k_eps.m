ref =     importdata('180_Re_1.dat');
rr_budg = importdata('180_RR_Budget.dat'); % El Khoury data
tt_budg = importdata('180_TT_Budget.dat'); % El Khoury data
zz_budg = importdata('180_ZZ_Budget.dat'); % El Khoury data

radius = 0.5;
% Remember: George's Data was generated on Pipe with radius = 1
nu = 1/5300 * (2.*radius); % Assuming U_bulk = 1

%u_tau = 6.1746647715901784E-02;%Re_tau = 361.21...
u_tau = 6.8322301823104711E-02; %Re_tau = 181.05..

lstar = nu/u_tau;

r = ref.data(:,1);

if (r ~= rr_budg.data(:,1) | ...
  r ~= zz_budg.data(:,1) | ...
  r ~= tt_budg.data(:,1))
    disp('Input data radii must agree')
    return
end

% Turn in to real radius, instead of (1-r). 
% Rescale to [0,0.5]
r = radius*(1-r);
r(1)=0; r(end)=radius+1e-4; 
% compensate for floating point precision of mesh radius
Umean = ref.data(:,4); 

ref_k = u_tau^2*.5*(ref.data(:,5).^2+ref.data(:,6).^2+ref.data(:,7).^2);
ref_eps = - u_tau^4/nu*.5*(rr_budg.data(:,8) + zz_budg.data(:,8) + tt_budg.data(:,8));

if (length(ref_k) ~= length(Umean) | ...
  length(Umean) ~= length(ref_eps) | ...
  length(ref_eps) ~= length(r))
    disp('Input data must be of same length')
    return
end


fileID = fopen('sem_input.txt','w');
fprintf(fileID, '#r   U_mean  k   eps\n',r,Umean);
fprintf(fileID, '%d\n',length(r));
fprintf(fileID, '%1.16E   %1.16E   %1.16E   %1.16E\n',[r,Umean,ref_k,ref_eps]');
fclose(fileID);


ref_r = ref.data(:,2); % Plus units!
sigma=ref_k.^(3./2.)./ref_eps;
[sigma_max,maxi] = max(sigma);
sigma_max

clf;
subplot(2,1,1);
semilogx(ref_r,sigma);
hold on;

semilogx(ref_r,min(lstar*ref_r,min(sigma,0.41*radius)));
semilogx(ref_r,ref_eps);
semilogx(ref_r, sqrt(2/3*ref_k));
l = legend('\sigma = k^{3/2}/\epsilon', 'min(\sigma, 0.41*\delta)', '\epsilon', 'sqrt(2/3 k)');
l.Location='northwest';
plot(ref_r(maxi),sigma_max,'o' );
text(ref_r(maxi)*1.05 ,sigma_max*1.05, ['\sigma_{max} = ' num2str(sigma_max)]);

xlabel('r^+'); 
ylabel('[SI]');
grid on;

hold off;

subplot(2,1,2);
hold on;
plot(ref_r,sigma);
plot(ref_r,min(ref_r*lstar,min(sigma,0.41*radius)));
plot(ref_r,ref_eps);
plot(ref_r, sqrt(2/3*ref_k));

legend({'\sigma = k^{3/2}/\epsilon', 'min(\sigma, 0.41*\delta)', '\epsilon', 'sqrt(2/3 k)'},'Position',[0.75,0.2,.1,.05]);
plot(ref_r(maxi),sigma_max,'o' );
text(ref_r(maxi)*1.05 ,sigma_max*1.05, ['\sigma_{max} = ' num2str(sigma_max)]);
xlabel('r^+'); 
ylabel('[SI]');
grid on;
%
hold off;

