% read Mass matrix
[mass_data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('../maspod0.f00001');
if status < 0 
  error('Could not read mass matrix');
end

mass_mat = reshape(mass_data(:,:,4),[size(mass_data,1)*size(mass_data,2) 1]);
mass_mat = [mass_mat; mass_mat; mass_mat];


snapshots = dir('../pipe_coarse0.f*');
if (length(snapshots) == 0)
    error('no files found')
end

nModes = 10;
snaps = 1:length(snapshots);

mat = [];
min_time = inf;
max_time = -inf;

for s=snaps
  [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['../' snapshots(s).name]);
  min_time = min(min_time, time);
  max_time = max(max_time, time);
  %data: nek5000 data ordered as {iel,inode,[x|y|z|u|v|w|p|T|s_i]}

  uvw = reshape(data(:,:,4:6),[size(data,1)*size(data,2)*3 1]);

  uvw = (mass_mat.^(0.5)).*uvw;
  mat = [mat, uvw];
end

[u,s,v] = svd(mat, 'econ');

%Plot 
h=figure('visible','off');
semilogy(1:length(s),diag(s));
title('Singular values');
saveas(h,'singular_values','fig');

time_coeffs = s*v';
h=figure('visible','off');
plot(time_coeffs(1:4,:)');
title('Time coefficients');
legend(['1';'2';'3';'4']);
saveas(h,'time_coeffs','fig');

for m=1:nModes

  data_new = data;
  data_new(:,:,7:end) = zeros(size(data(:,:,7:end))); %don't set pressure etc
  data_new(:,:,4:6) = reshape((mass_mat.^(-0.5)).*u(:,m), size(data(:,:,4:6)));

  writenek(['pipe_mode0.f' num2str(m,'%05d')],data_new,lr1,elmap,max_time-min_time,m-1,fields,emode,wdsz,etag);
end
