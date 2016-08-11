% read Mass matrix
[mass_data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('../maspod0.f00001');
if status < 0 
  error('Could not read mass matrix');
end

mass_mat = reshape(mass_data(:,:,4),[size(mass_data,1)*size(mass_data,2) 1]);
mass_mat = [mass_mat; mass_mat; mass_mat];

% some metadata
N=lr1(1);
nFpp=176;

%build index map only once
eps=1e-3;
tic

disp('Build element map...')
index_map1 = zeros(nFpp,1);
gll=(N+1)*N/2;
for e=1:nFpp
  for e2=e:nFpp
    if (index_map1(e)==0)
      for gll2=1:N^2
	if (abs(mass_data(e,gll,2)+mass_data(e2,gll2,2))<eps) 
	  if (abs(mass_data(e,gll,1)-mass_data(e2,gll2,1))<eps) 
	    index_map1(e)=e2;
	    index_map1(e2)=e;
	    break
	  end 
	end
      end
    end
  end
  if (index_map1(e)==0)
    error('Element %i not found',e);
  end
end
disp('...OK');

disp('Build point map...')
index_map2 = zeros(nFpp,N^2);
for e=1:nFpp
  for gll=1:N^2
    if (index_map2(e,gll)==0)
      for gll2=1:N^2     
	if (abs(mass_data(e,gll,2)+mass_data(index_map1(e),gll2,2))<eps) 
	  if (abs(mass_data(e,gll,1)-mass_data(index_map1(e),gll2,1))<eps) 
	    index_map2(e,gll) = gll2;
	    index_map2(index_map1(e),gll2) = gll;
	    break
	  end
	end
      end
    end
    if (index_map2(e,gll)==0)
      error('Element %i point %i not found',e,gll)
    end
  end
end
disp('...OK');
toc


snapshots = dir('../pipe_coarse0.f*');
%snapshots = dir('../../field0/pipe0.f*');
if (length(snapshots) == 0)
    error('no files found')
end

nModes = 10;
snaps = 1:length(snapshots);

mat = [];
min_time = inf;
max_time = -inf;

disp('pod_coarse');
for s=snaps
 [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['../' snapshots(s).name]);
% [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['../../field0/' snapshots(s).name]);

  disp(['read in file ' snapshots(s).name]);
  min_time = min(min_time, time);
  max_time = max(max_time, time);
  %data: nek5000 data ordered as {iel,inode,[x|y|z|u|v|w|p|T|s_i]}

  uvw = reshape(data(:,:,4:6),[size(data,1)*size(data,2)*3 1]);

  uvw = (mass_mat.^(0.5)).*uvw;
  mat = [mat, uvw];

  data_mirrored = data;
  for slab=0:floor(size(data,1)/nFpp)-1
    for z_gll=0:lr1(3)-1
    for e=1:nFpp

      data_mirrored(slab*nFpp+e,z_gll*N*N+1:(z_gll+1)*N*N,4) =  data(slab*nFpp+index_map1(e),z_gll*N*N+index_map2(e,:),4);
      data_mirrored(slab*nFpp+e,z_gll*N*N+1:(z_gll+1)*N*N,5) = -data(slab*nFpp+index_map1(e),z_gll*N*N+index_map2(e,:),5);
      data_mirrored(slab*nFpp+e,z_gll*N*N+1:(z_gll+1)*N*N,6) =  data(slab*nFpp+index_map1(e),z_gll*N*N+index_map2(e,:),6);
    end
  end
  end

  uvw = reshape(data_mirrored(:,:,4:6),[size(data,1)*size(data,2)*3 1]);

  uvw = (mass_mat.^(0.5)).*uvw;
  mat = [mat, uvw];
end

[u,s,v] = svd(mat, 'econ');
%opts.tol = 1e-8;
%opts.maxit = 150;
%[u,s,v] = svds(mat, 12,'L',opts);

%Plot 
h=figure('visible','off');
semilogy(1:length(s),diag(s));
title('Singular values');
saveas(h,'singular_values_sym','fig');

time_coeffs = s*v';
h=figure('visible','off');
plot(time_coeffs(1:4,:)');
title('Time coefficients');
legend(['1';'2';'3';'4']);
saveas(h,'time_coeffs_sym','fig');
save('time_coeffs_sym.mat','time_coeffs');

for m=1:nModes

  data_new = data;
  data_new(:,:,7:end) = zeros(size(data(:,:,7:end))); %don't set pressure etc
  data_new(:,:,4:6) = reshape((mass_mat.^(-0.5)).*u(:,m), size(data(:,:,4:6)));
% data_new(:,:,4:6) = reshape(*u(:,m), size(data(:,:,4:6)));

  writenek(['pipe_mode_sym0.f' num2str(m,'%05d')],data_new,lr1,elmap,max_time-min_time,m-1,fields,emode,wdsz,etag);
end
