clear all 
close all

disp('Reading data...')

[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('./pipe_coarse0.f00680');
if status < 0 
  error('Could not read file');
else
  disp('...OK');
end
%data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])

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
	if (abs(data(e,gll,2)+data(e2,gll2,2))<eps) 
	  if (abs(data(e,gll,1)-data(e2,gll2,1))<eps) 
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
	if (abs(data(e,gll,2)+data(index_map1(e),gll2,2))<eps) 
	  if (abs(data(e,gll,1)-data(index_map1(e),gll2,1))<eps) 
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


snapshots = dir('./pipe_coarse0.f*');
if (length(snapshots) == 0)
    error('no files found')
end

for s=1:length(snapshots)
  [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['./' snapshots(s).name]);

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

  writenek(['mirror0.f' num2str(s,'%05d')],data_mirrored,lr1,elmap,time,istep,fields,emode,wdsz,etag);
end
