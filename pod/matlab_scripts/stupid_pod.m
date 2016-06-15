if ~exist('fname')
  fname = '../new0.f00001';
end

%Glob files
snapshots = dir('../pipe_coarse0.f*');
if (length(snapshots) == 0)
    error('no files found')
end

nModes = 10;
snaps = 1:2:100;
if ~exist('nElPerFace')
  nElPerFace = 52; % TAKE THIS FROM ../../pipe.upar ..
end

mat = [];

for s=snaps
  [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['../' snapshots(s).name]);
  %data: nek5000 data ordered as {iel,inode,[x|y|z|u|v|w|p|T|s_i]}

  if mod(size(data,1),nElPerFace) > 0
    error('Number of Elements per face wrong!')
  end

  mat = [mat, reshape(data(:,:,4:6),[size(data,1)*size(data,2)*3 1])];


  %put u,v,w into matrix: 
  %  [u;v;w] of one pipe crossection form a column

  %Γ u_1       u_(lx1*ly1+1) .. |
  %| u_2       u_(lx1*ly1+2) .. |
  %| ..        ..            .. |
  %| u_lx1*ly1 u_(2*lx1*ly1) .. |
  %| v_1       v_(lx1*ly1+1) .. |
  %| v_2       v_(lx1*ly1+2) .. |
  %| ..        ..            .. |
  %| v_lx1*ly1 v_(2*lx1*ly1) .. |
  %| w_1       w_(lx1*ly1+1) .. |
  %| w_2       w_(lx1*ly1+2) .. |
  %| ..        ..            .. |
  %| v_lx1*ly1 w_(2*lx1*ly1) .. |
  % (below follows second element etc. until nElPerFace)
  %| ..        ..            .. |

  %lx1 = lr1(1);
  %ly1 = lr1(2);
  %lz1 = lr1(3);
  
  %nEl = size(data,1);
  %nNodes = size(data,2)/lz1;
  %nLin = 3*nElPerFace*nNodes;
  %nCol = nEl/nElPerFace*lz1;
  
  %mat = zeros(nLin,nCol);


  %for i=(nEl/nElPerFace) % Elements streamwise
  %  for j=1:lz1
  %    col =[];
  %    for k=1:nElPerFace
  %      el_nodes = [];
  %      el_nodes(         1 :   nNodes) = data((i-1)*nElPerFace+k,(j-1)*nNodes+1:j*nNodes,4); % u
  %      el_nodes(  nNodes+1 : 2*nNodes) = data((i-1)*nElPerFace+k,(j-1)*nNodes+1:j*nNodes,5); % v
  %      el_nodes(2*nNodes+1 : 3*nNodes) = data((i-1)*nElPerFace+k,(j-1)*nNodes+1:j*nNodes,6); % w
  %      col = [col;el_nodes'];
  %    end
  %    mat(:,(i-1)*lz1+j) = col';
  %  end
  %end

end

[u,s,v] = svds(mat, nModes);

%Plot 
semilogy(1:nModes,diag(s).^2/sum(diag(s).^2));
title('Relative energy content of mode');

for m=1:nModes

  data_new = data;
  data_new(:,:,4:6) = reshape(u(:,m),size(data(:,:,4:6)));

% col == u(:,m)

%data_new = data; data_new(:,:,4:end) = 0.;
%for k=1:nElPerFace
%  data_new(k,1:nNodes,4) = u( 3*(k-1)*nNodes+1 : (3*k-2)*nNodes,m); % u
%  data_new(k,1:nNodes,5) = u( (3*k-2)*nNodes+1 : (3*k-1)*nNodes,m); % v
%  data_new(k,1:nNodes,6) = u( (3*k-1)*nNodes+1 :  3*k*nNodes,m); % w
%end

  writenek(['pipe_mode0.f' num2str(m,'%05d')],data_new,lr1,elmap,time,istep,fields,emode,wdsz,etag);
end
