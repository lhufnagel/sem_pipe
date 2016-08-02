
[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('../../../field0/pipe0.f00001');
if status < 0 
  error('Could not read file');
end

x_val = -19.85 + 5*pi/2; % approx 12 diameters after exit
%x_val = -9.85 + 5*pi/2; % approx 2 diameters after exit
nFpp = 960;

e = (abs(data(:,1,1) - x_val)<1e-5); %build element index set only once

data2d = zeros(nFpp,lr1(1)*lr1(2),6);

%data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])

snapshots = dir('../../../field0/pipe0.f*');
if (length(snapshots) == 0)
    error('no files found')
end

for s=691:length(snapshots)
  [data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['../../../field0/' snapshots(s).name]);

  data2d(:,:,1) = data(e,1:(lr1(1)*lr1(2)),3); %x
  data2d(:,:,2) = data(e,1:(lr1(1)*lr1(2)),2); %y
  data2d(:,:,3) = data(e,1:(lr1(1)*lr1(2)),6); %u
  data2d(:,:,4) = data(e,1:(lr1(1)*lr1(2)),5); %v
  data2d(:,:,5) = -data(e,1:(lr1(1)*lr1(2)),4); %abuse p for w

  % center x axis
  data2d(:,:,1) = data2d(:,:,1) - .5*(max(max(data2d(:,:,1)))+min(min(data2d(:,:,1))));

  lr1(3)=1;
  elmap = int32([1:nFpp]);


  writenek(['slab12d0.f' num2str(s,'%05d')],data2d,lr1,elmap,time,istep,fields,emode,wdsz,etag);
end
