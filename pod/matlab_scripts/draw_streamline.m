
[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] =readnek('slab_mode_2d0.f00082'); % -2
%[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] =readnek('slab_mode_2d0.f00105'); % -5.33

x = reshape(data(:,:,1),size(data,1)*size(data,2),1);
y = reshape(data(:,:,2),size(data,1)*size(data,2),1);

u = reshape(data(:,:,3),size(data,1)*size(data,2),1);
v = reshape(data(:,:,4),size(data,1)*size(data,2),1);


pts = [x,y];
vel = [u,v];

[pts, idx, ~] = unique(pts,'rows');
vel = vel(idx,:);

x = pts(:,1);
y = pts(:,2);

u = vel(:,1);
v = vel(:,2);

[xq,yq] = meshgrid(-.5:.001:.5, -.5:.001:.5);
uq = griddata(x,y,u,xq,yq);
vq = griddata(x,y,v,xq,yq);

clf
h=streamslice(xq,yq,uq,vq,1,'linear');
for i=1:length(h)
  h(i).LineWidth = 2;
  h(i).Color = [0 0 0];
end
axis equal; 
axis off;
