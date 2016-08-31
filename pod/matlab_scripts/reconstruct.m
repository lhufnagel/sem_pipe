snapshots = 50:500;
used_modes = [1:3]; % 1 is mean!

modes = dir('./pipe_mode_sym0.f*');
if (length(modes) == 0)
    error('no files found')
end

if ~exist('time_coeffs')
  load('time_coeffs_sym.mat');
end

mode=cell(size(used_modes));

for m=used_modes
  [mode{m},lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(['./' modes(m).name]);
  %data: nek5000 data ordered as {iel,inode,[x|y|z|u|v|w|p|T|s_i]}
end

for s =snapshots
  disp(['snapshot ' num2str(s)]);

  data_new = mode{used_modes(1)};
  data_new(:,:,4:end) = time_coeffs(used_modes(1),(s-1)*2+1) * mode{used_modes(1)}(:,:,4:end);

  for m=used_modes(2:end)
    data_new(:,:,4:end) = data_new(:,:,4:end) + time_coeffs(used_modes(m),(s-1)*2+1) * mode{used_modes(m)}(:,:,4:end);
  end

  writenek(['reconstruction/reconstruction_sym0.f' num2str(s,'%05d')],data_new,lr1,elmap,s*0.25,s,fields,emode,wdsz,etag);
end

