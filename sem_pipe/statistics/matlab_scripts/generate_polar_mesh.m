close all
clear all
clc

double precision;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ret = 180;

L = .5; % pipe radius
C = 3; % Grading of points towards the wall

nt = 90; %number of azimuthal/theta sweeps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Generate the grid in the radial direction %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% seems to be some criterion for appropriate radial resolution
% TODO verify, probably a bound in plus-units
if (Ret == 180) 
  drp_min = 4.5;
  drw     = drp_min/Ret;
end

if (Ret == 550)
  drp_min = 0.1507;
  drw     = drp_min/Ret;
end

if (Ret == 1000)
  drp_min =  0.151527;
  drw     = drp_min/Ret;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nr = 150:1000

  [r, dh, d2h, deta] = gridr(L, C, nr);

  r10 = L-r';
  r01 =flipdim(r10,1);

  dre = r10(1) - r10(2);

  if dre < drw
    break    
  end

  clear r r10 r01 dh d2h deta

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Generate the grid in the azimuthal direction %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta(1:nt) = 0;
dtheta = 2*pi/nt;

for i = 1:nt
    theta(i) = (i-1)*dtheta; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save polar_mesh.bin %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do it as fortran binary
fid = fopen('../polar_mesh.bin', 'w','ieee-le.l64');

% first write two 4 bytes integers
data  = [C nr nt];
eor   = length(data)*4;
count = fwrite(fid,eor ,'int32');
count = fwrite(fid,data,'int32');
count = fwrite(fid,eor ,'int32');

% then write nr reals
data  = r10(1:nr); %TODO WARUM DURCH 2?!?!?
eor   = length(data)*8;
count = fwrite(fid,eor,'int32')   ;
count = fwrite(fid,data,'float64');
count = fwrite(fid,eor,'int32')   ;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xxp(1:2,1:(nr*nt)) = 0;

npoints = 0;

for i = 1:nt
for j = 1:nr 
  	 npoints = npoints + 1;	 
 	 xxp(1,npoints) = r10(j)*cos(theta(i));
  	 xxp(2,npoints) = r10(j)*sin(theta(i));
end
end

figure(11)
set(gcf,'PaperUnits','centimeters'); 
xSize = 15; ySize = 15;
set(gcf,'Position',[500 500 xSize*50 ySize*50])

plot(xxp(1,:),xxp(2,:),'b.','markersize',1)
axis tight
daspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








