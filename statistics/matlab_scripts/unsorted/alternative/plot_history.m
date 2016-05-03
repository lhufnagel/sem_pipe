%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Plot time history of statistical quantities %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   averaged on the whole domain       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
close all
clear all
clc

double precision;
format long;

his_files  = 0;
grep_files = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (his_files == 1) 

sf = 1 ;
ef = 0 ;
m  = 2 ;

for i=1:m
fname = sprintf('HIS_files/history%2.2d',i)
fid   = fopen(fname);

FX  = textscan(fid, '%f %f %f %f','headerlines',1);
fclose(fid);

Rtime    = FX{1}; 
Ruzbar   = FX{2}; 
Ruzuzbar = FX{3}; 
Rscale   = FX{4}; 

ss = length(Rtime);
ef = ef + ss;

time(sf:ef,:)     = Rtime;
uzbar(sf:ef,:)   = Ruzbar;
uzuzbar(sf:ef,:) = Ruzuzbar;
scale(sf:ef,:)   = Rscale;

sf = ef +1;

end

uzuzbar_rms =  sqrt(uzuzbar - uzbar.*uzbar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf,'PaperUnits','centimeters'); 
xSize = 10; ySize = 5;
set(gcf,'Position',[600 550 xSize*50 ySize*50])

plot(time,uzbar,'k','linewidth',1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
set(gcf,'PaperUnits','centimeters'); 
xSize = 10; ySize = 5;
set(gcf,'Position',[600 550 xSize*50 ySize*50])

plot(time,uzuzbar,'k','linewidth',1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
set(gcf,'PaperUnits','centimeters'); 
xSize = 10; ySize = 5;
set(gcf,'Position',[600 550 xSize*50 ySize*50])

plot(time,scale,'k','linewidth',1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (grep_files == 1) 
    
fid = fopen('HIS_FILES/time.dat');
FX  = textscan(fid, '%q %d %s %q %f %s %q %f %s %q %f %f %f');
fclose(fid);

v1  = FX{1}; 
v2  = FX{2}; 
v3  = FX{3}; 
v4  = FX{4}; 
v6  = FX{6}; 
v7  = FX{7}; 
v9  = FX{9}; 
v10 = FX{10}; 

t   = FX{5}; 
DT  = FX{8}; 
C   = FX{11}; 
tet = FX{12}; 
tpt = FX{13}; 

% tpt: time/time-step
% tet: total elapsed time        
clear FX

fid = fopen('HIS_FILES/volflowZ.dat');
FX  = textscan(fid, '%q %d %f %f %f %f %f %q %s');
fclose(fid);

char         = FX{1};
istep        = FX{2}; 
time         = FX{3}; 
scale        = FX{4}; 
delta_flow   = FX{5}; 
current_flow = FX{6};
flow_rate    = FX{7}; 
chv          = FX{8};
chv1         = FX{9};

clear FX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 1;
A = pi*R^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf,'PaperUnits','centimeters'); 
xSize = 10; ySize = 5;
set(gcf,'Position',[600 550 xSize*50 ySize*50])

plot(time(1:end),scale(1:end),'b-','linewidth',1); 
xlabel('time')
ylabel('pressure gradient')
axis ([0 62 0 0.008])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
set(gcf,'PaperUnits','centimeters'); 
xSize = 10; ySize = 5;
set(gcf,'Position',[600 550 xSize*50 ySize*50])

plot(time,tpt,'b-','linewidth',1); 
xlabel('time')
ylabel('time/timestep')
axis ([0 62 0 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% set(gcf,'PaperUnits','centimeters'); 
% xSize = 10; ySize = 5;
% set(gcf,'Position',[600 550 xSize*50 ySize*50])
% 
% plot(time(1:end),delta_flow(1:end),'b','linewidth',1); 
% xlabel('time')
% ylabel('delta flow')
% axis ([0 61 -1 0.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% set(gcf,'PaperUnits','centimeters'); 
% xSize = 10; ySize = 5;
% set(gcf,'Position',[600 550 xSize*50 ySize*50])
% 
% plot(time,current_flow,'b','linewidth',1); 
% xlabel('time')
% ylabel('current flow')
% axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% set(gcf,'PaperUnits','centimeters'); 
% xSize = 10; ySize = 5;
% set(gcf,'Position',[600 550 xSize*50 ySize*50])
% 
% plot(time,flow_rate,'b','linewidth',1); 
% xlabel('time')
% ylabel('flow rate')
% axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end


