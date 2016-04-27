function  [RR1]  = cyl_Vel_R(R1,theta,ng_r,ng_theta) 

% ----------------------------------------------------------------------- %
% This function evalutes the radial, circumferential and axial velocities %
% from U, V, W                                                            % 
% ----------------------------------------------------------------------- %

RR1(1:3,1:ng_r,1:ng_theta) = 0;

for j = 1:ng_theta
    for i = 1:ng_r
        
        Rot_m   = [ cos(theta(j)) sin(theta(j)) 0 ; ...
                   -sin(theta(j)) cos(theta(j)) 0 ; ...
                          0           0         1];

                      
% Rot_m   = [  1        0             0       ; ...
%              0  cos(theta(j)) sin(theta(j)) ; ...
%              0 -sin(theta(j)) cos(theta(j))];
% 

        RR1(:,i,j) = Rot_m*R1(:,i,j);
  
        
        clear  Rot_m
    end
end
