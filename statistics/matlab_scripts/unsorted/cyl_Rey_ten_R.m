function  [RR2] = cyl_Rey_ten_R(RQ2,theta,n_r,n_theta) 

% ----------------------------------------------------------------------- %
% This function evalutes the Reynolds stress tensor in cylindrical        %
% coordinates                                                             % 
% ----------------------------------------------------------------------- %

RR2(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
        
 %      Rot_m   = [ cos(theta(j)) sin(theta(j)) 0;
 %                -sin(theta(j)) cos(theta(j)) 0;
 %                      0           0          1];


 Rot_m   = [  1        0             0       ; ...
              0  cos(theta(j)) sin(theta(j)) ; ...
              0 -sin(theta(j)) cos(theta(j))];
 
        
        RR2(:,:,i,j) = Rot_m*RQ2(:,:,i,j)*Rot_m';
        
        clear Rot_m 
    end
end
