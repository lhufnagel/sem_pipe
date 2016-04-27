function  [dQdxj] = dQ1dxj_t(it,n_r,n_theta,r,theta,ssp,DM1,...
                            Dfr,tensRR1) 

% ----------------------------------------------------------------------- %
% This function calculates the derivaties dQ/dxj from Q on (r,theta)-grid.
% Returned values are dQ/dx, dQ/dy, dQ/dz
% This applies to the third order tensor
% ----------------------------------------------------------------------- %

dQdx_1(1:n_r,1:n_theta) = 0; 
dQdy_1(1:n_r,1:n_theta) = 0; 
dQdz_1(1:n_r,1:n_theta) = 0; 

dQdx_2(1:n_r,1:n_theta) = 0; 
dQdy_2(1:n_r,1:n_theta) = 0; 
dQdz_2(1:n_r,1:n_theta) = 0; 

for j = 1:n_theta
     
    c1 = cos(theta(j)); 
    c2 = sin(theta(j));
    
  Q = squeeze(tensRR1(it,ssp,j));
  
  dQdx_1(:,j) = flipdim(-DM1*Q',1)*c1;
   
  dQdy_1(:,j) = flipdim(-DM1*Q',1)*c2;
  

  clear Q

end

for i = 2:n_r
         
    Q = squeeze(tensRR1(it,i,:)); 
    
    dQdx_2(i,:) = Dfr*Q;
    
    dQdy_2(i,:) = Dfr*Q;
           
  clear Q
    
end

dQdxj(1:3,1:n_r,1:n_theta) = 0;
      
    for i = 1:n_r
    for j=1:n_theta
     
     c1 = -sin(theta(j))./r(i);
     c2 =  cos(theta(j))./r(i);
     
     dQdxj(1,i,j) = dQdx_1(i,j) + c1*dQdx_2(i,j);
     dQdxj(2,i,j) = dQdy_1(i,j) + c2*dQdy_2(i,j);
     dQdxj(3,i,j) = dQdz_1(i,j) +    dQdz_2(i,j);
    
    end
    end

end