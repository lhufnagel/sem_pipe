function  [dUidxj] = dUidxj_t(n_r,n_theta,r,theta,ssp,DM1,Dfr,RQ1) 

% ----------------------------------------------------------------------- %
% This function calculates the derivaties dUi/dxj from Ui on (r,theta)-grid
% ----------------------------------------------------------------------- %

dUdx_1(1:n_r,1:n_theta) = 0; 
dUdy_1(1:n_r,1:n_theta) = 0; 
dUdz_1(1:n_r,1:n_theta) = 0; 

dVdx_1(1:n_r,1:n_theta) = 0; 
dVdy_1(1:n_r,1:n_theta) = 0; 
dVdz_1(1:n_r,1:n_theta) = 0; 

dWdx_1(1:n_r,1:n_theta) = 0; 
dWdy_1(1:n_r,1:n_theta) = 0; 
dWdz_1(1:n_r,1:n_theta) = 0; 

dUdx_2(1:n_r,1:n_theta) = 0; 
dUdy_2(1:n_r,1:n_theta) = 0; 
dUdz_2(1:n_r,1:n_theta) = 0; 

dVdx_2(1:n_r,1:n_theta) = 0; 
dVdy_2(1:n_r,1:n_theta) = 0; 
dVdz_2(1:n_r,1:n_theta) = 0; 

dWdx_2(1:n_r,1:n_theta) = 0; 
dWdy_2(1:n_r,1:n_theta) = 0; 
dWdz_2(1:n_r,1:n_theta) = 0; 

for j = 1:n_theta
     
    c1 = cos(theta(j)); 
    c2 = sin(theta(j));
        
  U1 = squeeze(RQ1(1,ssp,j))';
  V1 = squeeze(RQ1(2,ssp,j))';
  W1 = squeeze(RQ1(3,ssp,j))';
  
  dUdx_1(:,j) = flipdim(-DM1*U1,1)*c1;
  dUdy_1(:,j) = flipdim(-DM1*U1,1)*c2;
  
  dVdx_1(:,j) = flipdim(-DM1*V1,1)*c1;
  dVdy_1(:,j) = flipdim(-DM1*V1,1)*c2;
  
  dWdx_1(:,j) = flipdim(-DM1*W1,1)*c1;
  dWdy_1(:,j) = flipdim(-DM1*W1,1)*c2;
     
  clear U1 V1 W1

end

for i = 1:n_r
     
    if (i~=1)
          
  U1 = squeeze(RQ1(1,i,:));
  V1 = squeeze(RQ1(2,i,:));
  W1 = squeeze(RQ1(3,i,:));
      
  dUdx_2(i,:) = Dfr*U1;
  dUdy_2(i,:) = Dfr*U1;
  
  dVdx_2(i,:) = Dfr*V1;
  dVdy_2(i,:) = Dfr*V1;
  
  dWdx_2(i,:) = Dfr*W1;
  dWdy_2(i,:) = Dfr*W1;

    end
    
  clear U1 V1 W1    
    
end

dUdx(1:n_r,1:n_theta) = 0; 
dUdy(1:n_r,1:n_theta) = 0; 
dUdz(1:n_r,1:n_theta) = 0; 

dVdx(1:n_r,1:n_theta) = 0; 
dVdy(1:n_r,1:n_theta) = 0; 
dVdz(1:n_r,1:n_theta) = 0; 

dWdx(1:n_r,1:n_theta) = 0; 
dWdy(1:n_r,1:n_theta) = 0; 
dWdz(1:n_r,1:n_theta) = 0; 

for i = 1:n_r
    for j=1:n_theta
     
     c1 = -sin(theta(j))./r(i); 
     c2 =  cos(theta(j))./r(i);
        
dUdx(i,j) = dUdx_1(i,j) + c1*dUdx_2(i,j);
dUdy(i,j) = dUdy_1(i,j) + c2*dUdy_2(i,j);
dUdz(i,j) = dUdz_1(i,j) +    dUdz_2(i,j);

dVdx(i,j) = dVdx_1(i,j) + c1*dVdx_2(i,j);
dVdy(i,j) = dVdy_1(i,j) + c2*dVdy_2(i,j);
dVdz(i,j) = dVdz_1(i,j) +    dVdz_2(i,j);

dWdx(i,j) = dWdx_1(i,j) + c1*dWdx_2(i,j);
dWdy(i,j) = dWdy_1(i,j) + c2*dWdy_2(i,j);
dWdz(i,j) = dWdz_1(i,j) +    dWdz_2(i,j);
    
    end
end

dUidxj(1:3,1:3,1:n_r,1:n_theta) = 0;

for j = 1:n_theta
   for i = 1:n_r
         dUidxj(:,:,i,j) = [dUdx(i,j) dUdy(i,j) dUdz(i,j);
                            dVdx(i,j) dVdy(i,j) dVdz(i,j); 
                            dWdx(i,j) dWdy(i,j) dWdz(i,j)];
    end
end    


end