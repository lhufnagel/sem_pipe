function  [d2Qdx2,d2Qdy2] = d2Chain(n_r,n_theta,r,theta,ssp,DM1,DM2,Dfr,RQ1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% d(W)/dx2 = + dW/dr*(sin(theta)^2)/r                       ( TERM1 )
%%%            + cos(theta)^2*d2(W)/dr2                       ( TERM2 )
%%%            - 2*cos(theta)*sin(theta)/r*d2(W)/drd(theta)   ( TERM3 )
%%%            + 2*cos(theta)*sin(theta)/r^2*d(W)/d(theta)    ( TERM4 )
%%%            + sin(theta)^2/r^2*d2W/d(theta)2               ( TERM5 ) 


%%% d(W)/dy2 = + dW/dr*(cos(theta)^2)/r                      ( TERM1 )
%%%            + sin(theta)^2*d2(W)/dr2                      ( TERM2 )
%%%            + 2*cos(theta)*sin(theta)/r*d2(W)/drd(theta)  ( TERM3 )
%%%            - 2/r^2*sin(theta)*cos(theta)*dW/d(theta)     ( TERM4 )
%%%            + cos(theta)^2/r^2*d2W/d(theta)2              ( TERM5 ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TERM1(1:n_r,1:n_theta) = 0; 
TERM2(1:n_r,1:n_theta) = 0; 
TERM3(1:n_r,1:n_theta) = 0; 

TERM4(1:n_r,1:n_theta) = 0; 
TERM5(1:n_r,1:n_theta) = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:n_theta

  U1 = RQ1(ssp,j);     
  
  TERM1(:,j) =  flipdim(-DM1*U1,1);  
  TERM2(:,j) =  flipdim(+DM2*U1,1);
    
  clear U1 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_r
     
    if (i~=1)
          
  U1 = RQ1(i,:);

  TERM3(i,:) = Dfr*TERM1(i,:)';
  TERM4(i,:) = Dfr*U1';
  TERM5(i,:) = Dfr*(Dfr*U1');

    end
    
  clear U1  
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d2Qdx2(1:n_r,1:n_theta) = 0;


for i = 1:n_r
    for j=1:n_theta
        
     c1x = 1./r(i);           
     c2x = 1; 
     c3x = 0; 
     c4x = 0;
     c5x = 1./(r(i).^2); 
     
     d2Qdx2(i,j) =  c1x.*TERM1(i,j) + ...
                    c2x.*TERM2(i,j) + ...
                    c3x.*TERM3(i,j) + ...
                    c4x.*TERM4(i,j) + ...
                    c5x.*TERM5(i,j);                
    end
end

% 
% for i = 1:n_r
%     for j=1:n_theta
%         
%      c1x = 0;%sin(theta(j)).^2./r(i);           
%      c2x = 0;%cos(theta(j)).^2; 
%      c3x = -2*sin(theta(j)).*cos(theta(j))./r(i); 
%      c4x = 0;%2*sin(theta(j)).*cos(theta(j))./(r(i).^2);
%      c5x = 0;%sin(theta(j)).^2/(r(i).^2); 
%      
%      d2Qdx2(i,j) =  c1x.*TERM1(i,j) + ...
%                     c2x.*TERM2(i,j) + ...
%                     c3x.*TERM3(i,j) + ...
%                     c4x.*TERM4(i,j) + ...
%                     c5x.*TERM5(i,j);                
%     end
% end

d2Qdy2(1:n_r,1:n_theta) = 0;

% for i = 1:n_r
%     for j=1:n_theta
%                         
%      c1y = 0;%cos(theta(j)).^2./r(i);    
%      c2y = 0;%sin(theta(j)).^2; 
%      c3y = 2*sin(theta(j)).*cos(theta(j))./r(i); 
%      c4y = 0;%-2*sin(theta(j)).*cos(theta(j))./(r(i).^2);
%      c5y = 0;%cos(theta(j)).^2./r(i).^2;
%       
%      d2Qdy2(i,j) =  c1y.*TERM1(i,j) + ...
%                     c2y.*TERM2(i,j) + ...
%                     c3y.*TERM3(i,j) + ...
%                     c4y.*TERM4(i,j) + ...
%                     c5y.*TERM5(i,j);              
%     
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
