function  [d2QdT2] = d2Chain(n_r,n_theta,r,ssp,DM1,DM2,Dfr,RQ1) 

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

%%% d(W)/dx2 + d(W)/dy2 =   d2(W)/dr2                        ( TERM1 )
%%%                       + (1/r)*d(W)/dr +                  ( TERM2 )
%%%                       + (1/r)^2*d2(W)/d(theta)2          ( TERM3 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TERM1(1:n_r,1:n_theta) = 0; 
TERM2(1:n_r,1:n_theta) = 0; 
TERM3(1:n_r,1:n_theta) = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:n_theta

  U1 = RQ1(ssp,j); 
  
  TERM1(:,j) =  flipdim(+DM2*U1,1);
  TERM2(:,j) =  flipdim(-DM1*U1,1);  
    
  clear U1 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_r    
  if (i~=1)
          
  U1 = RQ1(i,:);
  
  TERM3(i,:) = Dfr*(Dfr*U1');

  end    
  clear U1      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d2QdT2(1:n_r,1:n_theta) = 0;

for i = 1:n_r
    for j=1:n_theta
                        
     c1T =  1         ;    
     c2T =  1./r(i)   ; 
     c3T =  1./r(i).^2; 
        
     d2QdT2(i,j) =  c1T.*TERM1(i,j) + ...
                    c2T.*TERM2(i,j) + ...
                    c3T.*TERM3(i,j);           
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

