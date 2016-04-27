function  [Deven,Dodd,y] = DM_oddeven(N) 

% ----------------------------------------------------------------------- %
% Construct Chebyshev derivative matrices for even or odd functions
% Works for odd N
% ----------------------------------------------------------------------- %

 [y,DS]=chebdif(N,1);

 Deven(1:floor(N/2+1),:) = [ DS(1:floor(N/2+1),1:floor(N/2+1)-1) + ... 
                             DS(1:floor(N/2+1),N:-1:floor(N/2+1)+1), ...  
                             DS(1:floor(N/2+1),floor(N/2 + 1))];
 
 Dodd(1:floor(N/2+1),:) = [ DS(1:floor(N/2+1),1:floor(N/2+1)-1) - ... 
                             DS(1:floor(N/2+1),N:-1:floor(N/2+1)+1), ...  
                             DS(1:floor(N/2+1),floor(N/2 + 1))];
                                     