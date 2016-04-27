% Construct Chebyshev derivative matrices for
% even or odd functions
%
% Works at the moment only for odd N.

N=151;
[y,D]=chebdif(N,1);

r = rand(N,1);

% Sample odd function
xodd = y;

xodd = zeros(N,1);
xodd(1:floor(N/2)) = r(1:floor(N/2));
xodd(N:-1:floor(N/2)+2) = -r(1:floor(N/2));

% sample even function
xeven = y.^2;

xeven = zeros(N,1);
xeven(1:floor(N/2+1)) = r(1:floor(N/2+1));
xeven(N:-1:floor(N/2)+2) = r(1:floor(N/2));
xeven(floor(N/2+1))=0.;


% construct the new matrices
Deven = zeros(N,N);
Dodd = zeros(N,N);
for i=1:floor(N/2+1)
  Dodd(i,:)  = D(i,:) + D(N-i+1,:);
  Deven(i,:) = D(i,:) - D(N-i+1,:);
end
Deven(floor(N/2+1),:) = 0;
Dodd  = Dodd(1:floor(N/2+1),1:floor(N/2+1));
Deven = Deven(1:floor(N/2+1),1:floor(N/2+1));


% compute derivates of the test functions
xoddp  = D * xodd;
xoddpp = Dodd*xodd(1:floor(N/2+1));


xevenp  = D*xeven;
xevenpp = Deven*xeven(1:floor(N/2+1));


% plot the results
figure;
hold on
plot(y,xoddp)
plot(y(1:floor(N/2+1)),xoddpp,'r')

figure
hold on
plot(y,xevenp)
plot(y(1:floor(N/2+1)),xevenpp,'r')
