%======================================================
%  Matlab demonstration of solving Fisher social 
%  optimizarion by Lagrangian Method of Multipliers
%      maximize  5ln(2x1+x2)+8ln(3x3+x4)
%      s.t.      x1 + x3 =1 (y(1))
%                x2 + x4 =1 (y(2))
%                all variables >=0
%  
%  Output  x: good allocation, y: prices of goods
%
%  Detail is in Example 3 in Sect. 14.1 of
%  L&Y, Linear and nonlinear programming, 5th edition
%======================================================%  
% Set the initial data point
b=ones(2,1); % supply quantities
A=[1 0 1 0;0 1 0 1]; % constraint matrix
w1=5;u1=[2;1]; % first buyer information
w2=8;u2=[3;1]; % second buyer information
maxiter=500;
% Set the initial equilibrium prices
y=(13/2)*ones(2,1);
X=zeros(4,maxiter); % Puchase history
%
iter=0;
for k=1:maxiter;
 iter=iter+1;
 x=zeros(4,1);
 % Update first buyer's purchse:
 [u,i]=max(u1./y);
 x(i)=w1/y(i);
 % Update second buyer's purchse:
 [u,i]=max(u2./y);
 x(i+2)=w2/y(i);
 % Update the prices
 y=max(0.01, y+(1/sqrt(iter))*(A*x-b));
 %y=((w1+w2)/(b'*y))*y;
 X(:,iter)=x;
end;
sum(X')/maxiter
y,
%