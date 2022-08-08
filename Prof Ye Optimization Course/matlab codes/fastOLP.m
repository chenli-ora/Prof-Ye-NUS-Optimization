%======================================================
%  Matlab Demonstration of online LP of resource allocation
%      max   c'*x
%      s.t.  Ax <= b, 0 <= x <= 1.
%
%  Input 
%      A: inequality constraint matrix.
%      b(>0): inequality right-hand column vector
%      c: objective coefficient vector
%  Output
%          x: production level
%          p: dual prices of resources
%
%  Algorithm: dual stochastic sug-gradient method 
%             Sect. 8.8.
%  Problem can be found from Sect. 3.5 of Chap. 3 of 
%  L&Y, Linear and nonlinear programming, 5th edition
%======================================================%
function [x,p] = fastOLP(A,c,b);
% Set parameters
%
 [m,n] = size(A);
 br= b; % Set the initial remaining resources
 dk=b/n; % Set the avarage resouce inventory
 step=1/sqrt(n); % Set the step-size
%
 x = zeros(n,1); % Set initial primal solutions
 p = zeros(m,1); % Set initial prices
% p = (sum(c)/(m*sum(b)))*ones(m,1); % Set initial shadow prices
%
rp=randperm(n); %Randomly permute variable-order
%
% Start the loop
  for i=1:n,
   ii=rp(i);
   aa=A(:,ii);
   %
   % Set the primal increment
   xk = (sign(c(ii)-aa'*p)+1)/2;
   %
   % Update the dual soluton
   p=max(0,p+step*(xk*aa-dk)); 
   %
   % Update the remaining inventory and primal solution
   if min(br-xk*aa) >= 0,
       br=br-xk*aa;
       x(ii)=xk;
   end
   %
  end;
%