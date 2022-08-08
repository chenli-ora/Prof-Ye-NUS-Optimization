%======================================================
%  Matlab implementation of the project steepest descent
%  method for nonnegative-constraint quadratic minimization
%  
%      minimize    0.5x'Qx+c'x   s.t. x>=0
%
%  Input Data
%      Q: symmetric PSD matrix
%      c: vector
%      x0: any initial solution
%  Output
%      x: latest iterative solution
%
%   Algorithm details can be found in Sect. 12.1 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
%
if exist('beta') ~= 1
   beta=max(eigs(Q));   
end
if exist('maxiter') ~= 1 
   maxiter=200; 
end
if exist('x0') ~= 1 
  [m,n] = size(Q);
  x0=ones(n,1);
end
x=x0;
norm(Q*x+c)
for k=1:maxiter,
  % compute the gradient
  g=(Q*x+c);
  % make steepest-descent and then projection 
  x=max(0, x-(1/beta)*g);
end;
y=max(0,g);
norm(x.*y)
% 