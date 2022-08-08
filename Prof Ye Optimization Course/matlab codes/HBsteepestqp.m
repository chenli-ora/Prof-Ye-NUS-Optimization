%======================================================
%  Matlab implementation of the Heavy-Ball steepest 
%  descent method for strongly convex quadratic minim.
%  
%      minimize    0.5x'Qx+c'x
%
%  Input Data
%      Q: symmetric PD matrix
%      c: vector
%      x0: any initial solution
%  Output
%      x: latest iterative solution
%
%   Details can be found in Sect. 8.4 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
if exist('maxiter') ~= 1 
   maxiter=100; 
end
lambdan=max(eig(Q));
lambda1=min(eig(Q));
x=x0;
xx=0*x;
beta1=sqrt(lambdan)-sqrt(lambda1);
beta2=sqrt(lambdan)+sqrt(lambda1);
norm(Q*x+c)
for k=1:maxiter,
  g=(Q*x+c);
  d=-(4/beta2^2)*g+(beta1/beta2)*(x-xx);
  xx=x;
  x=x+d;
end;
norm(Q*x+c)
% 