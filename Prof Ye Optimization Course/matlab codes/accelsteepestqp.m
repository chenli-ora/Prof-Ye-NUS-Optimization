%======================================================
%  Matlab implementation of the accelerated steepest 
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
if exist('beta') ~= 1 
   beta=max(eigs(Q)); 
end
if exist('maxiter') ~= 1 
   maxiter=100; 
end
x=x0;
lambda0=0;
lambda1=1;
xt1=x;
norm(Q*x+c)
for k=1:maxiter
  lambda2=(1+sqrt(1+4*lambda1^2))/2;
  alpha=(1-lambda1)/lambda2;
  g=(Q*x+c);
  xt2=x-(1/beta)*g;
  x =(1-alpha)*xt2+alpha*xt1;
  lambda1=lambda2;
  xt1=xt2;
end;
norm(Q*x+c)
% 