%======================================================
%  Matlab implementation of the steepest descent method
%  with B-B stepsize for strongly convex quadratic minim.
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
   maxiter=50; 
end
[n,m]=size(x0);
x=x0;
xx=0.1*randn(n,1);
gg=Q*xx+c;;
norm(Q*x+c)
for k=1:maxiter,
  g=(Q*x+c);
  if norm(g) <= 1.e-12, norm(g), return 
  end; 
  deltax=x-xx;
  deltag=g-gg;
  xx=x;
  gg=g;
  alpha=(deltax'*deltag)/(deltag'*deltag);
  x=x-alpha*g;
end;
norm(Q*x+c)
% 