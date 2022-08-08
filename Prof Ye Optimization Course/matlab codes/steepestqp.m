%======================================================
%  Matlab implementation of the steepest descent method
%  for convex quadratic minimization
%  
%      minimize    0.5x'Qx+c'x
%
%  Input Data
%      Q: symmetric PSD matrix
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
norm(Q*x+c)
for k=1:maxiter,
  g=(Q*x+c);
  x=x-(1/beta)*g;
end;
norm(Q*x+c)
% 