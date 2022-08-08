%======================================================
%  Matlab implementation of the steepest descent method
%  with optimal stepsize for strongly convex quadratic min.
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
%   Details can be found in Sect. 8.2 of
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
  alpha=(g'*g)/(g'*(Q*g));
  x=x-alpha*g;
end;
norm(Q*x+c)
% 