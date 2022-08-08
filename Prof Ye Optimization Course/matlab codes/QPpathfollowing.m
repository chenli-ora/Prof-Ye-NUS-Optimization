%======================================================
%  Matlab demonstration of the path-following method
%  for non-strongly convex quadratic minimization
%  (assume that an KKT solution exists)
%  
%      minimize    0.5x'Qx+c'x
%
%  Input Data
%      Q: symmetric PSD matrix
%      c: vector
%  Output
%      x: latest iterative solution
%
%   Details can be found in Sect. 8.7 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
%
[m,n]=size(Q);
x=0.01*rand(n,1);
mu=10;
%
for k=1:20,
g=Q*x+c;
d=(Q+mu*speye(n))\(g+mu*x);
x=x-d;
mu=mu/2;
end
norm(Q*x+c)
% 