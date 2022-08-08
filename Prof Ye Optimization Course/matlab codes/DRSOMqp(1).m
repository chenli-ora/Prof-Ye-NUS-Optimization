%======================================================
%  Matlab implementation of the Dimension-Reduced Hessian 
%  2nd-order descent method for convex quadratic minim.
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
%   Details can be found in Sect. 9.7 of
%   L&Y, Linear and nonlinear programming, 5th edition
%   Lecture Note #10, Slides 2-3
%======================================================% 
if exist('beta') ~= 1 
   beta=max(eigs(Q)); 
end
if exist('maxiter') ~= 1 
   maxiter=50; 
end
x=x0;
obhis=0.5*(x0'*Q*x0)+c'*x0;
%Make a steepest-descent step
g=Q*x+c;
norm(g)
d=-g/beta;
x=x+d;
obhis=[obhis 0.5*(x'*Q*x)+c'*x];
%Start the process
for k=1:maxiter,
  g=(Q*x+c);
  gnorm=norm(g);
  gTd=g'*d;
  % the first-order direction using two descent directions
  %dnorm=norm(d);
  %[gnorm^2 -gTd;-gTd dnorm^2]\[gnorm^2/beta;-gTd/beta];
  %
  % the second-order direction using dimension-reduced Hessian
  Qd=Q*d;
  [g'*(Q*g) -g'*Qd;-g'*Qd d'*Qd]\[gnorm^2;-gTd];
  %
  alphag=ans(1);alpham=ans(2);
  d=-alphag*g+alpham*d;
  x=x+d;
  obhis=[obhis 0.5*(x'*Q*x)+c'*x];
end;
norm(Q*x+c)
% 