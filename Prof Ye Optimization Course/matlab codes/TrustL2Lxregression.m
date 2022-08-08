%======================================================
%  Matlab demonstration of the regularized
%  linear regression with 0.5-norm and 1-norm (LASSO) 
%  regularizations (compressed sensing)
%
%      min   0.5*\|Ax-b\|^2 + mu* sum |x(j)|^{0.5}
%
%  Algorithm: interior (affine-scaling) trust region
%  Input 
%      A: Sparse constraint matrix.
%      b: the right-hand vector
%      x0:initial solution
%      mu0: initial regularization weight
%  Output
%     xh  : solution of 1/2-norm
%     x1  : solution of 1-norm
%  
%   Problem can be found homework #7.17 in Sect. 7.2 
%   and Algorithm in Sect. 8.7 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
function [xh,x1]=TrustL2Lxregression(A,b,x0,mu,maxiter)
if exist('mu') ~= 1 
   mu=0.1; 
end
if exist('maxiter') ~= 1 
   maxiter=25; 
end
[m,n]=size(A);
ATA=A'*A;
ee=ones(n,1);
if exist('x0') ~= 1 
   x0=ones(n,1); 
end
% Compute L2Lhalf solution
x=x0;
lambda=10;
iter=0;
%  Repeatly solving the affine-scalling trust-region step 
while (iter < maxiter),
   iter=iter+1;
   % generate affine-scaled gradient and Hessian and 
   % apply the trust-region method
   dd = abs(x);
   DD = diag(dd);
   gg = dd.*(A'*(A*x-b))+(mu/2)*(sqrt(dd).*sign(x));
   HH = DD*ATA*DD-(mu/4)*diag(sqrt(dd));
   emin=min(eig(HH));
   if (emin > 0), emin=0, end;
   sx = ee -(HH+(abs(emin)+lambda)*eye(n))\gg;
   x=dd.*sx;
   lambda=lambda/1.5;
end;
xh=max(0,x-0.01);
% Compute L2L1 solution
x=x0;
lambda=10;
iter=0;
%  Repeatly solving the affine-scalling trust-region step 
while (iter < maxiter),
   iter=iter+1;
   % generate affine-scaled gradient and Hessian and 
   % apply the trust-region method
   dd = abs(x);
   DD = diag(dd);
   HH = DD*ATA*DD;
   gg = dd.*(A'*(A*x-b))+mu*x;
   sx = ee -(HH+lambda*eye(n))\gg;
   x=dd.*sx;
   lambda=lambda/1.5;
end;
x1=max(0,x-0.01);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the function