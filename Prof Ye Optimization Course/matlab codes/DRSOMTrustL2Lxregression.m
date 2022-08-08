%======================================================
%  Matlab demonstration of the regularized
%  linear regression with 0.5-norm and 1-norm (LASSO) 
%  regularizations (compressed sensing)
%
%      min   0.5*\|Ax-b\|^2 + mu* sum |x(j)|^{0.5}
%
%  Algorithm: DRSOM affine-scaling trust region
%  Input 
%      A: Sparse constraint matrix.
%      b: the right-hand vector
%      x0:initial solution
%      mu0: initial regularization weight
%  Output
%     xh2 : solution of 1/2-norm by 2-dimension SOM
%  
%   Problem can be found homework #7.17 in Sect. 7.2 
%   and Algorithm in Sect. 8.7 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
function [xh2]=DRSOMTrustL2Lxregression(A,b,x0,mu,maxiter)
if exist('mu') ~= 1 
   mu=0.5; 
end
if exist('maxiter') ~= 1 
   maxiter=30; 
end
[m,n]=size(A);
ATA=A'*A;
ee=ones(n,1);
if exist('x0') ~= 1 
   x0=ones(n,1); 
end
%
% Compute L2Lhalf solution with dimension-reduced SOM
%
x=x0;
lambda=20;
iter=0;
g =(A'*(A*x-b))+(mu/2)*((1./sqrt(abs(x))).*sign(x));
d=-g/lambda;
x=x+d;
%  Repeatly solving the trust-region step 
while (iter < maxiter),
   iter=iter+1;
   % apply the ellipsoidal trust-region method
   dd = abs(x);
   % compute the gradient
   g =(A'*(A*x-b))+(mu/2)*((1./sqrt(dd)).*sign(x));
   % compute Hessian
   HK = ATA-(mu/4)*diag(1./(dd.^(1.5))); %compute Hessian
   % construct the reduced 2-dimensonal Hessian
   Hd=HK*d;
   QK = [g'*(HK*g) -g'*Hd;-g'*Hd d'*Hd];
   % construct the ellisoidal matrix
   Dg = g./dd; %scale gradient
   Dd = d./dd; %scale the momentum direction
   DK = [Dg'*Dg -Dg'*Dd; -Dg'*Dd Dd'*Dd];
   %
   emin=min(eig(QK,DK));
   if (emin > 0), emin=0; end;
   %(QK+(abs(emin)+lambda)*DK)\[norm(g)^2;-g'*d];
   (QK+(abs(emin)+lambda)*DK)\[norm(g)^2;-g'*d];
   alphag=ans(1);alpham=ans(2);
   d=-alphag*g+alpham*d;
   x=x+d;
   lambda=max(1.e-8,lambda/1.5);
end;
xh2=max(0,x-0.05);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the function