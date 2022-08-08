%======================================================
%  Matlab demonstration of the regular SDM and Affine
%  Scaling SDM for linear regression with the 0.5-norm
%  regularizations (compressed sensing)
%
%      min   0.5*\|Ax-b\|^2 + mu* sum |x(j)|^{0.5}
%
%  Algorithm: abs(x)-affine-scaling and steepest descent
%  Input 
%      A: Sparse constraint matrix.
%      b: the right-hand vector
%      x0:initial solution
%      mu: regularization weight
%  Output
%     xa  : solution from affine-scaling sdm
%     xg   : solution from rgular firt-order sdm
%
%   Problem can be found in homework #7.17 of Sect. 7.2 
%   and Algorithm in Sect. 8.5 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
function [xa,xg]=affineL2Lxregression(A,b,x0,mu,maxiter,toler)
if exist('mu') ~= 1 
   mu=0.5; 
end
if exist('maxiter') ~= 1 
   maxiter=500; 
end
[m,n]=size(A);
ATA=A'*A;
beta=eigs(ATA,1);
ee=ones(n,1);
% 
x=x0;
iter=0;
%  Repeat the affine-scalling SDM step 
while (iter < maxiter),
   iter=iter+1;
   % generate affine-scaled gradient 
   dd = abs(x);
   gg = dd.*(A'*(A*x-b))+(mu/2)*(sqrt(dd).*sign(x));
   % update x
   x=dd.*(ee -(1/beta)*gg);
end;
xa=max(0,x-0.05);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute regular SDM solution
x=x0;
iter=0;
while (iter < maxiter),
   iter=iter+1;
   % generate gradient 
   gg = (A'*(A*x-b))+(mu/2)*((1./sqrt(dd)).*sign(x));
   % update x
   x=x -(1/beta)*gg;
end;
xg=max(0,x-0.05);
% end of the function