%======================================================
%  Matlab implementation of the Federated-Learning vs the
%  dimension-reduced 2nd-order method based on approximate
%  Hessian using Parallel-Tangent or Conjugate directions
%  for linear regressian at three data centers
%  
%  minimize (\|A_1x-b_1\|^2+\|A_2x-b_2\|^2+\|A_3x-b_3\|^2)/2
%
%  Input Data
%      x0: any initial solution
%  Output
%      xg: latest iterative solution of the gradient method
%       x: latest iterative solution of the DRHSOM
%
%   Details can be found in Lecture Note #10, Slides 3-4
%======================================================% 
if exist('beta') ~= 1 
   beta=max(eigs(A1'*A1+A2'*A2+A3'*A3)); 
end
if exist('maxiter') ~= 1 
   maxiter=100; 
end
x=x0;
% Start the Federated-Leaning
for k=1:maxiter,
%Compute gradient at each data center
g1=A1'*(A1*x-b1);
g2=A2'*(A2*x-b2);
g3=A3'*(A3*x-b3);
g=g1+g2+g3;
x=x-g/beta;
end;
xg=x;
% Start the DRH 2nd method
x=x0;
g1=A1'*(A1*x-b1);
g2=A2'*(A2*x-b2);
g3=A3'*(A3*x-b3);
g0=g1+g2+g3;
norm(g0)
d=-g0/beta;
x=x+d;
% using a half number of iteration limit
for k=1:maxiter/2,
  g1=A1'*(A1*x-b1);
  g2=A2'*(A2*x-b2);
  g3=A3'*(A3*x-b3);
  g=g1+g2+g3;
  %Construct the linear term of reduced objectve
  gnorm=norm(g);
  gTd=g'*d;
  %Construct the dimension-reduced approximate Hessian
  % Approximate Hg
  g1=A1'*(A1*(x+g)-b1);
  g2=A2'*(A2*(x+g)-b2);
  g3=A3'*(A3*(x+g)-b3);
  Hg=(g1+g2+g3)-g;
  % Approximate Hd
  Hd=g-g0;
  % Take the Newton step
  [g'*(Hg) -Hg'*d;-Hg'*d d'*Hd]\[gnorm^2;-gTd];
  alphag=ans(1);alpham=ans(2);
  d=-alphag*g+alpham*d;
  x=x+d;
  g0=g;
end;
[norm(A1'*(A1*xg-b1)+A2'*(A2*xg-b2)+A3'*(A3*xg-b3))... 
 norm(A1'*(A1*x-b1)+A2'*(A2*x-b2)+A3'*(A3*x-b3))]
% 