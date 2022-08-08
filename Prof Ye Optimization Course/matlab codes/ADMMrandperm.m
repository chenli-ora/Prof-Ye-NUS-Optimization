%======================================================
%  Matlab demonstration of the diverging example for
%  of the ADMM with three blocks:
%      min   0*x1 + 0*x2 + 0*x3 
%      s.t.    x1 + x2 + x3 =0
%              x1 + x2+2*x3 =0
%              x1+2*x2+2*x3 =0
%
%  Algorithm: Random-Permutation of multi-block ADMM 
%      
%  Input: any initial solution x0 and multipliers y0
%
%  Output  x: final slution, y: final multipliers
%         maxiter: the number of iterations
%
%  Details can be found in Example 3 of Sect. 14.6 
%  L&Y, Linear and nonlinear programming, 5th edition
%======================================================%  
function [x, y] = ADMMrandperm(x0,y0,maxiter)
if exist('maxiter') ~= 1 
   maxiter=500; 
end
%
x=x0;
y=y0;
a1=[1;1;1];
a2=[1;1;2];
a3=[1;2;2];
A=[a1 a2 a3];
history = [];
%
for k=1:maxiter,
  pp=randperm(3);
  b=-x(pp(2))*A(:,pp(2))-x(pp(3))*A(:,pp(3));
  x(pp(1))=A(:,pp(1))\(b+y);
  b=-x(pp(1))*A(:,pp(1))-x(pp(3))*A(:,pp(3));
  x(pp(2))=A(:,pp(2))\(b+y);
  b=-x(pp(1))*A(:,pp(1))-x(pp(2))*A(:,pp(2));
  x(pp(3))=A(:,pp(3))\(b+y);
  y=y-(A*x);
  history = [history norm(A*x)];
end;
plot(1:k, history(1:k), '-','linewidth',2);
xlabel('Iteration: Randomly Permuted Method');
ylabel('Absolute Error');
% 