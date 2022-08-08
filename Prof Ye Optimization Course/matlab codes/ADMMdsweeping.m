%======================================================
%  Matlab demonstration of the diverging example for
%  of the ADMM with three blocks:
%      min   0*x1 + 0*x2 + 0*x3 
%      s.t.    x1 + x2 + x3 =0
%              x1 + x2+2*x3 =0
%              x1+2*x2+2*x3 =0
%
%  Algorithm: Double-Sweep of multi-block ADMM 
%      
%  Input: any initial solution x0 and multipliers y0
%         maxiter: the number of iterations
%
%  Output  x: final slution, y: final multipliers
%
%  Details can be found in Example 2 of Sect. 14.6 
%  L&Y, Linear and nonlinear programming, 5th edition
%======================================================%  
function [x, y] = ADMMdsweeping(x0,y0,maxiter)
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
  b=-x(2)*a2-x(3)*a3;
  x(1)=a1\(b+y);
  b=-x(1)*a1-x(3)*a3;
  x(2)=a2\(b+y);
  b=-x(1)*a1-x(2)*a2;
  x(3)=a3\(b+y);
  b=-x(1)*a1-x(3)*a3;
  %sweeping back
  x(2)=a2\(b+y);
  b=-x(2)*a2-x(3)*a3;
  x(1)=a1\(b+y);
  %updating the dual
  y=y-(A*x);
  history = [history norm(A*x)];
end;
plot(1:k, history(1:k), '-','linewidth',2);
xlabel('Iteration: Double-Sweeping Method');
ylabel('Absolute Error');
% 