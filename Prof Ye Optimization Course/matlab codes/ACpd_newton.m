%======================================================
%   One-Iteration Demonstration of the Primal-Dual 
%   Algorithm for computing the analytic center of 
%      {x: A x = b, x >= 0}
%
%  Input 
%      A: m x n matrix, b: rhs vector, c: obj vector
%      x0: initial and interior feasible solution for primal
%      y0: initial and interior feasible solution for dual
%      iter (= 0): iteration count
%   
%   Details can be found in Sect. 5.4 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================%
if (iter == 0)
    [m,n]  = size(A);
    x=x0;
    y=y0;
    s=-A'*y;
end;
rk     = ones(n,1)-x.*s;
eta    =norm(rk)
Ak     = A*diag(x./s);
Mk     = Ak*A';
dy     = -Mk\(A*(rk./s));
ds     = -A'*dy;
dx     = rk./s - (x.*ds)./s;
y      = y+dy;
s      = s+ds;
x      = x+dx;
%
iter   = iter+1
%
