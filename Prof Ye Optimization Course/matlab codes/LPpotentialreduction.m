function [x,y]=LPpotentialreduction(A,b,c,x0,y0)
%
[m,n]  = size(A);
x=x0;
y=y0;
s=c-A'*y;
alpha = 0.99;
rho = 2*n;
potential=(n+rho)*log(x'*s)-ones(n,1)'*log(x.*s)
mu     = x'*s/n;
%
for k=1:10,
      gamma  = n/(n+rho);
      rk     = gamma*mu*ones(n,1)-x.*s;
      Ak     = A*diag(x./s);
      Mk     = Ak*A';
      dy     = -Mk\(A*(rk./s));
      ds     = -A'*dy;
      dx     = rk./s - (x.*ds)./s;
% Conservative step size
%theta  = min(x.*s)/((rk./x)'*(rk./s));
%theta  = alpha*sqrt(theta);
%
% Aggressive step size
      theta = min([dx./x;ds./s]);
      theta = abs(alpha/theta);
%
      y = y+theta*dy;
      s = s+theta*ds;
      x = x+theta*dx;
      mu = x'*s/n;
      potential=(n+rho)*log(x'*s)-ones(n,1)'*log(x.*s)
end;
%
%  The primal-dual potential reduction algorithm for solving LP:
%
%      maximize    c'*x
%      subject to  A x = b, x >= 0
%
%  Input 
%      A: m x n matrix, b: rhs vector, c: obj vector
%      x0: initial interior feasible point for the primal
%      y0: initial interior feasible point for the dual
