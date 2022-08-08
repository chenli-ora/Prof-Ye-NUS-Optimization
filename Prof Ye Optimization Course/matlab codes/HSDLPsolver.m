%======================================================
%   Complete Implementation of the Homogenerous and 
%   Self-Dual Primal-Dual Potential-Reduction Algorithm 
%   for LP 
%      minimize    c'*x
%      subject to  A*x=b, 
%                  x_f free, x_p >= 0, u >= x_b (>= 0).
%
%  Input 
%      A: m x n matrix
%      b: the right-hand-side vector
%      c: the objective coefficient vector
%  Optional Input
%      u: upp bound vector for a set of nonnegative variable x_b, whose
%         indices are represented by bindx.
%      bindx: the index set for those upp-bounded nonnegative variables
%             Default value: []
%      findx: the index set for free variables x_f
%             Default value: []
%      toler: stopping tolerance, the objective value close to the      
%             optimal one in the range of the tolerance. 
%             Default value: 1.0e-6
%      beta : step size, 0 < beta < 1. 
%             Default value: .995
%
%  Output
%     x  : optimal solution
%     y  : optimal dual solution (shadow price) for equality constraints
%     w  : optimal dual solution for uppbound constraints
%     s  : optimal dual slack solution
%     zh : (infeasible) primal-dual gap history vs iteration. 
%
%  How to use it? Just type "hsdLPsolver" and hit the "return"
%   Details can be found in Sect. 5.7 of
%   L&Y, Linear and nonlinear programming, 5th edition
%======================================================% 
function [x,y,s] = HSDLPsolver(A,b,c);
%
%  See user guide at the end of this file.
%
 if exist('toler') ~= 1
   toler = 1.e-8;
 end;
 if exist('beta') ~= 1
   beta = .995;
 end;
 if exist('bindx') ~= 1
   bindx = [];
 end;
 if exist('findx') ~= 1
   findx = [];
 end;
 [nb,x]  = size(bindx);
 [nf,x]  = size(findx);
 norbc = 1+max([c;b])-min([c;b]);
 c = c/norbc;
 b = b/norbc;
 if nb > .5
     u=u/norbc;
 end
 if nf > .5,
   A = [A  -A(:,findx)];
   c = [c ; - c(findx)];
 end;
 [m,n] = size(A);
 [i,j,v] = find(A);
 ee = ones(n,1);
%
% Initialization
%
% b    = b/(norm(b)+1);
% c    = c/(norm(c)+1);
 x    = ((norm(b)+1)/sqrt(n))*ones(n,1);
 s    = ((norm(c)+1)/sqrt(n))*ones(n,1);
 if nb > 0.5,
   z    = ((norm(b)+1)/sqrt(n))*ones(nb,1);
   w    = ((norm(c)+1)/sqrt(n))*ones(nb,1);
   eez  = ones(nb,1);
 else
   z   = 1;
   w   = 0;
   eez = 0;
   rb  = 0;
   u   = 0;
 end;
 y    = zeros(m,1);
 tau0  = 1;
 kappa0= (norm(c)+1)*(norm(b)+1)/n;
% tau0  = (norm(b)+1)/sqrt(n);
% kappa0= (norm(c)+1)/sqrt(n);
 tau   = tau0;
 kappa = kappa0;
%
% Compute initial residuals
%
 mu0  = (x'*s+z'*w+tau*kappa)/(n+nb+1);;
 mu   = mu0;
 rp  = tau*b - A*x;
%
 if nb > 0.5,  
   rb  = x(bindx) + z - tau*u; 
 end;
%
 rd  = tau*c - A'*y -s;
 rd(bindx) = rd(bindx) + w;
 obp = c'*x;
 obd = b'*y-u'*w;
 rg  = obp - obd + kappa;
 zh   = (obp-obd)/tau;
%
% Iteration
%
 gamma=1/sqrt(n+nb+1);
 go   = 1;
 iter=1;
 while go >= toler,
%
%  Call hsdLP sub-solver
%
      cvx = ((gamma*mu)*ee)./x - s;
      cvz = ((gamma*mu)*eez)./z - w;
      woz = w./z;
      r1  = cvx -rd;
      r1(bindx) = r1(bindx) - cvz - rb.*woz;
      r2  = c;
      r22 = c;
      r2(bindx)  = r2(bindx)  - u.*woz;
      r22(bindx) = r22(bindx) + u.*woz;
      d   = s./x;
      d(bindx) = d(bindx)+woz;
      d   = sqrt(d);
      AD = sparse(i,j,v./d(j), m, n);
%
% Find a scaling parameter
%
%alpha = (max([abs(max(vd));abs(min(vd))]))/1000;
%
      r1  = [  r1./d  ; rp];
      r2  = [ -r2./d  ;  b];
      r22 = [-r22./d  ;  b];
%
      ss=[speye(n) AD'; AD sparse(m,m)]\[r1 r2];
      ss(n+1:n+m,:)=-ss(n+1:n+m,:);
%
% get dtau
%
      dtau=(gamma*mu/tau - kappa + rg + u'*cvz + u'*(woz.*rb) - r22'*ss(:,1));
      dtau = dtau/(kappa/tau+u'*(woz.*u)+ r22'*ss(:,2));
      ss   = ss(:,1)+dtau*ss(:,2);
%
% get dx 
%
      dx  = ss(1:n)./d;
%
      dy  = ss(n+1:n+m);
%
% get ds 
%
      ds  = cvx - (s.*dx)./x;
%
% get dz and dw if bounds exist
%
      if nb > 0.5,
      dz  = dtau*u - dx(bindx) - rb;
      dw  = cvz - woz.*dz;
      end;
%
% get dkappa
%
      dkappa= gamma*mu/tau - kappa - kappa*dtau/tau;
%
      if nb > 0.5;
      ratp  = beta/abs(min([dx./x;dz./z;dtau/tau;dkappa/kappa]));
      ratd  = beta/abs(min([ds./s;dw./w;dtau/tau;dkappa/kappa]));
      else
      ratp  = beta/abs(min([dx./x;dtau/tau;dkappa/kappa]));
      ratd  = beta/abs(min([ds./s;dtau/tau;dkappa/kappa]));
      end;
      x    = x    + ratp*dx;
      y    = y    + ratd*dy;
      s    = s    + ratd*ds;
      if nb > 0.5,
      z    = z    + ratp*dz;
      w    = w    + ratd*dw;
      end;
      taup = tau  + ratp*dtau;
      taud = tau  + ratd*dtau;
      tau  = min([taup; taud]);
      if taup <= taud,
      kappa = kappa + ratp*dkappa;
      y    = (tau/taud)*y;
      s    = (tau/taud)*s;
      if nb > 0.5, w = (tau/taud)*w; end;
      else
      kappa = kappa + ratd*dkappa;
      x     = (tau/taup)*x;
      if nb > 0.5, z = (tau/taud)*z; end;
      end;
% End of the subsolver
%
%  Compute new residuals
%   tau/kappa
   mu  = (x'*s+z'*w+tau*kappa)/(n+nb+1);
   rp  = tau*b - A*x;
   if nb > 0.5,  
     rb  = x(bindx) + z - tau*u; 
   end;
   rd  = tau*c - A'*y -s;
   rd(bindx) = rd(bindx) + w;
   obp = c'*x;
   obd = b'*y-u'*w;
   rg  = obp - obd + kappa;
%   pause
   go = max([norm([rp;rb])/(tau+norm([x;z])); norm(rd)/(tau+norm([s;w]))]);
   go = max([go; abs(obp-obd)/(tau+abs(obd))]);
%
%  Check infeasibility
% 
   if (tau*kappa0/(tau0*kappa) < toler) & (mu/mu0 < toler/n),
     if (obp < 0) & (obd > 0),
      disp('Both primal and dual are infeasible.'); go=-1;
     end;
     if obp < 0, 
      disp('Dual is infeasible.')  ; go=-1; 
     end;
     if obd > 0,
      disp('Primal is infeasible.'); go=-1;
     end;
   end;
   zh   = [zh (obp-obd)/tau];
%   
%  Adjust centering weight
%
   gamma = 1/sqrt(n+nb+1);
   if nb > 0.5, 
     if min([x.*s;z.*w;tau*kappa])/mu >= (1-beta),
       gamma = 0;
     end;
   else
     if min([x.*s;tau*kappa])/mu >= (1-beta),
       gamma = 0;
     end;
   end;
%
   iter = iter + 1;
 end;
%
% Make an output
%
% b    = b/(1-norm(b));
% c    = c/(1-norm(c));
 s(bindx) = s(bindx) - w;
 if go >= 0,
   x = norbc*x/tau;
   z = norbc*z/tau;
   s = norbc*s/tau;
   y = norbc*y/tau;
   w = norbc*w/tau;
 end
 if nf > .5,
   n = n - nf;
   x(findx) = x(findx) - x(n+1:n+nf);
   x = x(1:n);
   s = s(1:n);
   A = A(:,1:n);
   c = c(1:n);
 end;
 c = c*norbc;
 b = b*norbc;
 if nb > .5
     u=u*norbc;
 end;
%
%  Technical References
%
%  Y. Ye, M. J. Todd and S. Mizuno, "An $O(\sqrt{n}L)$-iteration
%  homogeneous and self-dual linear programming algorithm," (1992)
%  
%  X. Xu, P. Hung and Y. Ye, "A simplified homogeneous and self-dual
%  linear programming algorithm and its implementation," (1996).
%