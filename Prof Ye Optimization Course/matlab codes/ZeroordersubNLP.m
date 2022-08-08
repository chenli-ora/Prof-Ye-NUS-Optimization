%------------------------------------------------------------------------------
%  Subroutine function called by ZeroorderNLP
%  Variable Glossary:
%------------------------------------------------------------------------------
%
%ALP (4):   vector of ??. Alp(4) is "alpha"
%CH:        "change"
%CONSTRAINT:vector of constraint values
%G:         gradient
%GAP (2):   lower and upper "gap"
%LPB (2):   vector flag which indicates the presence of parameter bounds
%           and or inequality constraints.
%             LPB(1) refers to parameter bounds, it is 0 if there are 
%                    none, 1 if there are one or more.  
%             LPB(2) refers to constraints of either type.
%msg:       long error message string
%NC:        total number of constraints (=NEC+NIC)
%NEC:       number of equality constraints
%NIC:       number of inequality constraints
%NP:        number of parameters
%om:        string for optimization messages
%P0:        parameter vector on input
%P:         parameter vector on output
%SOB (3):   vector of ??
%Z (5):     vector of ??
%
function [p,y,h,l]=ZeroordersubNLP(p0,op,yy,ob,pb,h,l,par)
%
rho=op(1); maxit=op(2); delta=op(3); tol=op(4);   nec=op(5); 
nic=op(6); np=op(7);    nc=nec+nic;  npic=np+nic; lpb=op(8:9); ch=1;
clear op
alp=[0 0 0];
%
% make the scale for the cost, the equality constraints, the inequality
% constraints, and the parameters
%
if nec>0,   
  scale=[ob(1);ones(nec,1)*max(abs(ob(2:nec+1)))];
else        
  scale=1;
end
if lpb(2)<=0, 
  scale=[scale; p0];
else        
  scale=[scale; ones(size(p0))];
end
scale=min(max(abs(scale),tol),1/tol);
%
% scale the cost, the equality constraints, the inequality constraints, 
% the parameters (inequality parameters AND actual parameters), 
% and the parameter bounds if there are any
% Also make sure the parameters are no larger than (1-tol) times their bounds
%
ob=ob./scale(1:nc+1);  
p0=p0./scale(nec+2:nc+np+1);
if lpb(2)>=0.5,
  if lpb(1)<=0.5,
    mm=nic;
  else 
    mm=npic;
  end
  pb=pb./[scale(nec+2:nec+mm+1) scale(nec+2:nec+mm+1)];
end
%
% scale the lagrange multipliers and the Hessian
%
if nc>0.5,
  yy=scale(2:nc+1).*yy/scale(1);
end
h=h.*(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)')/scale(1);
%
om=['SOLNP--> ';'         '];
msg=[om(1,:) 'Redundant constraints were found. Poor              '
     om(2,:) 'intermediate results may result.  Suggest that you  '
     om(2,:) 'remove redundant constraints and re-OPTIMIZE.       '];
%
j=ob(1);
a=[0*ones(nec,nic);-eye(nic)];
g=0*ones(npic,1);
%
if nc>0.5,
  constraint=ob(2:nc+1);
  for i=1:np,
    p0(nic+i)=p0(nic+i)+delta;
    ob=cost(p0(nic+1:npic).*scale(nc+2:nc+np+1),0,par)./scale(1:nc+1);
    g(nic+i)=(ob(1)-j)/delta;
    a(:,nic+i)=(ob(2:nc+1)-constraint)/delta;
    p0(nic+i)=p0(nic+i)-delta;
  end
  if nic>0.5,
    constraint(nec+1:nec+nic)=constraint(nec+1:nec+nic)-p0(1:nic);
  end
  if cond(a)>1/eps,
    disp(msg);
  end
  b=a*p0-constraint;
end
%
msg=[om(1,:) 'The linearized problem has no feasible     '
     om(2,:) 'solution.  The problem may not be feasible.'];
%
if nc>0.5,
  ch=-1;
  alp(1)=tol-max(abs(constraint));
  if alp(1)<=0,
    ch=1;
    if lpb(2)==0,
      p0=p0-a'*((a*a')\constraint);
      alp(1)=1;
    end
  end
  if alp(1)<=0,
    p0(npic+1)=1;
    a=[a, -constraint];
    c=[0*ones(1,npic), 1];
    dx=ones(npic+1,1);
    go=1; 
    minit=0;
    while go>=tol,
      minit=minit+1;
      gap=[p0(1:mm,1)-pb(:,1),pb(:,2)-p0(1:mm,1)];
      gap=sort(gap')';
      dx(1:mm)=gap(:,1);
      dx(npic+1)=p0(npic+1);
      if lpb(1)==0,
        dx(mm+1:npic)=max([dx(1:mm);100])*ones(npic-mm,1);
      end
      y=(a*diag(dx))'\(dx.*c');
      v=dx.*(dx.*(c'-a'*y));
      if v(npic+1)>0,
        z=p0(npic+1)/v(npic+1);
        for i=1:mm,
          if v(i)<0,
            z=min(z,-(pb(i,2)-p0(i))/v(i));
          elseif v(i)>0, 
            z=min(z,(p0(i)-pb(i,1))/v(i)); 
          end
        end
        if z>= p0(npic+1)/v(npic+1),
          p0=p0-z*v; 
        else
          p0=p0-0.9*z*v; 
        end
        go=p0(npic+1);
        if minit >= 10, 
          go=0; 
        end
      else
        go=0;
        minit=10;
      end
    end
    if minit>=10,
      disp(msg),
    end
    a=a(:,1:npic); 
    b=a*p0(1:npic);
  end
end
%
clear constraint c z v gap;
%
p=p0(1:npic); 
y=0; 
if ch>0,
  ob=cost(p(nic+1:npic).*scale(nc+2:nc+np+1),0,par)./scale(1:nc+1);
end
j=ob(1);
%
if nic>0.5,
  ob(nec+2:nc+1)=ob(nec+2:nc+1)-p(1:nic);
end
if nc>0.5,
  ob(2:nc+1)=ob(2:nc+1)-a*p+b;
  j=ob(1)-yy'*ob(2:nc+1)+rho*norm(ob(2:nc+1))^2;
end
%
minit=0;
while minit<maxit,
  minit=minit+1
  if ch>0,
    for i=1:np,
      p(nic+i)=p(nic+i)+delta;
      obm=cost(p(nic+1:npic).*scale(nc+2:nc+np+1),0,par)./scale(1:nc+1);
      if nic>0,
        obm(nec+2:nc+1)=obm(nec+2:nc+1)-p(1:nic);
      end
      if nc>0,
        obm(2:nc+1)=obm(2:nc+1)-a*p+b;
        obm=obm(1)-yy'*obm(2:nc+1)+rho*norm(obm(2:nc+1))^2;
      end
      g(nic+i)=(obm-j)/delta;
      p(nic+i)=p(nic+i)-delta;
    end
    if nic>0.5,
      g(1:nic)=0*yy(nec+1:nc);
    end
  end
  if minit>1,
    yg=g-yg;
    sx=p-sx;
    sc(1)=sx'*h*sx;
    sc(2)=sx'*yg;
    if sc(1)*sc(2)>0,
      sx=h*sx;
      h=h-(sx*sx')/sc(1)+(yg*yg')/sc(2);
    end
  end
  dx=0.01*ones(npic,1);
  if lpb(2)>0.5,
    gap=[p(1:mm,1)-pb(:,1),pb(:,2)-p(1:mm,1)];
    gap=sort(gap')';
    gap=gap(:,1)+sqrt(eps)*ones(mm,1);
    dx(1:mm,1)=ones(mm,1)./gap;
    if lpb(1)<=0,
      dx(mm+1:npic)=min([dx(1:mm);0.01])*ones(npic-mm,1);
    end
  end
  go=-1;
  l=l/10;
  while go<=0,
    c=chol(h+l*diag(dx.*dx));
    c=inv(c);
    yg=c'*g;
    if nc<=0,
      u=-c*yg;
    else 
      y=(c'*a')\yg;
      u=-c*(yg-(c'*a')*y);
    end
    p0=u(1:npic)+p;
    if lpb(2)<=0,
      go=1;
    else
      go=min([p0(1:mm)-pb(:,1);pb(:,2)-p0(1:mm)]);
      l=3*l;
    end
  end
  alp(1)=0;ob1=ob;ob2=ob1;sob(1)=j;sob(2)=j;
  pt(:,1:2)=[p p];alp(3)=1.0;pt(:,3)=p0;
  ob3=cost(pt(nic+1:npic,3).*scale(nc+2:nc+np+1),-minit,par)./scale(1:nc+1);
  sob(3)=ob3(1);
  if nic>0.5,
    ob3(nec+2:nc+1)=ob3(nec+2:nc+1)-pt(1:nic,3);
  end
  if nc>0.5,
    ob3(2:nc+1)=ob3(2:nc+1)-a*pt(:,3)+b;
    sob(3)=ob3(1)-yy'*ob3(2:nc+1)+rho*norm(ob3(2:nc+1))^2;
  end
  go=1;
  while go>tol,
    alp(2)=(alp(1)+alp(3))/2;
    pt(:,2)=(1-alp(2))*p+alp(2)*p0;
    ob2=cost(pt(nic+1:npic,2).*scale(nc+2:nc+np+1),0,par)./scale(1:nc+1);
    sob(2)=ob2(1);
    if nic>0.5,
      ob2(nec+2:nc+1)=ob2(nec+2:nc+1)-pt(1:nic,2);
    end
    if nc>0.5,
      ob2(2:nc+1)=ob2(2:nc+1)-a*pt(:,2)+b;
      sob(2)=ob2(1)-yy'*ob2(2:nc+1)+rho*norm(ob2(2:nc+1))^2;
    end
    obm=max(sob);
    if obm<j,
      obn=min(sob);
      go=tol*(obm-obn)/(j-obm);
    end
    if sob(2)>=sob(1),
      sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
    elseif sob(1)<=sob(3),
      sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
    else
      sob(1)=sob(2);ob1=ob2;alp(1)=alp(2);pt(:,1)=pt(:,2);
    end
    if go>=tol,
      go=alp(3)-alp(1);
    end
  end
  sx=p;yg=g;ch=1;
  obn=min(sob);
  if j<=obn,
    maxit=minit;
  end
  reduce=(j-obn)/(1+abs(j));
  if reduce<tol,
    maxit=minit;
  end
  if sob(1)<sob(2),
    j=sob(1);p=pt(:,1);ob=ob1;
  elseif sob(3)<sob(2),
    j=sob(3);p=pt(:,3);ob=ob3;
  else 
    j=sob(2);p=pt(:,2);ob=ob2;
  end
  clear ob1 ob2 ob3 pt;
end;
p=p.*scale(nec+2:nc+np+1);  % unscale the parameter vector
if nc>0.5,
  y=scale(1)*y./scale(2:nc+1); % unscale the lagrange multipliers
end
h=scale(1)*h./(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)');
%
if reduce>tol,
  disp([...
  om(1,:) 'Minor optimization routine did not converge in the ';...
  om(2,:) 'specified number of minor iterations.  You may need';...
  om(2,:) 'to increase the number of minor iterations.        '])
end
return
%