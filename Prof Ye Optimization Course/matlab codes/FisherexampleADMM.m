%======================================================
%  Matlab demonstration of solving Fisher social 
%  optimizarion by two-block ADMM
%      maximize  5ln(u1)+8ln(u2)
%      s.t.      x1 + x2 =1
%                x3 + x4 =1
%                2x1+x3 -u1=0
%                3x1+x4 -u2=0
%                all variables >=0
%
%  
%  Algorithm: Two-Block ADMM by coping variables
%             xi=si and requring si>=0. 
%  Output  x: good allocation, y: prices of goods
%
%  Problem can be found from Example 3 in Sect. 11.6
%  and ADMM method in Sect.14.5 of
%  L&Y, Linear and nonlinear programming, 5th edition
%======================================================%  
% Set constraint matrix for the first block variable x(1:4)
A1=[1 1 0 0;0 0 1 1;2 0 1 0;0 3 0 1;eye(4)];
% Set constraint matrix for the second block variable x(5:10), including
% two u and four s variables in the lecture note
A2=[zeros(2,6);-eye(6)];
% Precompute the inverse (A1'*A1) that would be repeatdly used later.
A1INV=inv(A1'*A1);
% Set the b vector in dimension 8
b=[1;1;zeros(6,1)];
% set the budgets of the two agent
w=[5;8];
beta=1;
% Set the initial solution where the two goods are evenly distributed
x=[1/2;1/2;1/2;1/2;3/2;2;1/2;1/2;1/2;1/2];
% Set the multipliers to zeros
y=zeros(8,1);
for k=1:100;
 %Update the first block variable:
 %   prepare the fixed coefficients
 b1=b-A2*x(5:10);
 %   minimize a quadratic objective
 x(1:4)=A1INV*(A1'*(b1+y/beta));
 %Update the second block variable: First two u and later 
 %four s variables
 %  prepare the fixed coefficientts 
 b2=A1*x(1:4)-b;
 tempb=beta*b2(3)-y(3);
 %  find the right root for u1
 x(5)=(tempb+sqrt(tempb^2+4*beta*w(1)))/(2*beta);
 %  find the right root for u2
 tempb=beta*b2(4)-y(4);
 x(6)=(tempb+sqrt(tempb^2+4*beta*w(2)))/(2*beta);
 %  Update the four s variables that need to be nonnegative
 x(7:10)=max(0,b2(5:8)-y(5:8)/beta);
 % Update the multipliers
 y=y-beta*(A1*x(1:4)+A2*x(5:10)-b);
end;
x=x(1:4),
y=-y(1:2),
%