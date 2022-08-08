% Builds and solves a simple linear program

n = 3;
c = [2;1;1];

fprintf('solve the small LP \n');
cvx_begin
   variable x(n)
   dual variables y z
   minimize( c' * x )
   subject to
      y : sum(x) == 1;
      z : x >= 0;
cvx_end
disp('The optimal soluiton x of the LP is ');
x


fprintf('solve the small SOCP \n');
cvx_begin
   variable x(n)
   dual variables y z
   minimize( c' * x )
   subject to
      y : sum(x) == 1;
      z : norm(x(2:3)) <= x(1);
cvx_end
disp('The optimal soluiton x of the SOCP is ');
x

cvx_begin sdp
  % X is a PSD symmetric matrix 
  variable X(2,2) symmetric;
  X >= 0;
  variable x(3);
  % constrained matrix entries.
  X(1,1) == x(1);
  X(1,2) == x(2);
  X(2,2) == x(3);
  sum(x) == 1;
  minimize( c' * x )
cvx_end
disp('The optimal soluiton X of the SDP is ');
X


