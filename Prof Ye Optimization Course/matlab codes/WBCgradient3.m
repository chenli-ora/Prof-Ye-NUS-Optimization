% Construct the WD transportation data
A=sparse([ones(1,4) zeros(1,4) zeros(1,4) zeros(1,4);
   zeros(1,4) ones(1,4) zeros(1,4) zeros(1,4);
   zeros(1,4) zeros(1,4) ones(1,4) zeros(1,4);
   zeros(1,4) zeros(1,4) zeros(1,4) ones(1,4);
   eye(4) eye(4) eye(4) eye(4)]);
% Remove the last row since it is redundant
A=A(1:7,:);
c=[0 1 1 2 1 0 2 1 1 2 0 1 2 1 1 0]';
% Three given distributions
bl=[3;3;3;0];
bm=[0;3;3;3];
br=[3;0;3;3];
x=9*ones(4,1)/4;
for k=1:200,
  % compute the gradient of three WD functions
  b=[x;bl(1:3)];
  [xij,y,s] = HSDLPsolver(A,b,c);
  yl=y(1:4);
  b=[x;bm(1:3)];
  [xij,y,s] = HSDLPsolver(A,b,c);
  ym=y(1:4);
  b=[x;br(1:3)];
  [xij,y,s] = HSDLPsolver(A,b,c);
  yr=y(1:4);
  % Compute the average of the gradients
  g=(yl+ym+yr)/3;
  % Gradient Projection
  %[x]=Convqp(eye(4),-(x-0.3*g),-eye(4),zeros(4,1),ones(1,4),9);
  x=max(0,x-0.1*g);
  x=(9/sum(x))*x;
end;
% Steepest Descent and Projection Method for computing the Barycenter of
% the WD Course Example of three given distributions.
% Input x0: initial distribution