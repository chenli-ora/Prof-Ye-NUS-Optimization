%======================================================
%  Matlab Demonstration of the LP formulation of the
%  4-City WBC example
%
% Lecture Note #1 
%
% Construct the transportation matrix for a single demand
% distribution 
A=sparse([ones(1,4) zeros(1,4) zeros(1,4) zeros(1,4);
   zeros(1,4) ones(1,4) zeros(1,4) zeros(1,4);
   zeros(1,4) zeros(1,4) ones(1,4) zeros(1,4);
   zeros(1,4) zeros(1,4) zeros(1,4) ones(1,4);
   eye(4) eye(4) eye(4) eye(4)]);
% Remove the last row since it is redundant
A=A(1:7,:);
% Input the objective coefficientvector for a single demand
c=[0 1 1 2 1 0 2 1 1 2 0 1 2 1 1 0]';
% Three given possible demand distributions
dl=[3;3;3;0];
dm=[0;3;3;3];
dr=[3;0;3;3];
% Construct the constraint matrix for three dist. scenario
AA=sparse([A zeros(7,16) zeros(7,16) [-eye(4);zeros(3,4)];
    zeros(7,16) A zeros(7,16) [-eye(4);zeros(3,4)];
    zeros(7,16) zeros(7,16) A [-eye(4);zeros(3,4)];
    zeros(1,48) ones(1,4)]);
% Construct the objective vector for three dist. scenario
cc=[c;c;c;zeros(4,1)];
% Construct the RHS vector
bb=[zeros(4,1);dl(1:3);zeros(4,1);dm(1:3);zeros(4,1);dr(1:3);9];
% Call LP solver
[x,y] = HSDLPsolver(AA,bb,cc);
% Final WBC solution
s=x(49:52)
% 