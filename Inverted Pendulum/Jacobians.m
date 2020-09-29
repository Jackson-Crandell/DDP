function [A,B] = Jacobians(x,u)

global g;
global m; 
global l; 
global I; 
global b;

x1 = x(1,1);
x2 = x(2,1);

u1 = u(1,1);

A = zeros(2,2);

A(1,1) = 0; 
A(1,2) = 1;
A(2,1) = -((m*g*l)/I)*cos(x1);
A(2,2) = -b/I;

B = zeros(2,1);
B(1,1) = 0;
B(2,1) = 1/I; 

