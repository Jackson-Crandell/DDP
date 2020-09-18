global M;
global m;
global b;
global I;
global g;
global l;

M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.81;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator in the A and B matrices

 B_ = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];

A = zeros(4,4);
A(1,2) = 1;
A(2,2) = -(I+m*l^2)*b/p;
A(2,3) = (m^2*g*l^2)/p;
A(4,2) = -(m*l*b)/p ;
A(3,4) = 1;
A(4,3) = m*g*l*(M+m)/p;
B = zeros(4,1);
B(2,1) = (I+m*l^2)/p;
B(4,1) = m*l/p;
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];



