function [A,B,C,c] = fnState_And_Control_Transition_Matrices(x,u,du,dt)

global m_c;
global m_p;
global g;
global l;

x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

u1 = u(1,1);
u2 = u(2,1);

A = zeros(4,4);
A(1,2) = 1;
A(2,3) = -(m_p*g)/m_c;
A(3,4) = 1;
A(4,3) = -((m_c + m_p)*g)/(l*m_c);

B = zeros(4,1);
B(2,1) = 1/m_c;
B(4,1) = 1/(l*m_c);

C = zeros(2,4);
C(1,1) = 1;
C(2,3) = 1;




