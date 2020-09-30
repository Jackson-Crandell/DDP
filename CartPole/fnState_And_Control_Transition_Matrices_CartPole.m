function [A,B] = fnState_And_Control_Transition_Matrices_CartPole(x,u)

global m_c;
global m_p;
global g;
global l;

x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

u1 = u(1,1);

A = zeros(4,4);
A(1,2) = x2; 
A(2,3) = (m_p*(cos(x3)*(l*(x4)^2+g*cos(x3))-g*sin(x3).^2)*(m_c+m_p*sin(x3).^2)-u1*m_p*sin(2*x3)*sin(x3)*(l*(x4).^2+g*cos(x3))) / (m_c+m_p*sin(x3).^2).^2;
A(3,4) = x4;
A(4,3) = (-m_p*sin(2*(x3))*(-l*m_p*(x4.^2)*sin(2*x3)-2*g*sin(x3)*(m_c+m_p)) + 2*(-l*m_p*(x4.^2)*cos(2*x3) - g*cos(x3)*(m_c + m_p))*(m_c+m_p*sin(x3).^2)) / (2*l*(m_c+m_p*sin(x3).^2));
A(2,4) = (2*m_p*sin(x3)*l*x4) / (m_c + m_p*sin(x3).^2);
A(4,4) = (-2 * m_p*x4*sin(2*x3)) / (m_c + m_p *sin(x3).^2);


B = zeros(4,1);
B(2,1) = 1/(m_c + m_p*sin(x3).^2);
B(4,1) = -cos(x3)/(l*(m_c + m_p *sin(x3).^2));




