function [A,B,C,c] = fnState_And_Control_Transition_Matrices(x,u)

global g; 
global m_c;
global m_p; 
global l; 


% Linearize along nominal trajectory x
x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

% Linearize along nominal control  u
u1 = u(1,1);


%% Compute the Jacobian of Fx 
A = zeros(4,4);

A(1,1) = 0; 
A(1,2) = 1; 
A(1,3) = 0; 
A(1,4) = 0;

A(2,1) = 0; 
A(2,2) = 0; 

a23_1 = ((m_p.*l.*(x4.^2).*cos(x3)).*(m_c + m_p.*(sin(x3).^2)) - (mp.*sin(x3).*l.*(x4.^2))*(2.*m_p.*sin(x3).*cos(x3)))./((m_c + m_p.*(sin(x3).^2)).^2); 
a23_2 = ((-2.*g.*cos(x3).*sin(x3)).*(m_c + m_p.*(sin(x3).^2)) - (g.*(cos(x3).^2)).*(2*m_p.*sin(x3).*cos(x3)))./((m_c + m_p.*(sin(x3).^2)).^2); 
a23_3 = (-u1.*2.*m_p.*sin(x3).*cos(x3))./((m_c + m_p.*(sin(x3).^2)).^2); 
A(2,3) = a23_1 + a23_2 + a23_3; 

A(2,4) = (2.*m_p.*sin(x3).*l.*x4)./((m_c + m_p.*(sin(x3).^2)).^2); 

A(3,1) = 0; 
A(3,2) = 0; 
A(3,3) = 0; 
A(3,4) = 1;

A(4,1) = 0; 
A(4,2) = 0; 

a43_1 = ((-m_p.*l.*(x4.^2).*(-sin(x3).*cos(x3) + (cos(x3).^2))).*(l.*(m_c + m_p.*(sin(x3).^2))) - (-m_p.*l.*(x4.^2).*cos(x3).*sin(x3)).*(2.*l.*m_p.*sin(x3).*cos(x3)))./((l.*(m_c + m_p.*(sin(x3).^2))).^2);
a43_2 = ((-(m_c + m_p).*g.*cos(x3)).*(l.*(m_c + m_p.*(sin(x3).^2))) - (-(m_p + m_c).*g*sin(x3)).*(2.*l.*m_p.*sin(x3).*cos(x3)))./((l.*(m_c + m_p.*(sin(x3).^2))).^2);                             
a43_3 = ((u1.*sin(x3)).*(l.*(m_c + m_p.*(sin(x3).^2))) - (-u1.*cos(x3)).*(2.*l.*m_p.*sin(x3).*cos(x3)))./((l.*(m_c + m_p.*(sin(x3).^2))).^2);
A(4,3) = a43_1 + a43_2 + a43_3; 

A(4,4) = (-2.*m_p.*l.*cos(x3).sin(x3).*x4)./((l.*(m_c + m_p.*(sin(x3).^2))).^2);


%% Compute the Jacobian of Fu 

B = zeros(4,4);

B(1,1) = 0;
B(2,1) = (1)./(m_c + m_p.*(sin(x3).^2))
B(3,1) = 0;
B(4,1) = (-cos(x3))./(l.*(m_c + m_p.*(sin(x3).^2))); 
