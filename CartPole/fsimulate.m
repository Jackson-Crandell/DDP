function [x] = fnsimulate(xo,u_new,Horizon,dt,dynamics,sigma)

global g; 
global m_c;
global m_p; 
global l;


x = xo;

for k = 1:(Horizon-1)
       
%       Fx(1,1) = x(2,k); 
%       Fx(2,1) = (1/(m_c + m_p*sin(x(3,k)).^2))*(m_p*sin(x(3,k))*(l*x(4,k).^2 + g*cos(x(3,k))));
%       Fx(3,1) = x(4,k);
%       Fx(4,1) = (1/(m_c + m_p*sin(x(3,k)).^2))*(-m_p*l*(x(4,k).^2)*cos(x(3,k))*sin(x(3,k))-(m_c +m_p)*g*sin(x(3,k)));
%          
%       G_x(2,1) = 1/(m_c + m_p*sin(x(3,k)).^2);
%       G_x(4,1) = -cos(x(3,k))/(l*(m_c + m_p *sin(x(3,k)).^2));
%       
% 
% 
% x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt  + G_x * u_new(:,k); %* sqrt(dt) * sigma * randn ;
x(:,k+1) = x(:,k) + dynamics.F(u_new(1, k), x(2,k), x(3,k), x(4,k)) * dt;
end