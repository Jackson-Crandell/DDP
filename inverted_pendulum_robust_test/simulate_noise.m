% Propogating nominal dynamics
function [x,u] = simulate_noise(xo,u_new,Horizon,dt,sigma,u_star,x_star,L_k)

global g;
global m; 
global l; 
global I; 
global b; 

x = xo;
u_new = zeros(1,Horizon-1);
for k = 1:(Horizon-1)
      Fx(1,1) = x(2,k); %Theta1_dot
      Fx(2,1) = ((-b/I)*x(2,k))-((m*g*l)/I)*sin(x(1,k));

      G_x(2,1) = 1/I;

    u_new(k) = u_star(:,k) + L_k(:,:,k)*(x(:,k) - x_star(:,k));
    
    x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt + G_x * u_new(:,k) * sqrt(dt) * sigma * randn; %Creates trajectory based on dynamics
    u = u_new;
end