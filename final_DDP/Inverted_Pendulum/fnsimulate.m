% Propogating nominal dynamics
function [x] = fnsimulate(xo,u_new,Horizon,dt,sigma)

global g;
global m; 
global l; 
global I; 
global b; 

x = xo;

for k = 1:(Horizon-1)
      Fx(1,1) = x(2,k); 
      Fx(2,1) = ((-b/I)*x(2,k))-((m*g*l)/I)*sin(x(1,k));

      G_x(2,1) = 1/I;

x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt + G_x * u_new(:,k) * sqrt(dt) * sigma * randn; %Creates trajectory based on dynamics
end