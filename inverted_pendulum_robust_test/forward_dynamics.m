function [dx] = forward_dynamics(x,u, dt)


global g;
global m; 
global l; 
global I; 
global b; 

Fx(1,1) = x(2,1); 
Fx(2,1) = (-b/I)*x(2,1)-((m*g*l)/I)*sin(x(1,1));

G_x(2,1) = (1/I)*u;

dx = Fx * dt + G_x *dt;