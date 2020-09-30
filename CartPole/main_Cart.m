%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE4803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% Variables for the cart-pole problem
global g; 
global m_c;
global m_p; 
global l; 

g = 9.81;  		% gravity 
m_c = 1.0; 		% mass of the cart
m_p = 0.01;		% mass of the pole
l = 0.25;  		% length of the pole

% Horizon 
Horizon = 800; % 1.5sec

% Number of Iterations
num_iter = 90;

% Obtain expressions for F, Fx, Fu, & Fb.
dynamics = fnDynamics(m_c, m_p, l, g);
dynamics;

% Discretization
dt = 0.01; % .01 * 300 = 3 seconds

% Weight in Final State: (part of terminal cost)
Q_f = eye(4,4);
Q_f(1,1) = 100;		%Penalize more w.r.t errors in position
Q_f(2,2) = 100;
Q_f(3,3) = 10;    	%Penalize less for errors in velocity
Q_f(4,4) = 10;

% Weight in the Control:
R = 1 * eye(1,1); 	% Weight control equally

% Initial Configuration: (Initial state)
xo = zeros(4,1);
xo(1,1) = 0;
xo(2,1) = 0;
xo(3,1) = 0;
xo(4,1) = 0;

% Initial Control:
u_k = zeros(1,Horizon-1);	% Horizon -1 b/c at last time step we have no control
du_k = zeros(1,Horizon-1);


% Initial trajectory:
x_traj = zeros(4,Horizon);
 

% Target: (Terminal States)
p_target(1,1) = 0;          % x         (position of cart)
p_target(2,1) = 0;          % x_dot     (velocity of the cart)
p_target(3,1) = pi;         % theta     (angular position of pole)
p_target(4,1) = 0;          % theta_dot (angular velocity of the pole)


% Learning Rate:c
gamma = 0.4; 
 
for k = 1:num_iter 	%Run for a certain number of iterations

%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
for  j = 1:(Horizon-1) 	%Discretize trajectory for each timestep
	 
	% linearization of dynamics (Jacobians dfx and dfu)
	[dfx,dfu] = fnState_And_Control_Transition_Matrices_CP(x_traj(:,j),u_k(:,j),du_k(:,j),dt, dynamics);

	% Quadratic expansion of the running cost around the x_trajectory (nominal) and u_k which is the nominal control
	[l0,l_x,l_xx,l_u,l_uu,l_ux] = fCost(x_traj(:,j), u_k(:,j), j,R,dt); %for each time step compute the cost
	L0(j) = dt * l0; 			% zero order term (scalar)
	Lx(:,j) = dt * l_x; 		% gradient of running cost w.r.t x (vector)
	Lxx(:,:,j) = dt * l_xx; 	% Hessian of running cost w.r.t x (matrix)
	Lu(:,j) = dt * l_u; 		% gradient of running cost w.r.t u (vector)
	Luu(:,:,j) = dt * l_uu; 	% Hessian of running cost w.r.t u (matrix)
	Lux(:,:,j) = dt * l_ux; 	% Hessian of running cost w.r.t ux (matrix)

	A(:,:,j) = eye(4,4) + dfx * dt; 	% This is PHI in notes (Identity matrix) + gradient of dynamics w.r.t x * dt
	B(:,:,j) = dfu * dt; 				%B matrix in notes is the linearized contols
end

%------------------------------------------------> Boundary Conditions
% initialize value function 
Vxx(:,:,Horizon)= Q_f; 																		% Initialize Hessian of value function (Matrix)                                                                                         
Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target);  										% Gradient of value function (Vector)
V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 	% Value function (scalar)


%------------------------------------------------> Backpropagation of the Value Function
for j = (Horizon-1):-1:1
		 
	 Q = L0(j) + V(:,j+1);
     Q_x = Lx(:,j) + A(:,:,j)'*Vx(:,j+1);
     Q_xx = Lxx(:,:,j) + A(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);
     Q_u  = Lu(:,j) + B(:,:,j)'*Vx(:,j+1);
     Q_uu = Luu(:,:,j) + B(:,:,j)'*Vxx(:,:,j+1)*B(:,:,j);
     Q_ux = Lux(:,:,j) + B(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);
     
     inv_Q_uu = inv(Q_uu);
	 L_k(:,:,j)= - inv_Q_uu*Q_ux;   % Feedback term
	 l_k (:,j) = - inv_Q_uu*Q_u;    % Feedforward term
	 
	 Vxx(:,:,j) = Q_xx - Q_ux'*inv_Q_uu*Q_ux;
	 Vx(:,j)= Q_x - Q_ux'*inv_Q_uu*Q_u;
	 V(:,j) = Q - 0.5*Q_u'*inv_Q_uu*Q_u;
     
end 


%----------------------------------------------> Find the controls
dx = zeros(4,1); % dx is initially zero because we start from the same point
for i=1:(Horizon-1)    
   du = l_k(:,i) + L_k(:,:,i) * dx; 		%Feedback Controller 
   dx = A(:,:,i) * dx + B(:,:,i) * du; 		% As we propagate forward, we use the linearized dynamics to approximate dx (this is the error from the nominal trajectory)
   u_new(:,i) = u_k(:,i) + gamma * du; 		% Update controls with gamma to prevent controls from updating too fast
end

u_k = u_new; 	%Update nominal trajectory (u_k) for new updated controls


%---------------------------------------------> Simulation of the Nonlinear System
[x_traj] = fsimulate(xo,u_new,Horizon,dt,dynamics,0);	 %Create new nominal trajectory based on new control (u_new)
[Cost(:,k)] =  fCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
x1(k,:) = x_traj(1,:);
 

fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end

time(1)=0;
for i= 2:Horizon
	time(i) =time(i-1) + dt;  
end

figure(1);
subplot(3,2,1)
hold on
plot(time,x_traj(1,:),'linewidth',4);  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta 1','fontsize',20); 
xlabel('Time in sec','fontsize',20)
hold off;
grid;


subplot(3,2,2);hold on;
plot(time,x_traj(2,:),'linewidth',4); 
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta 2','fontsize',20);
xlabel('Time in sec','fontsize',20)
hold off;
grid;

subplot(3,2,3);hold on
plot(time,x_traj(3,:),'linewidth',4); 
plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta 1 dot','fontsize',20)
 xlabel('Time in sec','fontsize',20)
hold off;
grid;

subplot(3,2,4);hold on
plot(time,x_traj(4,:),'linewidth',4); 
plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta 2 dot','fontsize',20)
 xlabel('Time in sec','fontsize',20)
hold off;
grid;

subplot(3,2,5);hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20)
title('Cost','fontsize',20);
   %save('DDP_Data');