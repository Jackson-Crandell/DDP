%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803  Fall  2018                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Evangelos Theodorou                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% Variables for the inverted pendulum
global g;
global m; 
global l; 
global I; 
global b; 


% Sample parameters of choice as used in
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling

g = 9.81;       % gravity 
m = 2;        % Mass of the pendulum 
l = 0.3;        % length of pendulum
b = 0.1;        % damping coefficient
I = m*(l.^2);   % inertia 


% Horizon 
% Best: 700
Horizon = 800; % 1.5sec

% Number of Iterations
% Best: 200
num_iter = 250;

% Discretization
dt = 0.001;	% .01 * 300 = 3 seconds


% Weight in Final State: (part of terminal cost)
Q_f = zeros(2,2);
Q_f(1,1) = 100; 	%Penalize more w.r.t errors in theta
Q_f(2,2) = 10;     %Penalize more w.r.t errors in theta_dot


% Weight in the Control: 
% Best: 1
R = 1 * eye(1,1); % Weight control equally


% Initial Configuration: (Initial state)
xo = zeros(2,1);
xo(1,1) = 0;
xo(2,1) = 0;


% Initial Control:
u_k = zeros(1,Horizon-1);   % Horizon -1 b/c at last time step we have no control
du_k = zeros(1,Horizon-1);


% Initial trajectory:
x_traj = zeros(2,Horizon);

% Target: (Terminal States)
p_target(1,1) = pi;     % theta
p_target(2,1) = 0;      % theta_dot

% Learning Rate .5
gamma = 0.5;

%Initialize Q Value Function
Q = zeros(1,Horizon);
Q_x = zeros(2,Horizon);
Q_u = zeros(1,Horizon);
Q_xx = zeros(2,2,Horizon);
Q_uu = zeros(1,1,Horizon);
Q_ux = zeros(1,2,Horizon);
 
 
for k = 1:num_iter % Run for a certain number of iterations

%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
for  j = 1:(Horizon-1) %Discretize trajectory for each timestep

	% linearization of dynamics (Jacobians dfx and dfu)
	[dfx, dfu] = Jacobians(x_traj(:,j),u_k(:,j));

	% Quadratic expansion of the running cost around the x_trajectory (nominal) and u_k which is the nominal control
	[l0,l_x,l_xx,l_u,l_uu,l_ux] = rCost(x_traj(:,j), u_k(:,j), j,R,dt);    % for each time step compute the cost

	q0(j) = dt * l0;            % zero order term (scalar)
	q_k(:,j) = dt * l_x;        % gradient of running cost w.r.t x (vector)
	Q_k(:,:,j) = dt * l_xx;     % Hessian of running cost w.r.t x (matrix)
	r_k(:,j) = dt * l_u;        % gradient of running cost w.r.t u (vector)
	R_k(:,:,j) = dt * l_uu;     % Hessian of running cost w.r.t u (matrix)
	P_k(:,:,j) = dt * l_ux;     % Hessian of running cost w.r.t ux (matrix)

	A(:,:,j) = eye(1,1) + dfx * dt;     % This is PHI in notes (Identity matrix) + gradient of dynamics w.r.t x * dt
	B(:,:,j) = dfu * dt;                % B matrix in notes is the linearized contols
end

%------------------------------------------------> Boundary Conditions
% initialize value function 
Vxx(:,:,Horizon)= Q_f;                                  % Initialize Hessian of value function (Matrix)                                                                                         
Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target);   % Gradient of value function (Vector)
V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);  %Value function (scalar)


%------------------------------------------------> Backpropagation of the Value Function
for j = (Horizon-1):-1:1
		 
	 Q(:,j) = l0 + V(:,j+1);
     Q_x(:,j) = l_x + A(:,:,j)'*Vx(:,j+1);
     Q_u(:,j) = l_u + B(:,:,j)'*Vx(:,j+1);
     Q_xx(:,:,j) = l_xx + A(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);
     Q_uu(:,:,j) = l_uu + B(:,:,j)'*Vxx(:,:,j+1)*B(:,:,j);
     Q_ux(:,:,j) = l_ux + B(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);
     
	 L_k(:,:,j)= - inv(Q_uu(:,:,j))*Q_ux(:,:,j);   % Feedback term
	 l_k (:,j) = - inv(Q_uu(:,:,j))*Q_u(:,j);    % Feedforward term
	 
	 Vxx(:,:,j) = Q_xx(:,:,j) - Q_ux(:,:,j)'*inv(Q_uu(:,:,j))*Q_ux(:,:,j);
	 Vx(:,j)= Q_x(:,j) - Q_ux(:,:,j)'*inv(Q_uu(:,:,j))*Q_u(:,j);
	 V(:,j) = Q(:,j) - 0.5*Q_u(:,j)'*inv(Q_uu(:,:,j))*Q_u(:,j);
end 


%----------------------------------------------> Find the controls
dx = zeros(2,1);    % dx is initially zero because we start from the same point

for i=1:(Horizon-1)    
	 du = l_k(:,i) + L_k(:,:,i) * dx;   	%Feedback Controller 
	 dx = A(:,:,i) * dx + B(:,:,i) * du;    %As we propagate forward, we use the linearized dynamics to approximate dx (this is the error from the nominal trajectory)
	 u_new(:,i) = u_k(:,i) + gamma * du;    %Update controls with gamma to prevent controls from updating too fast
end

u_k = u_new;    %Update nominal trajectory (u_k) for new updated controls


%---------------------------------------------> Simulation of the Nonlinear System
[x_traj] = simulate(xo,u_new,Horizon,dt,0);   %Create new nominal trajectory based on new control (u_new)
[Cost(:,k)] =  CostComputation(x_traj,u_k,p_target,dt,Q_f,R);
x1(k,:) = x_traj(1,:);
 

fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end

time(1)=0;
for i= 2:Horizon
	time(i) =time(i-1) + dt;  
end


figure(1);
subplot(2,2,1)
hold on
plot(time,x_traj(1,:),'linewidth',4);  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta 1','fontsize',20); 
xlabel('Time in sec','fontsize',20)
hold off;
grid;

subplot(2,2,2);hold on;
plot(time,x_traj(2,:),'linewidth',4); 
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
title('Theta dot','fontsize',20);
xlabel('Time in sec','fontsize',20)
hold off;
grid;

subplot(2,2,3);hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20)
title('Cost','fontsize',20);
%save('DDP_Data');