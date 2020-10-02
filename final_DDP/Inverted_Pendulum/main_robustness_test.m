%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% Variables for the inverted pendulum
global g;
global m; 
global l; 
global I; 
global b; 


% Inverted Pendulum Parameters
g = 9.81;       % gravity 
m = 2;          % Mass of the pendulum 
l = 1;          % length of pendulum
b = 0.1;        % damping coefficient
I = m*(l.^2);   % inertia 


% Horizon 
Horizon = 300; % 1.5sec

% Number of Iterations
num_iter = 100;

% Discretization
dt = 0.01;	

% Weight in Final State: (part of terminal cost)
Q_f = zeros(2,2);
Q_f(1,1) = 100; 	%Penalize more w.r.t errors in theta
Q_f(2,2) = 10;      %Penalize more w.r.t errors in theta_dot

%  Weight in the state (for running cost)
P = zeros(2,2); 
P(1,1) = 0;
P(2,2) = 0;

% Weight in the Control: 
R = .1* eye(1,1); 

% Initial Configuration: (Initial state)
xo = zeros(2,1);
xo(1,1) = 0;
xo(2,1) = 0;

% Initial Control:
% Horizon -1 b/c at last time step we have no control
u_k = zeros(1,Horizon-1);   
du_k = zeros(1,Horizon-1);

% Initial trajectory:
x_traj = zeros(2,Horizon);

% Target: (Terminal States)
p_target(1,1) = pi;     % theta
p_target(2,1) = 0;     % theta_dot

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

    % Linearization of the dynamics
    [dfx, dfu] = Jacobians(x_traj(:,j),u_k(:,j));
    
    A(:,:,j) = eye(2) + dfx * dt;    
    B(:,:,j) = dfu * dt;     
    
    % Quadratic expansion of the running cost around the x_traj (nominal trajectory) and u_k (nominal control)
    [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt); 

    L0(j) = dt * l0;            
    Lx(:,j) = dt * l_x;         
    Lxx(:,:,j) = dt * l_xx;     

    Lu(:,j) = dt * l_u;         
    Luu(:,:,j) = dt * l_uu;     
    Lux(:,:,j) = dt * l_ux;    
    
end

%------------------------------------------------> Boundary Conditions
% Initialize value function at the boundary conditions
Vxx(:,:,Horizon)= Q_f;                                                                                                                                                               
Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target);                                       
V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 


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
% dx is initially zero because we start from the same point
dx = zeros(2,1);    

for i=1:(Horizon-1)    
	 du = l_k(:,i) + L_k(:,:,i) * dx;   	% Feedback Controller 
	 dx = A(:,:,i) * dx + B(:,:,i) * du;    % As we propagate forward, we use the linearized dynamics to approximate dx (this is the error from the nominal trajectory)
	 u_new(:,i) = u_k(:,i) + gamma * du;    % Update controls with gamma to prevent controls from updating too fast
end

%Update nominal trajectory (u_k) for new updated controls
u_k = u_new;    

%---------------------------------------------> Simulation of the Nonlinear System
%Create new nominal trajectory based on new control (u_new)
[x_traj] = fnsimulate(xo,u_new,Horizon,dt,0);   
[Cost(:,k)] =  fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
x1(k,:) = x_traj(1,:);
 

fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end

time(1)=0;
for i= 2:Horizon
	time(i) =time(i-1) + dt;  
end
%%
%---------------------------------------------> Implementation of Robustness Test Algorithm
x_star = x_traj; 
u_star = u_new; 
K_star = L_k;

x_new = zeros(2,Horizon);
u = zeros(1,Horizon-1);
num_traj = 20; 

for i = 1:num_traj
    
    [x_new,u_neww] = fnsimulate_noise(x_new,u_star,x_star,K_star,Horizon,dt,1.0);
    [Cost(:,k)] =  fnCostComputation(x_new,u_neww,p_target,dt,Q_f,R);

    figure(1);

    subplot(2,2,1);
    set(gca,'FontSize',32)

    hold on;
    plot(time,x_new(1,:),'linewidth',4);  
    plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4);
    title('$\theta$','Interpreter','latex','fontsize',32);
    xlabel('Time in sec','fontsize',32);
    ylabel('Rad/s','fontsize',32);
    hold off;
    grid;   

    subplot(2,2,2);
    set(gca,'FontSize',32)
    hold on;
    plot(time,x_new(2,:),'linewidth',4); 
    plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4);
    title('$\dot{\theta}$','Interpreter','latex','fontsize',32);
    xlabel('Time in sec','fontsize',32);
    ylabel('Rad/s','fontsize',32);

    hold off;
    grid;
    
    subplot(2,2,3);
    set(gca,'FontSize',32)
    hold on
    plot(Cost,'linewidth',2); 
    xlabel('Iterations','fontsize',32);
    title('Cost','fontsize',32);
    
    subplot(2,2,4);
    set(gca,'FontSize',32)
    hold on;
    title('Animation','fontsize',32);
    O = [0 0];
    axis(gca,'equal');
    axis([-1.5 1.5 -1.5 1.5]);
    P=1*[sin(x_new(1,end)) -cos(x_new(1,end))];
    plot([O(1) P(1)],[O(2) P(2)], 'LineWidth', 4);
    grid;
end 









