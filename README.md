# Differential Dynamic Programming (DDP)
This repository implements Differential Dynamic Programming as a optimal control algorithm for three dynamical systems: the simple inverted pendulum, the cartpole system, and a two link arm. 

## Demo of Inverted Pendulum
  ![Demo](/media/InvPend.gif)

## State Convergence Graphs for Inverted Pendulum
  ![Convergence Graph](/media/InvPendGraph.png)

## Demo of the CartPole System
  ![Demo](/media/CartPoleVideo.gif)

## State Convergence Graphs for CartPole System
  ![Convergence Graph](/media/CartPole.png)

## Testing the Robustness of DDP on Inverted Pendulum

To test the robustness of our DDP policy against stochastic forces that act as disturbances in our dynamics, we first run the DDP optimization until it converges. This results in an optimal trajectory, optimal control input, and optimal feeback gains. 
We then initialize a new trajectory and propagate the dynamics forward, perturbed by a stochastic force term, using the controller detailed in the following algorithm. We run this for 20 trajectories in the results below. 
	![Robustness Test Algorithm](/media/robustness_test_algorithm.png)
	![Robustness Test Results](/media/robustness_test_results.png)