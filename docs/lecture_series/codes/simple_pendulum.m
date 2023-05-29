% Numerical simulation of simple pendulum
% Rishav (2023-05-28)

clc
clear
close all

% Simulation params
start_time = 0;
stop_time = 25;
dt = 0.05;
time = start_time:dt:stop_time;

% System params
L = 3;     % Length of pendulum
d = 0.004; % Damping coefficient
g = 9.8;   % Gravity
m = 0.01;  % Mass

% Initial conditions
theta = 0.01;    % Angular position
theta_dot = 0.0; % Angular velocity

% Memory allocations
state = zeros(2, length(time));
state(:,1) = [theta, theta_dot]';

% Numerical integration
for t = 1:length(time)-1
    fn = @(y)pendulum(y, L, d, g, m);
    state(:,t+1) = RK4(fn, state(:,t), dt);
end

% Plot
plot(time,state(1,:), 'lineWidth', 2); hold on; 
plot(time,state(2,:), 'lineWidth', 2);
xlabel('Time'); 
ylabel('State');
legend('\theta','d\theta/dt');
title('Simple Pendulum');
grid on;

%%%%%%%%%%%%%
% Functions %
%%%%%%%%%%%%%

% Simple pendulum differential equation
function [state_dot] = f(state, L, d, g, m)
    theta = state(1);
    theta_dot = state(2);
    
    theta_dot_dot = - (g / L) * sin(theta) / L  - (d / m) * theta_dot;
    state_dot = [theta_dot, theta_dot_dot]';
end

% RK4 integrator
function y_update = RK4(f,y,dt)
k1 = f(y);
k2 = f(y + 0.5 * dt * k1);
k3 = f(y + 0.5 * dt * k2);
k4 = f(y + dt * k3);
    
K = (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4); 
y_update = y + K * dt; 
end
