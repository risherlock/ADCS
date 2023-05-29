% Numerical simulation of simple pendulum
% Rishav (2023-05-28)

clc
clear
close all

% Simulation params
start_time = 0;
stop_time = 50;
dt = 0.05;
time = start_time:dt:stop_time;

% Control params
K = 0.01;

% Moment of inertial matrix
I = diag([0.05, 0.05, 0.05]);

% Initial conditions
q0 = [1, 0, 0, 0]';    % [q0, q1, q2, q3]
w0 = [0.2, 0.3, 0.4]'; % Rate, rad/s

% Memory allocations
state = zeros(7, length(time));
state(:,1) = [q0; w0];

% Numerical integration
for t = 1:length(time)-1
    % Control law
    w = state(5:7,t);
    tau = -K * w;
    
    % Numerical integration
    fn = @(y)satellite(y, I, tau);
    state(:,t+1) = RK4(fn, state(:,t), dt);
    
    % Normalize quaternion
    state(1:4,t+1) = state(1:4,t+1) ./ norm(state(1:4,t+1));
end

% Plot
subplot(4,1,1);
plot(time,state(1,:), 'lineWidth', 1.5);
xlabel('Time (s)'); ylabel("q_{0}"); grid on;
subplot(4,1,2);
plot(time,state(2,:), 'lineWidth', 1.5);
xlabel('Time (s)'); ylabel("q_{1}"); grid on;
subplot(4,1,3);
plot(time,state(3,:), 'lineWidth', 1.5);
xlabel('Time (s)'); ylabel("q_{2}"); grid on;
subplot(4,1,4);
plot(time,state(4,:), 'lineWidth', 1.5);
xlabel('Time (s)'); ylabel("q_{3}"); grid on;
suptitle("Quaternion profile");

figure
subplot(3,1,1);
plot(time,state(5,:), 'lineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel("\omega_{1} (rad/s)");
subplot(3,1,2);
plot(time,state(6,:), 'lineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel("\omega_{1} (rad/s)");
subplot(3,1,3);
plot(time,state(7,:), 'lineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel("\omega_{1} (rad/s)");
suptitle("Angular rate profile");

%%%%%%%%%%%%%
% Functions %
%%%%%%%%%%%%%

% Rotational dynamics of satellite
function [state_dot] = satellite(state, I, tau)
q = state(1:4);
w = state(5:7);

% Dynamics matrix
M = [   0, -w(1), -w(2), -w(3); 
     w(1),     0,  w(3), -w(2); 
     w(2), -w(3),     0,  w(1);
     w(3),  w(2), -w(1),    0];
       
q_dot = 0.5 * M * q;                  % Kinematics
w_dot = I \ (-cross(w, I * w) + tau); % Kinetics

state_dot = [q_dot; w_dot];
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
