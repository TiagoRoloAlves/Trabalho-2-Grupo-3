% Mars Lander + Rover Nonlinear Simulation
clc; clear; close all;

% Constants
gMars = 3.73; % Mars gravity (m/s^2)
rhoMars = 0.02; % Mars atmospheric density (kg/m^3)
Cd = 1.05; % Drag coefficient
Am = 9.5; % Approximate average cross-sectional area (m^2)
mCombined = 1803.24 + 1025; % Combined mass of lander + rover (kg)
T = 1000; % Thrust of each retrorocket (N)
alpha = deg2rad(0); % Angle of retrorockets (radians)
zcg = -0.26; % z-coordinate of thrusters relative to combined CoG (m)
J_combined = diag([4445.10, 3775.64, 6105.02]); % Combined moment of inertia

% Initial Conditions
p0 = [0; 0; 2000]; % Initial position (m) [x, y, z]
v0 = [0; 0; 0]; % Initial velocity (m/s) [vx, vy, vz]
lambda0 = [0; 0; 0]; % Initial Euler angles (radians) [phi, theta, psi]
omega0 = [0; 0; 0]; % Initial angular velocity (rad/s) [wx, wy, wz]
state0 = [p0; v0; lambda0; omega0]; % Initial state vector

% Simulation Parameters
tSpan = [0, 50]; % Time span for simulation (s)
dt = 0.01; % Time step (s)

% Position vectors of the retrorockets relative to CoM
p1 = [1.5; -1.9; zcg];
p2 = [-1.5; -1.9; zcg];
p3 = [-1.5; 1.9; zcg];
p4 = [1.5; 1.9; zcg];

% Dynamics function
dynamics = @(t, state) nonlinear_dynamics(state, T, alpha, mCombined, ...
    J_combined, gMars, rhoMars, Cd, Am, p1, p2, p3, p4);

% Solve ODE
[t, sol] = ode45(dynamics, tSpan, state0);

% Extract states
pos = sol(:, 1:3); % Position [x, y, z]
vel = sol(:, 4:6); % Velocity [vx, vy, vz]

% Plot results
figure;
subplot(2, 1, 1);
plot(t, pos, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');
title('Position States');
grid on;

subplot(2, 1, 2);
plot(t, vel, '--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');
title('Velocity States');
grid on;

function dx = nonlinear_dynamics(state, T, alpha, m, J, gMars, ...
    rho, Cd, A, p1, p2, p3, p4)
    % Extract states
    p = state(1:3); % Position
    v = state(4:6); % Velocity
    lambda = state(7:9); % Euler angles
    omega = state(10:12); % Angular velocity
    
    % Rotation matrix from body to inertial frame
    R = eul2rotm(lambda', 'ZYX');
    
    % Forces
    Fg = R' * [0; 0; -m * gMars]; % Gravity in body frame
    Fd = -0.5 * rho * Cd * A * norm(v) * v; % Drag force
    Ft = R' * (thrust_forces(T, alpha)); % Thrust forces
    
    % Total force
    total_force = Ft + Fd + Fg;
    
    % Moments
    moments = cross(p1, thrust_forces(T, alpha)) + ...
              cross(p2, thrust_forces(T, alpha)) + ...
              cross(p3, thrust_forces(T, alpha)) + ...
              cross(p4, thrust_forces(T, alpha));
          
    % Equations of motion
    dp = R * v; % Translational kinematics
    dv = total_force / m; % Translational dynamics
    dlambda = omega; % Euler kinematics (approximation)
    domega = J \ (moments - cross(omega, J * omega)); % Rotational dynamics
    
    % State derivative
    dx = [dp; dv; dlambda; domega];
end

function Ft = thrust_forces(T, alpha)
    % Calculate thrust forces
    Ft = [0; sin(alpha) * T; -cos(alpha) * T];
end
