clearvars; close all; clc;
set(0, "DefaultFigureWindowStyle", "docked")
%%
clear;close all;clc;

% --- 1. Define System Parameters ---
m = 1;      % Mass (kg)
k = 10;     % Spring stiffness (N/m)
c = 0;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 1000;           % End time (s)
fs_output = 100;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
t_output_points = linspace(t_init, t_final, t_final * fs_output);


% --- 6. Calculate Energy Over Time ---
KE = 0.5 * m * v_history.^2;         % Kinetic Energy
PE = 0.5 * k * x_history.^2;    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy

% --- 7. Plot the Energy Over Time ---
figure;
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_sol, TotalMechanicalEnergy, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Total Mechanical Energy (J)');
title('Energy of Spring-Mass-Damper System Over Time');
grid on;
legend('Total Mechanical Energy');

% % Optional: Plot individual energy components
% nexttile
% plot(t_sol, KE, 'r--', 'LineWidth', 1); hold on;
% plot(t_sol, PE, 'g:', 'LineWidth', 1);
% plot(t_sol, TotalMechanicalEnergy, 'b-', 'LineWidth', 1.5); hold off;
% xlabel('Time (s)');
% ylabel('Energy (J)');
% title('Energy Components of Spring-Mass-Damper System');
% grid on;
% legend('Kinetic Energy', 'Potential Energy', 'Total Mechanical Energy');

nexttile
plot(KE, PE)
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')

%%
clear;close all;clc;
clear integrate_VelocityVerlet integrate_verlet

% --- 1. Define System Parameters ---
m = 1;      % Mass (kg)
k = 10;     % Spring stiffness (N/m)
c = 0.5;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 100;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

x = zeros(length(t_output_points), 6);
v = zeros(length(t_output_points), 6);
time_to_solve = zeros(length(t_output_points), 6);
% --- 3. Define Initial State ---
% The state vector X is [delta_x; vx] (displacement; velocity)
x(1, :) = 0.3165;  % Initial displacement from equilibrium (m)
v(1, :) = 0;         % Initial velocity (m/s)


for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    %%%%%%%%%%%%%%%%%% RK5-Fehlberg Method %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i, 1); v(i, 1)];
    tic;
    X_next = evaluateRKF5(dt, t, X_vec, m, k, c);
    time_to_solve(i, 1) = toc;
    x(i+1, 1) =  X_next(1);
    v(i+1, 1) =  X_next(2);

    %%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i, 2); v(i, 2)];
    tic;
    X_next = evaluateYoshida(dt, t, X_vec, m, k, c);
    time_to_solve(i, 2) = toc;
    x(i+1, 2) =  X_next(1);
    v(i+1, 2) =  X_next(2);

    %%%%%%%%%%%%%%%%%%%% Velocity Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i, 3); v(i, 3)];
    tic;
    X_next = integrate_VelocityVerlet(dt, t, X_vec, m, k, c);
    time_to_solve(i, 3) = toc;
    x(i+1, 3) =  X_next(1);
    v(i+1, 3) =  X_next(2);

    %%%%%%%%%%%%%%%%%%%% Normal Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i+1, 4); v(i+1, 4)];
    tic;
    X_next = integrate_verlet(dt, t, X_vec, m, k, c);
    time_to_solve(i, 4) = toc;
    x(i+1, 4) =  X_next(1);
    v(i+1, 4) =  X_next(2);

    %%%%%%%%%%%%%%%%%%%% RK 4 Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i, 5); v(i, 5)];
    tic;
    X_next = evaluateRK4(dt, t, X_vec, m, k ,c);
    time_to_solve(i, 5) = toc;
    x(i+1, 5) =  X_next(1);
    v(i+1, 5) =  X_next(2);

    %%%%%%%%%%%%%%%%%%%% Eulers Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x(i, 6); v(i, 6)];  
    tic;
    X_next = evaluateEuler(dt, t, X_vec, m, k, c);
    time_to_solve(i, 6) = toc;
    x(i+1, 6) =  X_next(1);
    v(i+1, 6) =  X_next(2);

end
% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
initial_state = [x(1, 1); v(1, 1)];

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% ode45 takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-16, "AbsTol",1e-10,"Stats","on");
[t_sol, X_sol] = ode45(@(t, X) springMassDamperDynamics(t, X, m, k, c), [t_init, t_final], initial_state, options);

% --- 5. Extract Solution Components ---
% x(:, 7) = X_sol(:, 1); % All displacement values over time
% v(:, 7) = X_sol(:, 2);      % All velocity values over time

% --- 6. Calculate Energy Over Time ---
KE = 0.5 * m * v.^2;         % Kinetic Energy
PE = 0.5 * k * x.^2;    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy

KE_ode = 0.5 * m * X_sol(:, 2).^2;         % Kinetic Energy
PE_ode = 0.5 * k * X_sol(:, 1).^2;    % Potential Energy (elastic)
TotalMechanicalEnergy_ode = KE_ode + PE_ode;      % Total Mechanical Energy


% --- 7. Plot the Energy Over Time ---
lgd = {'RK5', 'Yoshida 4', 'Velocity Verlet', 'Normal Verlet', 'RK4', 'Euler', 'ODE45'};
figure;
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, TotalMechanicalEnergy);
hold on
plot(t_sol, TotalMechanicalEnergy_ode, 'k');
xlabel('Time (s)');
ylabel('Total Mechanical Energy (J)');
title('Energy of Spring-Mass-Damper System Over Time');
grid on;
legend(lgd);
ylim([0, 1]);
% % Optional: Plot individual energy components
% nexttile
% plot(t_sol, KE, 'r--', 'LineWidth', 1); hold on;
% plot(t_sol, PE, 'g:', 'LineWidth', 1);
% plot(t_sol, TotalMechanicalEnergy, 'b-', 'LineWidth', 1.5); hold off;
% xlabel('Time (s)');
% ylabel('Energy (J)');
% title('Energy Components of Spring-Mass-Damper System');
% grid on;
% legend('Kinetic Energy', 'Potential Energy', 'Total Mechanical Energy');

nexttile
plot(KE, PE)
hold on
plot(KE_ode, PE_ode, 'k');
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')
legend(lgd);
ylim([0, 1]);
xlim([0, 1]);

nexttile
semilogy(t_output_points, time_to_solve)
xlabel('Time (s)');
ylabel('log_{10}(CPU Time) [s]');
title('CPU Time for each iteration');
grid on;
legend(lgd{1:end-1});

figure
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x)
hold on
plot(t_sol, X_sol(:, 1), 'k');
xlabel('Time [s]')
ylabel('Displacement [m]')
legend(lgd);
xlim([98, 100])

nexttile
plot(t_output_points, v)
hold on
plot(t_sol, X_sol(:, 2), 'k');
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend(lgd);
xlim([98, 100])

%%
clear;close all;clc;

% --- System Parameters ---
L1 = 0.5;  % Length of first link (m)
L2 = 1.0;  % Length of second link (m)
m1 = 0.5;  % Mass of first bob (kg)
m2 = 2.0;  % Mass of second bob (kg)
g = 9.81;  % Acceleration due to gravity (m/s^2)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 10;           % End time (s)
fs_output = 100;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

theta1 = zeros(length(t_output_points), 3);
theta2 = zeros(length(t_output_points), 3);
omega1 = zeros(length(t_output_points), 3);
omega2 = zeros(length(t_output_points), 3);
alpha1 = zeros(length(t_output_points), 3);
alpha2 = zeros(length(t_output_points), 3);

time_to_solve = zeros(length(t_output_points), 3);
% --- 3. Define Initial State ---
% The state vector X is [delta_x; vx] (displacement; velocity)
theta1(1, :) = pi/2;    % Initial angle of first pendulum (90 degrees, horizontal)
theta2(1, :) = pi/2;    % Initial angle of second pendulum (90 degrees, horizontal relative to first)
omega1(1, :) = 0;       % Initial angular velocity of first pendulum
omega2(1, :) = 0;       % Initial angular velocity of second pendulum

initial_state = [theta1(1, 1); theta2(1, 1); omega1(1, 1); omega2(1, 1)];
args = {L1, L2, m1, m2, g};
% func = @(t, X) doublePendulumDynamics(t, X, args);

for i = 1%:length(t_output_points)-1
    t = t_output_points(i);
    %%%%%%%%%%%%%%%%%% RK5-Fehlberg Method %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 1); theta2(i, 1); omega1(i, 1); omega2(i, 1)];
    tic;
    X_next = evaluateRKF5(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 1) = toc;
    theta1(i+1, 1) = X_next(1);
    theta2(i+1, 1) = X_next(2);
    omega1(i+1, 1) = X_next(3);
    omega2(i+1, 1) = X_next(4);

    % %%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%%%
    % X_vec = [theta1(i, 2); theta2(i, 2); omega1(i, 2); omega2(i, 2)];
    % tic;
    % X_next = evaluateYoshida(dt, t, X, L1, L2, m1, m2, g);
    % time_to_solve(i, 2) = toc;
    % theta1(i+1, 2) = X_next(1);
    % theta2(i+1, 2) = X_next(2);
    % omega1(i+1, 2) = X_next(3);
    % omega2(i+1, 2) = X_next(4);

    %%%%%%%%%%%%%%%%%%%% Velocity Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 3); theta2(i, 3); omega1(i, 3); omega2(i, 3);alpha1(i, 3); alpha2(i, 3);];
    tic;
    X_next = integrate_VelocityVerlet(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 3) = toc;
    theta1(i+1, 3) = X_next(1);
    theta2(i+1, 3) = X_next(2);
    omega1(i+1, 3) = X_next(3);
    omega2(i+1, 3) = X_next(4);
    alpha1(i+1, 3) = X_next(5);
    alpha2(i+1, 3) = X_next(6);
    

    % %%%%%%%%%%%%%%%%%%%% Normal Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    % X_vec = [theta1(i, 4); theta2(i, 4); omega1(i, 4); omega2(i, 4)];
    % tic;
    % X_next = integrate_verlet(dt, t, X, L1, L2, m1, m2, g);
    % time_to_solve(i, 4) = toc;
    % theta1(i+1, 4) = X_next(1);
    % theta2(i+1, 4) = X_next(2);
    % omega1(i+1, 4) = X_next(3);
    % omega2(i+1, 4) = X_next(4);

    %%%%%%%%%%%%%%%%%%%% RK 4 Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 2); theta2(i, 2); omega1(i, 2); omega2(i, 2)];
    tic;
    X_next = evaluateRK4(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 2) = toc;
    theta1(i+1, 2) = X_next(1);
    theta2(i+1, 2) = X_next(2);
    omega1(i+1, 2) = X_next(3);
    omega2(i+1, 2) = X_next(4);

    %%%%%%%%%%%%%%%%%%%% Eulers Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 3); theta2(i, 3); omega1(i, 3); omega2(i, 3)];
    tic;
    X_next = evaluateEuler(func, dt, t, X_vec, args);
    time_to_solve(i, 3) = toc;
    theta1(i+1, 3) = X_next(1);
    theta2(i+1, 3) = X_next(2);
    omega1(i+1, 3) = X_next(3);
    omega2(i+1, 3) = X_next(4);

end

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% ode45 takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-16, "AbsTol",1e-10,"Stats","on"); %[t_init t_final]
[t_sol, X_sol] = ode23s(@(t, X) doublePendulumDynamics(t, X, args),t_output_points , initial_state, options);
% --- 5. Extract Solution Components ---
% x(:, 7) = X_sol(:, 1); % All displacement values over time
% v(:, 7) = X_sol(:, 2);      % All velocity values over time


% --- 6. Calculate Energy Over Time ---
KE = 0.5 * ((1 + m2) * omega1.^2 + m2 * (L2 * omega2).^2  + 2 * m2 * L2 * omega1 .* omega2 .* cos(theta1 .* theta2) );         % Kinetic Energy
PE = -g * (cos(theta1) + m2 * L2 * (cos(theta1) + cos(theta2) ) );    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy

theta_1_ode = X_sol(:, 1);
theta_2_ode = X_sol(:, 2);
omega_1_ode = X_sol(:, 3);
omega_2_ode = X_sol(:, 4);

KE_ode = 0.5 * ((1 + m2) * omega_1_ode.^2 + m2 * (L2 * omega_2_ode).^2  + 2 * m2 * L2 * omega_1_ode .* omega_2_ode .* cos(theta_1_ode .* theta_2_ode) );         % Kinetic Energy
PE_ode = -g * (cos(theta_1_ode) + m2 * L2 * (cos(theta_1_ode) + cos(theta_2_ode) ) );    % Potential Energy (elastic)
TotalMechanicalEnergy_ode = KE_ode + PE_ode;      % Total Mechanical Energy

% --- 7. Plot the Energy Over Time ---
lgd = {'RK5', 'RK4', 'Euler', 'ODE23s'};
figure;
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, TotalMechanicalEnergy);
hold on
plot(t_sol, TotalMechanicalEnergy_ode, 'k--');
xlabel('Time (s)');
ylabel('Total Mechanical Energy (J)');
title('Energy of Double Pendulum System Over Time');
grid on;
legend(lgd);
ylim([-50, 200]);

% % Optional: Plot individual energy components
% nexttile
% plot(t_sol, KE, 'r--', 'LineWidth', 1); hold on;
% plot(t_sol, PE, 'g:', 'LineWidth', 1);
% plot(t_sol, TotalMechanicalEnergy, 'b-', 'LineWidth', 1.5); hold off;
% xlabel('Time (s)');
% ylabel('Energy (J)');
% title('Energy Components of Spring-Mass-Damper System');
% grid on;
% legend('Kinetic Energy', 'Potential Energy', 'Total Mechanical Energy');

nexttile
plot(KE, PE, '-.')
hold on
plot(KE_ode, PE_ode, 'k-.');
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')
legend(lgd);
ylim([-60, 60]);
xlim([-5, 300]);

nexttile
semilogy(t_output_points, time_to_solve)
xlabel('Time (s)');
ylabel('log_{10}(CPU Time) [s]');
title('CPU Time for each iteration');
grid on;
legend(lgd{1:end-1});

figure
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, theta1)
hold on
plot(t_sol, X_sol(:, 1), 'k--');
grid on
xlabel('Time [s]')
ylabel('Angular Displacement [rad]')
legend(lgd);
ylim([-5, 5])

nexttile
plot(t_output_points, omega1)
hold on
plot(t_sol, X_sol(:, 3), 'k--');
grid on
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
legend(lgd);
ylim([-20, 20])

%%
X_RKF5 =    [theta1(:, 1), theta2(:, 1)];
X_RK4 =     [theta1(:, 2), theta2(:, 2)];
X_Euler =   [theta1(:, 3), theta2(:, 3)];

plot_pendulums(t_output_points, X_RK4, L1, L2);
           
%%

function plot_pendulums(t_sol, X_sol, L1, L2)
    % --- Optional: Visualize the Pendulum Motion ---
    % This requires converting angles to Cartesian coordinates
    x1 = L1 * sin(X_sol(:, 1));
    y1 = -L1 * cos(X_sol(:, 1)); % Y-axis pointing down
    x2 = x1 + L2 * sin(X_sol(:, 2));
    y2 = y1 - L2 * cos(X_sol(:, 2));
    
    
    % You can also animate it:
    figure;
    h_link1 = plot([0 x1(1)], [0 y1(1)], 'k-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    h_link2 = plot([x1(1) x2(1)], [y1(1) y2(1)], 'k-o', 'LineWidth', 2, 'MarkerSize', 8);
    axis([-3 3 -3 3]); axis equal; grid on;
    title('Double Pendulum Animation');
    pause(0.5)
    for k = 1:5:length(t_sol) % Animate every 5th frame
        set(h_link1, 'XData', [0 x1(k)], 'YData', [0 y1(k)]);
        set(h_link2, 'XData', [x1(k) x2(k)], 'YData', [y1(k) y2(k)]);
        % drawnow;
        pause(0.1)
    end

end

function dX = doublePendulumDynamics(t, X, L1, L2, m1, m2, g)
    % doublePendulumDynamics: Defines the differential equations for a double pendulum.
    %
    % Source: Neumann, E., 2023. myPhysicsLab Double Pendulum 
    % [WWW Document]. URL https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html 
    % (accessed 7.28.25).
    %
    % Inputs:
    %   t:    Current time (scalar) - passed by ODE solver.
    %   X:    Current state vector [theta1; theta2; omega1; omega2]
    %         theta1: Angle of first pendulum (radians, from vertical)
    %         theta2: Angle of second pendulum (radians, from vertical)
    %         omega1: Angular velocity of first pendulum (rad/s)
    %         omega2: Angular velocity of second pendulum (rad/s)
    %   L1:   Length of the first pendulum link (m)
    %   L2:   Length of the second pendulum link (m)
    %   m1:   Mass of the first pendulum bob (kg)
    %   m2:   Mass of the second pendulum bob (kg)
    %   g:    Acceleration due to gravity (m/s^2)
    %
    % Output:
    %   dX:   Time derivatives of the state vector [omega1; omega2; alpha1; alpha2]
    %         theta1_dot = d(theta1)/dt
    %         theta2_dot = d(theta2)/dt
    %         omega1_dot = d(omega1)/dt (angular acceleration of first pendulum)
    %         omega2_dot = d(omega2)/dt (angular acceleration of second pendulum)

    % Extract state variables
    theta1 = X(1);
    theta2 = X(2);
    omega1 = X(3);
    omega2 = X(4);

    theta1_dot = omega1;
    theta2_dot = omega2;

    % Pre-calculate common trigonometric terms for efficiency
    s1 = sin(theta1);
    c1 = cos(theta1);
    s1_2 = sin(theta1 - theta2); % sin(theta1 - theta2)
    c1_2 = cos(theta1 - theta2); % cos(theta1 - theta2)
    den = (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2));

    % Equation for ω1'
    omega1_dot_top_term1 = -g * (2 * m1 + m2) * s1;
    omega1_dot_top_term2 = - m2 * g * sin(theta1 - 2 * theta2);
    omega1_dot_top_term3 = - 2 * s1_2 * m2;
    omega1_dot_top_term4 = (omega2^2 * L2 + omega1^2 * L1 * c1_2);
    
    omega1_dot_top = omega1_dot_top_term1 + omega1_dot_top_term2 + omega1_dot_top_term3 * omega1_dot_top_term4;

    omega1_dot_bot = L1 * den;
    omega1_dot = omega1_dot_top / omega1_dot_bot;

    % Equations for ω2' 
    omega2_dot_top_term1 = 2 * s1_2;
    omega2_dot_top_term2 = omega1^2 * L1 * (m1 + m2);
    omega2_dot_top_term3 = g * (m1 + m2) * c1;
    omega2_dot_top_term4 = omega2^2 * L2 * m2 * c1_2;
    
    omega2_dot_top = omega2_dot_top_term1 * (omega2_dot_top_term2 + omega2_dot_top_term3 + omega2_dot_top_term4);
    omega2_dot_bot = L2 * den;
    omega2_dot = omega2_dot_top / omega2_dot_bot;
    
    % Assemble the derivatives vector dX/dt = [d(theta1)/dt; d(theta2)/dt; d(omega1)/dt; d(omega2)/dt]
    dX = [theta1_dot; theta2_dot; omega1_dot; omega2_dot];
end

function F = F_ext(t)
    if t < 10
        F = 0.1 * sin(2 * pi * t);
    else
        F = 0;
    end
    % F = 0;
end

function dX = springMassDamperDynamics(t, X, m, k, c)

    dx = X(2, :);
    dvx = (F_ext(t) - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end

function X_next = evaluateRKF5(func, dt, t, X_vec, args)

    k1 = dt * func(t, X_vec, args{:});
    k2 = dt * func(t, X_vec + 0.25 *  k1, args{:});
    k3 = dt * func(t, X_vec + 3/32 *  k1 + 9/32 *  k2, args{:});
    k4 = dt * func(t, X_vec + 1932/2197 *  k1 - 7200/2197 *  k2 + 7296/2197 *  k3, args{:});
    k5 = dt * func(t, X_vec + 439/216 *  k1 - 8 *  k2 + 3680/513 *  k3 - 845/4104 *  k4, args{:});
    k6 = dt * func(t, X_vec - 8/27 * k1 + 2 *  k2 - 3544/2565 *  k3 + 1859/4104 *  k4 - 11/40 *  k5, args{:});
    
     X_next =  X_vec + 1/16 *  k1 + 6656/12825 *  k3 + 28561/56430 *  k4 - 9/50 *  k5 + 2/55 *  k6;
end

function X_next = evaluateYoshida(dt, t, X_vec, m, k, c)
    w0 = -(2^(1/3)) / (2 - (2^(1/3)));
    w1 = -(1) / (2 - (2^(1/3)));
    
    c1 = w1 / 2; c4 = c1;
    c2 = (w0 + w1) / 2 ;c3 = c2;
    
    d1 = w1; d3 = d1;
    d2 = w0;

    x = X_vec(1);
    v = X_vec(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = x + c1 .* v .* dt;
    
    % Update acceleration with temporary values
    a = (F_ext(t) - k .* x - c .* v) ./ m;
    
    % Update temporary velocity with new accel
    v = v + d1 .* a .* dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = x + c2 .* v .* dt;
    
    % Update acceleration with temporary values
    a = (F_ext(t) - k .* x - c .* v) ./ m;
    
    % Update temporary velocity with new accel
    v = v + d2 .* a .* dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = x + c3 .* v .* dt;
    
    % Update acceleration with temporary values
    a = (F_ext(t) - k .* x - c .* v) ./ m;
    
    % Update temporary velocity with new accel
    v = v + d3 .* a .* dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update displacement with correctors
    x = x + c4 .* v .* dt;

    X_next = [x;v;a];
end

function acceleration = get_acceleration(func, t_in, X_vec, args)
    % This function adapts the state-space dynamics to return only acceleration.
    % It's designed to be passed as 'func' to integrate_VelocityVerlet_general.
    % Call the original state-space dynamics function
    dX = func(t_in, X_vec, args{:});
    
    [N, ~] = size(dX);

    if mod(N, 2) == 0
        if N > 2
            % Extract the acceleration component (dvx is the second element of dX)
            acceleration = dX(3:end); 
            acceleration = acceleration(:);
        else
            acceleration = dX(2);
        end
    else
        acceleration = [];
        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end
end

function X_next = integrate_VelocityVerlet(func, dt, t, X_vec, varargin)
    % Generalized Velocity Verlet integration method.
    % This function computes one step of the integration.
    %
    % Inputs:
    %   accel_func: Function handle to compute acceleration.
    %               Expected signature: a = accel_func(t, x, v, args_for_accel_func{:})
    %               where x is position and v is velocity (potentially half-step velocity).
    %   dt:         Time step.
    %   t:          Current time.
    %   X_vec:      Current state vector [x_current; v_current; a_old].
    %               x_current is position at time t.
    %               v_current is velocity at time t.
    %               a_old is acceleration at time t (from previous step's calculation).
    %   varargin:   Optional additional arguments to pass directly to accel_func.
    %
    % Output:
    %   X_next:     Next state vector [x_new; v_new; a_new].
    %               x_new is position at time t+dt.
    %               v_new is velocity at time t+dt.
    %               a_new is acceleration at time t+dt.

    % Extract current state variables
    [N, ~] = size(X_vec);
    if mod(N, 2) == 0
        if N == 6
            % Extract the acceleration component (dvx is the second element of dX)
            x_current = X_vec(1:2); 
            x_current = x_current(:);
            v_current = X_vec(3:4); 
            v_current = v_current(:);
            a_old = X_vec(5:6); 
            a_old = a_old(:);
            
        else
            x_current = X_vec(1); % x(t)
            v_current = X_vec(2); % v(t)
            a_old     = X_vec(3); % a(t) - acceleration from the previous step

        end
    else
        x_current = []; % x(t)
        v_current = []; % v(t)
        a_old     = []; % a(t) - acceleration from the previous step

        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end
    % 1. Calculate velocity at half-step
    v_half_step = v_current + a_old * (dt / 2);

    % 2. Update position
    x_new = x_current + v_half_step * dt;

    % 3. Calculate new acceleration using the provided accel_func
    % Pass the new position (x_new) and the half-step velocity (v_half_step)
    % along with any additional parameters to the acceleration function.
    X_vec_for_func = [x_new, v_half_step];
    a_new = get_acceleration(func, t_in, X_vec_for_func, varargin);% a(t + dt

    % 4. Calculate full-step velocity
    v_new = v_half_step + a_new * (dt / 2); % v(t + dt)
    
    % X_next contains [new_position; new_velocity; new_acceleration]
    X_next = [x_new; v_new; a_new]; 
end

function X_next = integrate_verlet(dt, t, X_vec, m, k, c)
    % Position Verlet integration method for dynamics.
    % X_vec should be [current_position; previous_position]
    
    persistent x_previous firstRun
    if isempty(firstRun)
        x_previous = X_vec(1);
        firstRun = 1;
    end

    x_current = X_vec(1); % x(t)
    % x_previous = X_vec(2); % x(t - dt)

    % Estimate velocity for damping force calculation (if damping is desired)
    % This makes it a "velocity-dependent" Verlet, deviating from pure position Verlet.
    v_estimated = (x_current - x_previous) / dt; % Simple central difference approx. v(t - dt/2) or similar

    % 1. Calculate acceleration at current position
    a_current = (F_ext(t) - k * x_current - c * v_estimated) / m; % a(t)

    % 2. Update position
    x_new = 2 * x_current - x_previous + a_current * dt.^2; % x(t + dt)

    % 3. Calculate velocity for output (not used for propagation in core Verlet)
    % A better estimate for v(t+dt) for analysis is often (x(t+dt) - x(t-dt)) / (2*dt)
    v_new_estimated = (x_new - x_previous) / (2 * dt); % v(t + dt)

    % X_next contains [new_position; current_position (for next step's "previous"); estimated_new_velocity]
    x_previous = x_current;
    X_next = [x_new; v_new_estimated]; % Store x_current as x_previous for the next call
end

function X_next = evaluateRK4(func, dt, t, X, args)

    k1 = dt * func(t, X, args{:});
    k2 = dt * func(t, X + 0.5 * k1, args{:});
    k3 = dt * func(t, X + 0.5 * k2, args{:});
    k4 = dt * func(t, X + k3, args{:});
    
    X_next = X + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function X_next = evaluateEuler(func, dt, t, X_vec, args)
    k1 = dt * func(t, X_vec, args{:});
    
    X_next = X_vec + k1;
end