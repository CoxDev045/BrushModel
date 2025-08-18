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
m = 7.64e-10;      % Mass (kg)
k = 0.37;     % Spring stiffness (N/m)
c = 1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 11;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

x1 = zeros(length(t_output_points), 10);
v1 = zeros(length(t_output_points), 10);
a1 = zeros(length(t_output_points), 7);
e1 = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), 10);

% --- 3. Define Initial State ---
% The state vector X is [delta_x; vx] (displacement; velocity)
x1(1, :) = 0.1;    % Initial angle of first pendulum (90 degrees, horizontal)
v1(1, :) = 0.0;       % Initial angular velocity of first pendulum

initial_state = [x1(1, 1); v1(1, 1)];
args = {m, k, c};

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% odeXX takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","on"); 
tic;
[t_sol, X_sol] = ode23t(@(t, X) springMassDamperDynamics(t, X, args{:}),t_output_points , initial_state, options);
toc


% For first iteration of adaptive RK(4)5
t_current = t_init;
h_current = dt;

options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off"); 
for i = 1:length(t_output_points)-1
    t = t_output_points(i);

    %%%%%%%%%%%%%%%%%% Adaptive RK45-Sarafyan Method %%%%%%%%%%%%%%%%%%%%%%%%
    % Advance the solution until we reach the target time
    y_current = [x1(i, 9); v1(i, 9)];
    t_target = t_output_points(i+1);
    tic;
    while t_current < t_target
        % Call the adaptive step function
        [y_current, h_current] = adaptive_ODE(@springMassDamperDynamics, h_current, t_current,t_final, y_current, args);

        % Update current time based on the step taken
        t_current = t_current + h_current;
    end
    time_to_solve(i, 9) = toc;
    x1(i+1, 9) = y_current(1);
    v1(i+1, 9) = y_current(2);
    e1(i+1, 9) = norm(X_sol(i, :).' - y_current, inf);

    %%%%%%%%%%%%%%%%%% RK5-Fehlberg Method %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1)];
    tic;
    X_next = evaluateRKF5(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    e1(i+1, 1) = norm(X_sol(i, :).' - X_next);

    %%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2); a1(i, 2);];
    tic;
    X_next = evaluateYoshida(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    a1(i+1, 2) = X_next(3);
    e1(i+1, 2) = norm(X_sol(i, :) - X_next(1:2));

    % %%%%%%%%%%%%%%%%%%%% Velocity Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 3); v1(i, 3); a1(i, 3);];
    tic;
    X_next = integrate_VelocityVerlet(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 3) = toc;
    x1(i+1, 3) = X_next(1);
    v1(i+1, 3) = X_next(2);
    a1(i+1, 3) = X_next(3);
    e1(i+1, 3) = norm(X_sol(i, :) - X_next(1:2));

    %%%%%%%%%%%%%%%%%%%% Normal Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    if i~= 1
        X_vec = [x1(i, 4); v1(i, 4); a1(i, 4);];
        tic;
        X_next = integrate_verlet(@springMassDamperDynamics, dt, t, X_vec, args);
        time_to_solve(i, 4) = toc;
        x1(i+1, 4) = X_next(1);
        v1(i+1, 4) = X_next(2);
        a1(i+1, 4) = X_next(3);
        e1(i+1, 4) = norm(X_sol(i, :) - X_next(1:2));

    else
        X_vec = [x1(1, 4); v1(1, 4);];
        X_next = evaluateEuler(@springMassDamperDynamics, dt, t_output_points(1), X_vec, args);
        a1(i+1, 4) = x1(1, 4);
        x1(i+1, 4) = X_next(1);
        v1(i+1, 4) = X_next(2);
        e1(i+1, 4) = norm(X_sol(i, :).' - X_next);

    end

    %%%%%%%%%%%%%%%%%%%% RK 4 Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 5); v1(i, 5);];
    tic;
    X_next = evaluateRK4(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 5) = toc;
    x1(i+1, 5) = X_next(1);
    v1(i+1, 5) = X_next(2);
    e1(i+1, 5) = norm(X_sol(i, :).' - X_next);

    %%%%%%%%%%%%%%%%%%%% Eulers Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 6); v1(i, 6);];
    tic;
    X_next = evaluateEuler(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 6) = toc;
    x1(i+1, 6) = X_next(1);
    v1(i+1, 6) = X_next(2);
    e1(i+1, 6) = norm(X_sol(i, :).' - X_next);


    %%%%%%%%%%%%%%% Implicit Euler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 7); v1(i, 7);];
    tic;
    X_next = evaluateImplicitEuler(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 7) = toc;
    x1(i+1, 7) = X_next(1);
    v1(i+1, 7) = X_next(2);
    e1(i+1, 7) = norm(X_sol(i, :).' - X_next);


    %%%%%%%%%%%%%% TR-BDF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 8); v1(i, 8);];
    tic;
    X_next = trbdf2_step(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 8) = toc;
    x1(i+1, 8) = X_next(1);
    v1(i+1, 8) = X_next(2);
    e1(i+1, 8) = norm(X_sol(i, :).' - X_next);

    %%%%%%%%% ODE23t at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 10); v1(i, 10);];
    tic;
    [~, X_next] = ode23t(@(t, X) springMassDamperDynamics(t, X, args{:}),[t, t_target] , X_vec, options);
    time_to_solve(i, 10) = toc;
    x1(i+1, 10) = X_next(end, 1);
    v1(i+1, 10) = X_next(end, 2);
    e1(i+1, 10) = norm(X_sol(i, :) - X_next(end, :));

end
%%
x_adapt = [];
v_adapt = [];

x_adapt(1) = x1(1, 1);    % Initial angle of first pendulum (90 degrees, horizontal)
v_adapt(1) = v1(1, 1);       % Initial angular velocity of first pendulum

t_adapt = [];
t_adapt(1) = t_init;

time_to_solve_adapt = [];
time_to_solve_adapt(1) = 1;

i = 1;
while t_adapt(i) < t_final
    %%%%%%%%%%%%%%%%%%%% Adaptive time stepping Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x_adapt(i); v_adapt(i);];
    if i == 1
        dt_next = dt;
    end
    tic;
    [X_next, dt_next] = adaptive_ODE12t(@springMassDamperDynamics, dt_next, t_adapt(i), X_vec, args);
    time_to_solve_adapt(i+1) = toc;
    x_adapt(i+1) = X_next(1);
    v_adapt(i+1) = X_next(2);
    t_adapt(i+1) = t_adapt(i) + dt_next;
    i = i + 1;

end

%%
% --- 5. Extract Solution Components ---
x_1_ode = X_sol(:, 1); % All displacement values over time
v_1_ode = X_sol(:, 2);      % All velocity values over time

% --- 6. Calculate Energy Over Time ---
p = m * v1;

KE = 0.5 * m * v1.^2;         % Kinetic Energy
PE = 0.5 * k * x1.^2;    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy
Lagrangian = KE - PE;
Hamiltonian = p.^2 / (2 * m) + PE;


p_ode = m * v_1_ode;

KE_ode = 0.5 * m * v_1_ode.^2;         % Kinetic Energy
PE_ode = 0.5 * k * x_1_ode.^2;    % Potential Energy (elastic)
TotalMechanicalEnergy_ode = KE_ode + PE_ode;      % Total Mechanical Energy
Lagrangian_ode = KE_ode - PE_ode;
Hamiltonian_ode = p_ode.^2 / (2 * m) + PE_ode;

% For adaptive method
p_adapt = m * v_adapt;

KE_adapt = 0.5 * m * v_adapt.^2;         % Kinetic Energy
PE_adapt = 0.5 * k * x_adapt.^2;    % Potential Energy (elastic)
TotalMechanicalEnergy_adapt = KE_adapt + PE_adapt;      % Total Mechanical Energy
Lagrangian_adapt = KE_adapt - PE_adapt;
Hamiltonian_adapt= p_adapt.^2 / (2 * m) + PE_adapt;

%%
% --- 7. Plot the Energy Over Time ---
% lgd = {'RK5','Yoshida 4', 'Velocity Verlet', 'Normal Verlet', 'RK4', 'Euler', 'Implicit Euler', 'TR-BDF2', 'adaptiveRK45', 'adaptive12t',  'ODE23t'};
% lgd = {'Implicit Euler', 'TR-BDF2', 'adaptiveRK45', 'ODE23s at each time step', 'adaptive12t' ,  'ODE23s'};
lgd = {'Implicit Euler', 'TR-BDF2', 'adaptiveRK45', 'ODE23s at each time step' ,  'ODE23s'};


figure;
T = tiledlayout(2, 2);
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, TotalMechanicalEnergy(:, end-3:end));
hold on
% plot(t_adapt, TotalMechanicalEnergy_adapt, 'r--');
plot(t_sol, TotalMechanicalEnergy_ode, 'k--');
xlabel('Time (s)');
ylabel('Total Mechanical Energy (J)');
title('Energy of SDOF System Over Time');
grid on;
legend(lgd);
% ylim([-0.02, 0.02]);

nexttile(2,[2, 1])
loglog(KE(:, end-3:end), PE(:, end-3:end), '-.')
hold on
% loglog((KE_adapt), (PE_adapt), 'r-.');
loglog((KE_ode), (PE_ode), 'k-.');
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')
legend(lgd);
% axis([0, 0.02, 0, 0.02]);

nexttile
semilogy(t_output_points, time_to_solve(:, end-3:end))
hold on
% semilogy(t_adapt, time_to_solve_adapt)
xlabel('Time (s)');
ylabel('log_{10}(CPU Time) [s]');
title('CPU Time for each iteration');
grid on;
legend(lgd{1:end-1});
% ylim([-1, 1]);

figure
T = tiledlayout('vertical');
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1(:, end-3:end))
hold on
% plot(t_adapt, x_adapt, 'r--');
plot(t_sol, X_sol(:, 1), 'k--');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
legend(lgd);
ylim([-1, 1]);

nexttile
plot(t_output_points, v1(:, end-3:end))
hold on
% plot(t_adapt, v_adapt, 'r--');
plot(t_sol, X_sol(:, 2), 'k--');
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend(lgd);
ylim([-2, 2]);

nexttile
semilogy(t_output_points, e1(:, end-3:end))
grid on
xlabel('Time [s]')
ylabel('log_{10}(Absolute Error) [m/s]')
legend(lgd{1:end-1});
 
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
t_final = 4;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

theta1 = zeros(length(t_output_points), 8);
theta2 = zeros(length(t_output_points), 8);
omega1 = zeros(length(t_output_points), 8);
omega2 = zeros(length(t_output_points), 8);
alpha1 = zeros(length(t_output_points), 8);
alpha2 = zeros(length(t_output_points), 8);

time_to_solve = zeros(length(t_output_points), 3);
% --- 3. Define Initial State ---
% The state vector X is [delta_x; vx] (displacement; velocity)
theta1(1, :) = pi/2;    % Initial angle of first pendulum (90 degrees, horizontal)
theta2(1, :) = pi/2;    % Initial angle of second pendulum (90 degrees, horizontal relative to first)
omega1(1, :) = 0;       % Initial angular velocity of first pendulum
omega2(1, :) = 0;       % Initial angular velocity of second pendulum

initial_state = [theta1(1, 1); theta2(1, 1); omega1(1, 1); omega2(1, 1)];
args = {L1, L2, m1, m2, g};

X_vec = [theta1(1, 4); theta2(1, 4); omega1(1, 4); omega2(1, 4);];
X_next_verlet = evaluateEuler(@doublePendulumDynamics, dt, t_output_points(1), X_vec, args);
alpha1(1, 4) = theta1(1, 4);
alpha2(1, 4) = theta2(1, 4);
theta1(1, 4) = X_next_verlet(1);
theta2(1, 4) = X_next_verlet(2);
omega1(1, 4) = X_next_verlet(3);
omega2(1, 4) = X_next_verlet(4);

for i = 1:length(t_output_points)-1
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

    %%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 2); theta2(i, 2); omega1(i, 2); omega2(i, 2);alpha1(i, 2); alpha2(i, 2);];
    tic;
    X_next = evaluateYoshida(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 2) = toc;
    theta1(i+1, 2) = X_next(1);
    theta2(i+1, 2) = X_next(2);
    omega1(i+1, 2) = X_next(3);
    omega2(i+1, 2) = X_next(4);
    alpha1(i+1, 2) = X_next(5);
    alpha2(i+1, 2) = X_next(6);

    % %%%%%%%%%%%%%%%%%%%% Velocity Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 3); theta2(i, 3); omega1(i, 3); omega2(i, 3); alpha1(i, 3); alpha2(i, 3);];
    tic;
    X_next = integrate_VelocityVerlet(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 3) = toc;
    theta1(i+1, 3) = X_next(1);
    theta2(i+1, 3) = X_next(2);
    omega1(i+1, 3) = X_next(3);
    omega2(i+1, 3) = X_next(4);
    alpha1(i+1, 3) = X_next(5);
    alpha2(i+1, 3) = X_next(6);
   

    %%%%%%%%%%%%%%%%%%%% Normal Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 4); theta2(i, 4); omega1(i, 4); omega2(i, 4); alpha1(i, 4); alpha2(i, 4);];
    tic;
    X_next = integrate_verlet(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 4) = toc;
    theta1(i+1, 4) = X_next(1);
    theta2(i+1, 4) = X_next(2);
    omega1(i+1, 4) = X_next(3);
    omega2(i+1, 4) = X_next(4);
    alpha1(i+1, 4) = X_next(5);
    alpha2(i+1, 4) = X_next(6);

    %%%%%%%%%%%%%%%%%%%% RK 4 Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 5); theta2(i, 5); omega1(i, 5); omega2(i, 5)];
    tic;
    X_next = evaluateRK4(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 5) = toc;
    theta1(i+1, 5) = X_next(1);
    theta2(i+1, 5) = X_next(2);
    omega1(i+1, 5) = X_next(3);
    omega2(i+1, 5) = X_next(4);

    %%%%%%%%%%%%%%%%%%%% Eulers Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 6); theta2(i, 6); omega1(i, 6); omega2(i, 6)];
    tic;
    X_next = evaluateEuler(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 6) = toc;
    theta1(i+1, 6) = X_next(1);
    theta2(i+1, 6) = X_next(2);
    omega1(i+1, 6) = X_next(3);
    omega2(i+1, 6) = X_next(4);

    %%%%%%%%%%%%%%% Implicit Euler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 7); theta2(i, 7); omega1(i, 7); omega2(i, 7)];
    tic;
    X_next = evaluateImplicitEuler(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 7) = toc;
    theta1(i+1, 7) = X_next(1);
    theta2(i+1, 7) = X_next(2);
    omega1(i+1, 7) = X_next(3);
    omega2(i+1, 7) = X_next(4);

    %%%%%%%%%%%%%%% TR-BDF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [theta1(i, 8); theta2(i, 8); omega1(i, 8); omega2(i, 8)];
    tic;
    X_next = trbdf2_step(@doublePendulumDynamics, dt, t, X_vec, args);
    time_to_solve(i, 8) = toc;
    theta1(i+1, 8) = X_next(1);
    theta2(i+1, 8) = X_next(2);
    omega1(i+1, 8) = X_next(3);
    omega2(i+1, 8) = X_next(4);

end

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% ode45 takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-16, "AbsTol",1e-10,"Stats","on"); %[t_init t_final]
[t_sol, X_sol] = ode23t(@(t, X) doublePendulumDynamics(t, X, args{:}),t_output_points , initial_state, options);
% --- 5. Extract Solution Components ---
% x(:, 7) = X_sol(:, 1); % All displacement values over time
% v(:, 7) = X_sol(:, 2);      % All velocity values over time


% --- 6. Calculate Energy Over Time ---
p1 = (1 + m2) * omega1 + m2 * L2 * omega2 .* cos(theta1 - theta2);
p2 = m2 * L2^2 * omega2 + m2 * L2 * omega1 .* cos(theta1 - theta2);

KE = 0.5 * ((1 + m2) * omega1.^2 + m2 * (L2 * omega2).^2  + 2 * m2 * L2 * omega1 .* omega2 .* cos(theta1 .* theta2) );         % Kinetic Energy
PE = -g * (cos(theta1) + m2 * L2 * (cos(theta1) + cos(theta2) ) );    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy
Lagrangian = KE - PE;
Hamiltonian = p1 .* omega1 + p2 .* omega2 - Lagrangian;

% For MATLAB ODE Solver
theta_1_ode = X_sol(:, 1);
theta_2_ode = X_sol(:, 2);
omega_1_ode = X_sol(:, 3);
omega_2_ode = X_sol(:, 4);


p1_ode = (1 + m2) * omega_1_ode + m2 * L2 * omega_2_ode .* cos(theta_1_ode - theta_2_ode);
p2_ode = m2 * L2^2 * omega_2_ode + m2 * L2 * omega_1_ode .* cos(theta_1_ode - theta_2_ode);

KE_ode = 0.5 * ((1 + m2) * omega_1_ode.^2 + m2 * (L2 * omega_2_ode).^2  + 2 * m2 * L2 * omega_1_ode .* omega_2_ode .* cos(theta_1_ode .* theta_2_ode) );         % Kinetic Energy
PE_ode = -g * (cos(theta_1_ode) + m2 * L2 * (cos(theta_1_ode) + cos(theta_2_ode) ) );    % Potential Energy (elastic)
TotalMechanicalEnergy_ode = KE_ode + PE_ode;      % Total Mechanical Energy
Lagrangian_ode = KE_ode - PE_ode;
Hamiltonian_ode = p1_ode .* omega_1_ode + p2_ode .* omega_2_ode - Lagrangian_ode;
%%
% --- 7. Plot the Energy Over Time ---
lgd = {'RK5','Yoshida 4', 'Velocity Verlet', 'Normal Verlet', 'RK4', 'Euler', 'Implicit Euler', 'TR-BDF2',  'ODE23s'};
figure;
T = tiledlayout(2, 2);
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

nexttile(2,[2, 1])
plot(KE, PE)
hold on
plot(KE_ode, PE_ode, 'k.');
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')
legend(lgd);

% ylim([-60, 60]);
% xlim([-5, 300]);

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
title('Angular Displacement of Bob 1')
legend(lgd);
ylim([-5, 5])

nexttile
plot(t_output_points, omega1)
hold on
plot(t_sol, X_sol(:, 3), 'k--');
grid on
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity of Bob 1')
legend(lgd);
ylim([-20, 20])

%%
X_RKF5 =    [theta1(:, 1), theta2(:, 1)];
X_Verlet =  [theta1(:, 2), theta2(:, 2)];
X_RK4 =     [theta1(:, 3), theta2(:, 3)];
X_Euler =   [theta1(:, 4), theta2(:, 4)];
X_TR    =   [theta1(:, 8), theta2(:, 8)];
X_ODE =     [theta_1_ode, theta_2_ode];

plot_pendulums(t_output_points, X_ODE, X_TR, L1, L2);
           
%%

function plot_pendulums(t_sol, X_sol, Y_sol, L1, L2)
    % --- Optional: Visualize the Pendulum Motion ---
    % This requires converting angles to Cartesian coordinates
    x1 = L1 * sin(X_sol(:, 1));
    y1 = -L1 * cos(X_sol(:, 1)); % Y-axis pointing down
    x2 = x1 + L2 * sin(X_sol(:, 2));
    y2 = y1 - L2 * cos(X_sol(:, 2));

    a1 = L1 * sin(Y_sol(:, 1));
    b1 = -L1 * cos(Y_sol(:, 1)); % Y-axis pointing down
    a2 = a1 + L2 * sin(Y_sol(:, 2));
    b2 = b1 - L2 * cos(Y_sol(:, 2));
    
    
    % You can also animate it:
    figure;
    h_link1 = plot([0 x1(1)], [0 y1(1)], 'k-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    h_link2 = plot([x1(1) x2(1)], [y1(1) y2(1)], 'k-o', 'LineWidth', 2, 'MarkerSize', 8);

    b_link1 = plot([0 a1(1)], [0 b1(1)], 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    b_link2 = plot([a1(1) a2(1)], [b1(1) b2(1)], 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    axis([-3 3 -3 3]); axis equal; grid on;
    title('Double Pendulum Animation');
    legend('ODE23s', '', 'Verlet', '')
    pause(0.5)
    for k = 1:10:length(t_sol) % Animate every 5th frame
        set(h_link1, 'XData', [0 x1(k)], 'YData', [0 y1(k)]);
        set(h_link2, 'XData', [x1(k) x2(k)], 'YData', [y1(k) y2(k)]);

        set(b_link1, 'XData', [0 a1(k)], 'YData', [0 b1(k)]);
        set(b_link2, 'XData', [a1(k) a2(k)], 'YData', [b1(k) b2(k)]);
        % drawnow;
        pause(0.01)
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
    den = ( 2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2) );

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
    
    X_next =  X_vec + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4 - (9/50) * k5 + (2/55) * k6;
end

function X_next = evaluateYoshida(func, dt, t, X_vec, args)
    % 4th Order Yoshida Integrator for systems with explicit acceleration function.
    % This function performs one integration step from t to t+dt.
    %
    % Inputs:
    %   func:       Function handle to compute acceleration.
    %               Expected signature: a = func(t, x, v, args_for_func{:})
    %               where x is position and v is velocity (potentially half-step velocity).
    %   dt:         Time step.
    %   t:          Current time.
    %   X_vec:      Current state vector [x_current; v_current; a_current].
    %               x_current: Position(s) at time t.
    %               v_current: Velocity(s) at time t.
    %               a_current: Acceleration(s) at time t (calculated from previous step).
    %   varargin:   Optional additional arguments to pass directly to func.
    %
    % Output:
    %   X_next:     Next state vector [x_new; v_new; a_new].
    %               x_new: Position(s) at time t+dt.
    %               v_new: Velocity(s) at time t+dt.
    %               a_new: Acceleration(s) at time t+dt.

    % --- 1. Define Correct Yoshida Coefficients ---
    % These are for the common A-B-A composition (position-velocity-position)
    % which requires 4 position updates and 3 velocity updates.
    p = 2^(1/3);
    w0 = -p / (2 - p); % approx -0.7657
    w1 = 1 / (2 - p);  % approx 1.3512

    c1 = w1 / 2; c4 = c1;
    c2 = (w0 + w1) / 2; c3 = c2;
    
    d1 = w1; d3 = d1;
    d2 = w0;

    % --- 2. Robust State Variable Extraction ---
    % Assuming X_vec = [x_components; v_components; a_old_components]
    N = numel(X_vec); % Use numel for total elements
    num_dims = N / 3; % Number of dimensions (e.g., 1 for 1D, 2 for 2D)

    if mod(N, 3) ~= 0 || N < 3
        error('Input state vector X_vec must have a total number of elements divisible by 3 (position, velocity, and acceleration components for each dimension). For 1D, N=3; for 2D, N=6, etc.');
    end
    
    x = X_vec(1:num_dims);            % Position(s) at time t
    v = X_vec(num_dims+1 : 2*num_dims); % Velocity(s) at time t
    a_old = X_vec(2*num_dims+1 : 3*num_dims); % Acceleration(s) at time t

    % --- 3. Integration Steps (Composition) ---
    % Here, 'A' step is x update, 'B' step is v update, 'func' provides 'a' (acceleration)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c1 * dt)
    x = x + c1 .* v .* dt;
    
    % B(d1 * dt) - calculate 'a' at new position (x) and current velocity (v)
    % The 'v' here (v_current) is often used for damping in the 'a' calculation
    a_temp1 = get_acceleration(func,t, [x; v], args); % Calculate a(t+c1*dt) approximately
    v = v + d1 .* a_temp1 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c2 * dt)
    x = x + c2 .* v .* dt;

    % B(d2 * dt)
    a_temp2 = get_acceleration(func,t, [x; v], args); % Calculate a(t+(c1+c2)*dt) approximately
    v = v + d2 .* a_temp2 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c3 * dt)
    x = x + c3 .* v .* dt;

    % B(d3 * dt)
    a_temp3 = get_acceleration(func,t, [x; v], args); % Calculate a(t+(c1+c2+c3)*dt) approximately
    v = v + d3 .* a_temp3 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c4 * dt) - Final position update
    x = x + c4 .* v .* dt;

    % --- 4. Recalculate Final Acceleration for Output Consistency ---
    % After all position and velocity updates, calculate the acceleration
    % corresponding to the final (t+dt) state.
    a_final = get_acceleration(func,t+dt, [x; v], args); % a(t+dt)

    % --- 5. Prepare Output State Vector ---
    % X_next contains [new_position; new_velocity; new_acceleration]
    X_next = [x; v; a_final];
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

function X_next = integrate_VelocityVerlet(func, dt, t, X_vec, args)
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
    if mod(N, 3) == 0
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
    X_vec_for_func = [x_new; v_half_step];
    a_new = get_acceleration(func, t, X_vec_for_func, args);% a(t + dt

    % 4. Calculate full-step velocity
    v_new = v_half_step + a_new * (dt / 2); % v(t + dt)
    
    % X_next contains [new_position; new_velocity; new_acceleration]
    X_next = [x_new; v_new; a_new]; 
end

function X_next = integrate_verlet(func, dt, t, X_vec, args)
    % Position Verlet integration method for dynamics.
    % X_vec should be [current_position; previous_position]
   
    % Extract current state variables
    [N, ~] = size(X_vec);
    if mod(N, 3) == 0
        if N == 6
            % Extract the acceleration component (dvx is the second element of dX)
            x_current = X_vec(1:2); 
            x_current = x_current(:);
            v_estimated = X_vec(3:4); 
            v_estimated = v_estimated(:);
            x_previous = X_vec(5:6); 
            x_previous = x_previous(:);
            
        else
            x_current   = X_vec(1); % x(t)
            v_estimated = X_vec(2); % v(t)
            x_previous  = X_vec(3); % x(t-dt) - displacement from the previous step

        end
    else
        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end

    % 1. Calculate acceleration at current position
    X_vec_for_func = [x_current;v_estimated];
    a_current = get_acceleration(func, t, X_vec_for_func, args); % a(t)

    % 2. Update position
    x_new = 2 * x_current - x_previous + a_current * dt.^2; % x(t + dt)

    % 3. Calculate velocity for output (not used for propagation in core Verlet)
    % A better estimate for v(t+dt) for analysis is often (x(t+dt) - x(t-dt)) / (2*dt)
    v_new_estimated = (x_new - x_previous) / (2 * dt); % v(t + dt)

    % X_next contains [new_position;estimated_new_velocity; current_position (for next step's "previous"); ]
    X_next = [x_new; v_new_estimated; x_current]; 
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

function X_next = evaluateImplicitEuler(func, dt, t, X_vec, args)
    % evaluateImplicitEuler_Newton: Implements the Implicit Euler method using MATLAB's fsolve
    % (which uses a trust-region or quasi-Newton algorithm).
    % This approach is more robust for stiff ODEs compared to fixed-point iteration.
    %
    % Inputs:
    %   func:       Function handle for the ODE. Signature: dXdt = func(time, state_vector, args{:})
    %   dt:         Time step size.
    %   t:          Current time (t_n).
    %   X_vec:      Current state vector (X_n).
    %   varargin:   Optional additional arguments to pass to func.
    %
    % Output:
    %   X_next:     Approximation of the state vector at time t + dt (X_{n+1}).

    % Ensure 'optimoptions' and 'fsolve' are available in the MATLAB environment.
    % This code requires the Optimization Toolbox.

    t_next = t + dt;

    % --- Define the function to find the root of ---
    % We want to solve for X_next such that:
    % X_next - X_vec - dt * func(t_next, X_next, varargin{:}) = 0
    % Define this as an anonymous function for fsolve.
    % The 'args_cell' captures 'varargin' to pass to 'func' correctly.
    args_cell = args; 
    
    equation_to_solve = @(X_guess) X_guess - X_vec - dt * func(t_next, X_guess, args_cell{:});

    % --- Initial Guess for X_next ---
    % A good initial guess is crucial for fsolve convergence.
    % The Explicit Euler step is often used as a starting point.
    initial_guess = X_vec + dt * func(t, X_vec, args_cell{:});

    % --- Set fsolve options (optional, but recommended for control) ---
    % You can adjust these based on your specific problem's needs.
    % 'Display': 'off', 'final', 'iter'
    % 'FunctionTolerance': Tolerance on the function value
    % 'StepTolerance': Tolerance on the change in X_next
    options = optimoptions('fsolve', ...
                           'Display', 'off', ... % Suppress verbose output from fsolve
                           'FunctionTolerance', 1e-8, ... % Tolerance for F(X_next) close to zero
                           'StepTolerance', 1e-8);     % Tolerance for change in X_next

    % --- Call fsolve to solve for X_next ---
    % fsolve returns X_next, and potentially fval (value of the function at solution),
    % exitflag, output structure, and Jacobian. We only need X_next for this function.
    [X_next, ~, exitflag, output] = fsolve(equation_to_solve, initial_guess, options);

    % --- Check fsolve exit flag (optional, but good practice) ---
    if exitflag <= 0 % exitflag < 1 typically means no convergence
        warning('evaluateImplicitEuler_Newton:fsolveNoConvergence', ...
                'fsolve did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
                t, dt, exitflag, output.message);
    end
end

function [y_next] = trbdf2_step(func, dt, t, X_vec, args)
    % TRBDF2_STEP performs one step of the TR-BDF2 integration scheme.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   J: Function handle for the Jacobian of f with respect to y: J(t, y) = df/dy
    %   t: Current time
    %   y: Current solution vector
    %   h: Timestep size
    %   tol: Tolerance for the Newton's method
    %   max_iter: Maximum iterations for Newton's method
    %
    %   y_next: Solution at time t+h
    %   t_next: Time at the end of the step (t+h)
    
    % The optimal parameter gamma for L-stability
    gamma = 2 - sqrt(2);
    
    % --- Stage 1: Trapezoidal Rule (Implicit Step) ---
    % Solve for y_intermediate at t_gamma = t + gamma*h
    % The equation to solve is: y_intermediate - y - (gamma*h/2)*(f(t,y) + f(t_gamma, y_intermediate)) = 0
    y = X_vec;
    y_intermediate = X_vec; % Initial guess for Newton's method
    h = dt;
    % Calculate intermediate timestep
    t_intermediate = t + gamma * h;
    
    F1 = @(y_i_gamma) y_i_gamma - ( y + (gamma * h / 2) * ( func(t, y, args{:}) + func(t_intermediate, y_i_gamma, args{:}) ) );

    % --- Set fsolve options (optional, but recommended for control) ---
    % You can adjust these based on your specific problem's needs.
    % 'Display': 'off', 'final', 'iter'
    % 'FunctionTolerance': Tolerance on the function value
    % 'StepTolerance': Tolerance on the change in X_next
    options = optimoptions('fsolve', ...
                           'Display', 'off', ... % Suppress verbose output from fsolve
                           'FunctionTolerance', 1e-8, ... % Tolerance for F(X_next) close to zero
                           'StepTolerance', 1e-8);     % Tolerance for change in X_next

    % --- Call fsolve to solve for X_next ---
    % fsolve returns X_next, and potentially fval (value of the function at solution),
    % exitflag, output structure, and Jacobian. We only need X_next for this function.
    [y_intermediate, ~, exitflag, output] = fsolve(F1, y_intermediate, options);

    % --- Check fsolve exit flag (optional, but good practice) ---
    if exitflag <= 0 % exitflag < 1 typically means no convergence
        warning('evaluateImplicitEuler_Newton:fsolveNoConvergence', ...
                'fsolve did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
                t, dt, exitflag, output.message);
    end
    
    % --- Stage 2: BDF2 (Implicit Step) ---
    % Solve for y_next at t_next = t + h
    % The equation to solve is: y_next - (1/(2-gamma))*( ((1-gamma)^2/gamma)*y + (1/gamma)*y_intermediate) - (1-gamma)/(2-gamma)*h*f(t_next, y_next) = 0
    
    y_next = y_intermediate; % Initial guess for Newton's method
    F2 = @(y_ip1) y_ip1 - (1/(2-gamma)) * ( (1/gamma)*y_intermediate - ((1-gamma)^2/gamma)*y + (1-gamma)*h*func(t + h, y_ip1, args{:}) );
   
     % --- Call fsolve to solve for X_next ---
    % fsolve returns X_next, and potentially fval (value of the function at solution),
    % exitflag, output structure, and Jacobian. We only need X_next for this function.
    [y_next, ~, exitflag, output] = fsolve(F2, y_next, options);

    % --- Check fsolve exit flag (optional, but good practice) ---
    if exitflag <= 0 % exitflag < 1 typically means no convergence
        warning('evaluateImplicitEuler_Newton:fsolveNoConvergence', ...
                'fsolve did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
                t, dt, exitflag, output.message);
    end

end


function [y_next, h] = adaptive_ODE12t(func, dt, t, X_vec, args)
    % adaptive_ODE12t performs one step of numerical integration.
    % An adaptive timestep with euler as first step and trapezoidal as
    % finer resolution step.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   dt: Current time step
    %   t: Current time
    %   X_vec: Current solution vector
    %   args: Arguments for function func(t, x, args)
    %
    %   y_next: Solution at time t+h
    %   h: Time step at the end of the current step 
    y = X_vec;
    h = dt;
    hasToRunAgain = true;
    rtol = 1e-6;
    atol = 1e-6;
    h_min = 100 * eps;
    h_max = 1;
    max_iters = 5;

    % Calculate tolerance value
    tolerance = max(rtol * norm(y), atol);

    % Start counter
    iters = 0;
    while hasToRunAgain
        % --- Stage 1: RK1 Rule (Euler Step) ---
        % Solve for y_n+1 at t_n+1 = t + h
        k1 = h * func(t, y, args{:});
        
        y_euler = X_vec + k1;
    
        % --- Stage 2: Trapezoidal rule (Heun's Method) ---
        % Solve for y_next at t_next = t + h
        k2 = h * func(t + h, y_euler, args{:});
        
        y_heun = y + 0.5 * (k1 + k2);
        % Estimate error betweem two schemes
        error = norm(y_heun - y_euler);
 
        % Calculate scaling factor
        S = 0.9 * (tolerance / error)^(1/2);
        
        if error <= tolerance
            % Accept current timestep
            y_next = y_heun;
            % Increase the step size
            h = min(h_max, S * h);
            hasToRunAgain = false;
        else
            % Reject the step
            h = max(h_min, S * h);
        end
        % Increment counter
        iters = iters + 1;
        % break out of loop to avoid infinite while loop
        % Can break out of loop if h == h_min because no refinement will
        % happen
        if iters >= max_iters || h == h_min
            hasToRunAgain = false;
            % Assume trapezoidal rule to be more accurate
            y_next = y_heun;
        end
    end

end

