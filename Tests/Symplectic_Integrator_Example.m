clearvars; close all; clc;
set(0, "DefaultFigureWindowStyle", "docked")
%%
clear;close all;clc;

% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 500;%0.37;     % Spring stiffness (N/m)
c = 5;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7;           % End time (s)
fs_output = 2000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
% F = 0.001 * sin(2 * pi * t_output_points) .* (t_output_points <= 5);% + randn(size(t_output_points)) * 0.0001;
% Sawtooth sweep parameters
F_amplitude = 1e3;    % Maximum force amplitude
f_start = 0;      % Starting frequency in Hz
f_end = 10.0;        % Ending frequency in Hz

F = F_amplitude * sawtooth(2 * pi * (f_start + (f_end - f_start) * (t_output_points / t_final) / 2) .* t_output_points);

% F_ext = @forcing_func;
F_ext = @(t, x) interp1(t_output_points, F, t, "linear");

x1 = zeros(length(t_output_points), 13);
v1 = zeros(length(t_output_points), 13);
a1 = zeros(length(t_output_points), 7);
e1 = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), 13);

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.0;    % Initial displacement of mass
v1(1, :) = 0.0;    % Initial velocity of mass

initial_state = [x1(1, 1); v1(1, 1)];
args = {m, k, c, F_ext};

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% odeXX takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
my_dynamics = @(t, X) springMassDamperDynamics(t, X, args{:});
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","on"); 
tic;
[t_sol, X_sol] = ode23s(my_dynamics,t_output_points , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)
%%
% X_sol = zeros(length(t_output_points), 2);
% t_sol = t_output_points;
% For first iteration of adaptive RK(4)5
tRK_current = t_init;
tTR_current = t_init;
hRK_current = dt;
hTR_current = dt;


options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off"); 
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);

    %%%%%%%%%%%%% Adaptive Heun's Method %%%%%%%%%%%%%
    X_vec = [x1(i, 11); v1(i, 11)];
    while tTR_current < t_target
        % Calculate required step
        hTR_current = min(hTR_current,  t_target - tTR_current);

        % Call the adaptive step function
        [X_next, hTR_next] = adaptive_ODE12t(my_dynamics, hTR_current, tTR_current, X_vec);

        % Update current time based on the step taken
        tTR_current = tTR_current + hTR_current;
        % Update time step
        hTR_current = hTR_next;
        % Update solution
        X_vec = X_next;
    end
    time_to_solve(i, 11) = toc;
    x1(i+1, 11) = X_next(1);
    v1(i+1, 11) = X_next(2);
    e1(i+1, 11) = rms(X_sol(i, :).' - X_next);

    %%%%%%%%%%%%%%%%%% Adaptive RKF45 - III method %%%%%%%%%%%%%%%%%%%%%%%%
    % Advance the solution until we reach the target time
    X_vec = [x1(i, 9); v1(i, 9)];
    tic;
    while tRK_current < t_target
        % Calculate required step
        hRK_current = min(hRK_current,  t_target - tRK_current);

        % Call the adaptive step function
        [X_next, hRK_next] = adaptiveRK45(my_dynamics, hRK_current, tRK_current, X_vec);

        % Update current time based on the step taken
        tRK_current = tRK_current + hRK_current;
        % Update time step
        hRK_current = hRK_next;
        % Update solution
        X_vec = X_next;
    end
    time_to_solve(i, 9) = toc;
    x1(i+1, 9) = X_next(1);
    v1(i+1, 9) = X_next(2);
    e1(i+1, 9) = rms(X_sol(i, :).' - X_next);

    % %%%%%%%%%%%%%%%%%% RKF5 Method %%%%%%%%%%%%%%%%%%%%%%%%
    % X_vec = [x1(i, 1); v1(i, 1)];
    % tic;
    % X_next = evaluateRKF5(@springMassDamperDynamics, dt, t, X_vec, args);
    % time_to_solve(i, 1) = toc;
    % x1(i+1, 1) = X_next(1);
    % v1(i+1, 1) = X_next(2);
    % e1(i+1, 1) = norm(X_sol(i, :).' - X_next);

    % % %%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%%%
    % % X_vec = [x1(i, 2); v1(i, 2); a1(i, 2);];
    % % tic;
    % % X_next = evaluateYoshida(@springMassDamperDynamics, dt, t, X_vec, args);
    % % time_to_solve(i, 2) = toc;
    % % x1(i+1, 2) = X_next(1);
    % % v1(i+1, 2) = X_next(2);
    % % a1(i+1, 2) = X_next(3);
    % % e1(i+1, 2) = norm(X_sol(i, :) - X_next(1:2));
    % % 
    % % % %%%%%%%%%%%%%%%%%%%% Velocity Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    % % X_vec = [x1(i, 3); v1(i, 3); a1(i, 3);];
    % % tic;
    % % X_next = integrate_VelocityVerlet(@springMassDamperDynamics, dt, t, X_vec, args);
    % % time_to_solve(i, 3) = toc;
    % % x1(i+1, 3) = X_next(1);
    % % v1(i+1, 3) = X_next(2);
    % % a1(i+1, 3) = X_next(3);
    % % e1(i+1, 3) = norm(X_sol(i, :) - X_next(1:2));
    % % 
    % % %%%%%%%%%%%%%%%%%%%% Normal Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%
    % % if i~= 1
    % %     X_vec = [x1(i, 4); v1(i, 4); a1(i, 4);];
    % %     tic;
    % %     X_next = integrate_verlet(@springMassDamperDynamics, dt, t, X_vec, args);
    % %     time_to_solve(i, 4) = toc;
    % %     x1(i+1, 4) = X_next(1);
    % %     v1(i+1, 4) = X_next(2);
    % %     a1(i+1, 4) = X_next(3);
    % %     e1(i+1, 4) = norm(X_sol(i, :) - X_next(1:2));
    % % 
    % % else
    % %     X_vec = [x1(1, 4); v1(1, 4);];
    % %     X_next = evaluateEuler(@springMassDamperDynamics, dt, t_output_points(1), X_vec, args);
    % %     a1(i+1, 4) = x1(1, 4);
    % %     x1(i+1, 4) = X_next(1);
    % %     v1(i+1, 4) = X_next(2);
    % %     e1(i+1, 4) = norm(X_sol(i, :).' - X_next);
    % % 
    % % end
    % % 
    % % %%%%%%%%%%%%%%%%%%%% RK 4 Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % X_vec = [x1(i, 5); v1(i, 5);];
    % % tic;
    % % X_next = evaluateRK4(@springMassDamperDynamics, dt, t, X_vec, args);
    % % time_to_solve(i, 5) = toc;
    % % x1(i+1, 5) = X_next(1);
    % % v1(i+1, 5) = X_next(2);
    % % e1(i+1, 5) = norm(X_sol(i, :).' - X_next);
    % % 
    % % %%%%%%%%%%%%%%%%%%%% Eulers Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % X_vec = [x1(i, 6); v1(i, 6);];
    % % tic;
    % % X_next = evaluateEuler(@springMassDamperDynamics, dt, t, X_vec, args);
    % % time_to_solve(i, 6) = toc;
    % % x1(i+1, 6) = X_next(1);
    % % v1(i+1, 6) = X_next(2);
    % % e1(i+1, 6) = norm(X_sol(i, :).' - X_next);

    %%%%%%%%%%%%%%% Implicit Euler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 7); v1(i, 7);];
    tic;
    X_next = evaluateImplicitEuler(my_dynamics, dt, t, X_vec);
    time_to_solve(i, 7) = toc;
    x1(i+1, 7) = X_next(1);
    v1(i+1, 7) = X_next(2);
    e1(i+1, 7) = rms(X_sol(i, :).' - X_next);


    %%%%%%%%%%%%%% TR-BDF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 8); v1(i, 8);];
    tic;
    X_next = evaluateTRBDF2(my_dynamics, dt, t, X_vec);
    time_to_solve(i, 8) = toc;
    x1(i+1, 8) = X_next(1);
    v1(i+1, 8) = X_next(2);
    e1(i+1, 8) = rms(X_sol(i, :).' - X_next);

    %%%%%%%%% ODE23tb at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 10); v1(i, 10);];
    tic;
    [~, X_next] = ode23tb(my_dynamics,[t, t_target] , X_vec, options);
    time_to_solve(i, 10) = toc;
    x1(i+1, 10) = X_next(end, 1);
    v1(i+1, 10) = X_next(end, 2);
    e1(i+1, 10) = rms( X_sol(i, :) - X_next(end, :) );

    %%%%%%%%% ODE23t at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 12); v1(i, 12);];
    tic;
    [~, X_next] = ode23t(my_dynamics,[t, t_target] , X_vec, options);
    time_to_solve(i, 12) = toc;
    x1(i+1, 12) = X_next(end, 1);
    v1(i+1, 12) = X_next(end, 2);
    e1(i+1, 12) = rms( X_sol(i, :) - X_next(end, :) );

    %%%%%%%%% ODE23s at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 13); v1(i, 13);];
    tic;
    [~, X_next] = ode23s(my_dynamics,[t, t_target] , X_vec, options);
    time_to_solve(i, 13) = toc;
    x1(i+1, 13) = X_next(end, 1);
    v1(i+1, 13) = X_next(end, 2);
    e1(i+1, 13) = rms( X_sol(i, :) - X_next(end, :) );

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

%%
% --- 7. Plot the Energy Over Time ---
% lgd = {'RK5','Yoshida 4', 'Velocity Verlet', 'Normal Verlet', 'RK4', 'Euler', 'Implicit Euler', 'TR-BDF2', 'adaptiveRK45', 'ODE23t at each time step', 'AdaptiveODE12t' ,  'ODE23t'};
% lgd = {'Implicit Euler', 'TR-BDF2', 'adaptiveRK45', 'ODE23t at each time step', 'adaptive12t' ,  'ODE23t'};
lgd = {'Implicit Euler','TR-BDF2', 'adaptiveRK45','ODE23tb at each time step','Adaptive Heun','ODE23t at each time step','ODE23s at each time step', 'ODE23s'};


figure;
T = tiledlayout(2, 2);
T.Padding ="compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, TotalMechanicalEnergy(:, end-5:end));
hold on
% plot(t_output_points, TotalMechanicalEnergy(:,1))
plot(t_sol, TotalMechanicalEnergy_ode, 'k--');
xlabel('Time (s)');
ylabel('Total Mechanical Energy (J)');
title('Energy of SDOF System Over Time');
grid on;
legend(lgd);
% ylim([-0.02, 0.02]);

nexttile(2,[2, 1])
plot(KE(:, end-5:end), PE(:, end-5:end), '-.')
hold on
% plot(KE(:,1), PE(:,1))
plot((KE_ode), (PE_ode), 'k-.');
grid on
xlabel('Kinetic Energy [J]')
ylabel('Potential Energy [J]')
legend(lgd);
% axis([0, 0.02, 0, 0.02]);

nexttile
semilogy(t_output_points, time_to_solve(:, end-5:end))
hold on
% semilogy(t_output_points, time_to_solve(:,1))
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
plot(t_output_points, x1(:, end-5:end))
hold on
% plot(t_output_points, x1(:,1))
plot(t_sol, X_sol(:, 1), 'k--');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
legend(lgd);
ylim([-1, 1]);

nexttile
plot(t_output_points, v1(:, end-5:end))
hold on
% plot(t_output_points, v1(:,1))
plot(t_sol, X_sol(:, 2), 'k--');
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend(lgd);
ylim([-2, 2]);

nexttile
semilogy(t_output_points, e1(:, end-5:end))
hold on
% semilogy(t_output_points, e1(:, 1))
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
    X_next = evaluateTRBDF2(@doublePendulumDynamics, dt, t, X_vec, args);
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

% % function F = F_ext(t)
% %     if t < 10
% %         F = 0.1 * sin(2 * pi * t);
% %     else
% %         F = 0;
% %     end
% %     % F = 0;
% % end

function dX = springMassDamperDynamics(t, X, m, k, c, Forcing)

    dx = X(2, :);
    dvx = (Forcing(t, X) - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end

function F = forcing_func(t, X)
    % This is your state-dependent forcing function.
    % It calculates the forces based on time and the current system state.
    v = X(2, :);

    % Example: Forcing that depends on velocity
    F_magnitude = 10 * X(1, :); % Example magnitude
    F = F_magnitude * sign(v); % Example: a force that opposes velocity
end



















