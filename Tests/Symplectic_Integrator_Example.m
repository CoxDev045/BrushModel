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

function X_next = evaluateRKF5(dt, t, X_vec, m, k, c)
    k1 = dt * springMassDamperDynamics(t, X_vec, ...
                                       m, k, c);
    k2 = dt * springMassDamperDynamics(t, X_vec + 0.25 *  k1, ...
                                       m, k, c);
    k3 = dt * springMassDamperDynamics(t, X_vec + 3/32 *  k1 + 9/32 *  k2, ...
                                       m, k, c);
    k4 = dt * springMassDamperDynamics(t, X_vec + 1932/2197 *  k1 - 7200/2197 *  k2 + 7296/2197 *  k3, ...
                                       m, k, c);
    k5 = dt * springMassDamperDynamics(t, X_vec + 439/216 *  k1 - 8 *  k2 + 3680/513 *  k3 - 845/4104 *  k4, ...
                                       m, k, c);
    k6 = dt * springMassDamperDynamics(t, X_vec - 8/27 * k1 + 2 *  k2 - 3544/2565 *  k3 + 1859/4104 *  k4 - 11/40 *  k5, ...
                                       m, k, c);
    
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

function X_next = integrate_VelocityVerlet(dt, t, X_vec, m, k, c)
    % Velocity Verlet integration method for spring-mass-damper system.
    % X_vec should be [current_position; current_velocity; previous_acceleration]
    
    persistent a_old firstRun
    if isempty(firstRun)
        a_old = 0;
        firstRun = 1;
    end


    x_current = X_vec(1);     % x(t)
    v_current = X_vec(2);     % v(t)
    % a_old = X_vec(3);         % a(t) - acceleration from the previous step

    % 1. Calculate velocity at half-step
    v_half_step = v_current + a_old * (dt / 2);

    % 2. Update position
    x_new = x_current + v_half_step * dt;

    % 3. Calculate new acceleration based on new position
    % (Assuming F_ext(t) is defined and accessible, e.g., as a function handle
    % or global. For general use, it's better to pass it as an argument if it's dynamic.)
    % Using v_half_step for damping force is common here.
    a_new = (F_ext(t) - k * x_new - c * v_half_step) / m; % a(t + dt)

    % 4. Calculate full-step velocity
    v_new = v_half_step + a_new * (dt / 2); % v(t + dt)
    
    % X_next contains [new_position; new_velocity; new_acceleration]
    a_old = a_new;
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

function X_next = evaluateRK4(dt, t, X, m, k ,c)

    k1 = dt * springMassDamperDynamics(t, X, ...
                                       m, k, c);
    k2 = dt * springMassDamperDynamics(t, X + 0.5 * k1, ...
                                       m, k, c);
    k3 = dt * springMassDamperDynamics(t, X + 0.5 * k2, ...
                                       m, k, c);
    k4 = dt * springMassDamperDynamics(t, X + k3, ...
                                       m, k, c);
    
    X_next = X + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function X_next = evaluateEuler(dt, t, X_vec, m, k, c)
    
    x = X_vec(1);
    v = X_vec(2);

    a = (F_ext(t) - k * x - c * v) / m;
    %%%% Calculate velocity%%%%%%%%%%%%
    v = v + a * dt;
    %%%% Calculate displacement %%%%%%%%%%%%
    x = x + v * dt;

    X_next = [x; v];
end