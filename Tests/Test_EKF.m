set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear;close all;clc
% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 10;%0.37;     % Spring stiffness (N/m)
c = 0;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
F = 0.1 * sin(2 * pi * t_output_points) .* (t_output_points <= 5) + randn(size(t_output_points)) * 0.01;
F_ext = @(t) interp1(t_output_points, F, t, "linear", "extrap");

x1 = zeros(length(t_output_points), 3);
v1 = zeros(size(x1));
a1 = zeros(size(x1));
e1 = zeros(size(x1));
normP = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), 13);

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.1;    % Initial displacement of mass
v1(1, :) = 0.0;    % Initial velocity of mass

initial_state = [x1(1, 1); v1(1, 1)];
args = {F_ext, m, k, c};

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% odeXX takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","on"); 
my_dynamics = @(t, X) springMassDamperDynamics(t, X, args{:});
tic;
[t_sol, X_sol] = ode23tb(my_dynamics,t_output_points , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)

options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
numStates = 2; % Displacement and velocity
numStatesMeasure = 1; % Displacement
EKF_opts = struct('ThreshScal', sqrt(eps) * ones(numStates, 1),...
                  'facF', [], ...
                  'facH', [], ...
                  'ProcessNoiseCov', eye(numStates) * 1e-6, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.01, ...
                  'ErrorCov', zeros(numStates), ...
                  'dt', 1/fs_output, ...
                  'ControlVec', 1);
args_ekf = {m, k, c};
my_dynamics_ekf = @(t, X) springMassDamperDynamics(t, X, F_ext,args_ekf{:});
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);
    %%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1);];
    % Y_vec = X_sol(i, 1).';
    tic;
    [X_next, ~, P, EKF_opts] = extendedKalmanFilter(my_dynamics_ekf, t, X_vec, [], EKF_opts);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    e1(i+1, 1) = rms(X_sol(i, :).' - X_next);
    normP(i+1) = norm(P);

    %%%%%%%%%%%%%% TR-BDF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2);];
    tic;
    X_next = trbdf2_step(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    e1(i+1, 2) = rms(X_sol(i, :).' - X_next);

    %%%%%%%%% ODE23tb at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 3); v1(i, 3);];
    tic;
    [~, X_next] = ode23tb(@(t, X) springMassDamperDynamics(t, X, args{:}),[t, t_target] , X_vec, options);
    time_to_solve(i, 3) = toc;
    x1(i+1, 3) = X_next(end, 1);
    v1(i+1, 3) = X_next(end, 2);
    e1(i+1, 3) = rms( X_sol(i, :) - X_next(end, :) );
end

lgd = {'EKF', 'TR-BDF2', 'ode23tb'};

figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Displacement [m]')

nexttile
plot(t_output_points, v1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Velocity [m/s]')

nexttile
semilogy(t_output_points, e1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('RMS Error [-]')

nexttile
semilogy(t_output_points, normP)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Norm of Error Covariance [-]')


nexttile
semilogy(t_output_points, time_to_solve)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('log_{10}(CPU Time) [s]')

%%
clear;close all;clc
% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 10;%0.37;     % Spring stiffness (N/m)
c = 0;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 20;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
F = 0.1 * sin(2 * pi * t_output_points) .* (t_output_points <= t_final * 0.67);% + randn(size(t_output_points)) * 0.01;
F_ext = @(t) interp1(t_output_points, F, t, "linear", "extrap");

x1 = zeros(length(t_output_points), 2);
v1 = zeros(size(x1));
M = zeros(size(x1));
K = zeros(size(x1));
C = zeros(size(x1));
e1 = zeros(size(x1));
normP = zeros(size(x1));

time_to_solve = zeros(size(x1));

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.1;    % Initial displacement of mass
v1(1, :) = 0.0;    % Initial velocity of mass
M(1, :) = m;    % Initial displacement of mass
K(1, :) = k;    % Initial velocity of mass
C(1, :) = c;    % Initial displacement of mass

initial_state = [x1(1, 1); v1(1, 1)];
args = {F_ext, m, k, c};

% --- 4. Solve the Ordinary Differential Equation (ODE) ---
% odeXX takes:
%   - @(t, X) brush_dynamics(t, X, m, k, c): An anonymous function that calls
%     your dynamics function, passing system parameters.
%   - t_output_points: The specific time points at which you want the solution.
%   - initial_state: The initial values of your state variables.
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","on"); 
my_dynamics = @(t, X) springMassDamperDynamics(t, X, args{:});
tic;
[t_sol, X_sol] = ode23tb(my_dynamics,t_output_points , initial_state, options);
X_sol = X_sol + 0.05 * randn(size(X_sol));
fprintf('Elapsed time: %gs \n', toc)

options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
numStates = 5; % Displacement and velocity
numStatesMeasure = 1; % Displacement
EKF_opts = struct('ThreshScal', sqrt(eps) * ones(numStates, 1),...
                  'facF', [], ...
                  'facH', [], ...
                  'ProcessNoiseCov', eye(numStates) * 0.01, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', zeros(numStates), ...
                  'dt', 1/fs_output, ...
                  'ControlVec', []);

args_ekf = {F_ext, m, k, c};
my_dynamics_ekf = @(t, X) springMassDamperParamEst(t, X, dt, @springMassDamperDynamics, args_ekf);
my_dynamics_ekf_ode = @(t, X) springMassDamperParamEst_ODE(t, X, dt, @springMassDamperDynamics, args_ekf);
%%
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);
    %%%%%%%%%%%%%% TR-BDF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1); M(i, 1); K(i, 1); C(i, 1)];
    Y_vec = X_sol(i, :).';
    tic;
    [X_next, ~, P, EKF_opts] = extendedKalmanFilter(my_dynamics_ekf, t, X_vec, Y_vec, EKF_opts);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    M(i+1, 1) = X_next(3);
    K(i+1, 1) = X_next(4);
    C(i+1, 1) = X_next(5);
    
    e1(i+1, 1) = rms(X_sol(i, :).' - X_next(1:2));
    normP(i+1, 1) = norm(P);

    %%%%%%%%%%%%%% ODE23tb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2); M(i, 2); K(i, 2); C(i, 2)];
    Y_vec = X_sol(i, :).';
    tic;
    [X_next, ~, P, EKF_opts] = extendedKalmanFilter(my_dynamics_ekf_ode, t, X_vec, Y_vec, EKF_opts);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    M(i+1, 2) = X_next(3);
    K(i+1, 2) = X_next(4);
    C(i+1, 2) = X_next(5);

    e1(i+1, 2) = rms(X_sol(i, :).' - X_next(1:2));
    normP(i+1, 2) = norm(P);
end
%%
lgd = {'Measured', 'TR-BDF2', 'ode23tb'};

figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
hold on
plot(t_sol, X_sol(:, 1))
plot(t_output_points, x1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Displacement [m]')

nexttile
hold on
plot(t_sol, X_sol(:, 2))
plot(t_output_points, v1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Velocity [m/s]')

lgd = {'TR-BDF2', 'ode23tb'};

nexttile
semilogy(t_output_points, e1)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('RMS Error [-]')

nexttile
semilogy(t_output_points, time_to_solve)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('log_{10}(CPU Time) [s]')


figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
semilogy(t_output_points, normP)
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Norm of Error Covariance [-]')

nexttile
plot(t_output_points, (M))
grid on
title('Mass Parameter')

nexttile
plot(t_output_points, (K))
grid on
title('Stiffness Parameter')

nexttile
plot(t_output_points, (C))
grid on
title('Damping Parameter')


%%
function dX = springMassDamperDynamics(t, X, Forcing, m, k, c)

    dx = X(2, :);
    dvx = (Forcing(t) - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end

function states = springMassDamperParamEst(t, X_vec, dt, func, args)

    % % options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
    % % t_target = t + dt;

    % Extract states to integrate
    X_states = X_vec(1:2, :);
    % Update parameters of system
    args{2} = X_vec(3);
    args{3} = X_vec(4);
    args{4} = X_vec(5);
    % Apply integration method
    X_next = trbdf2_step(func, dt, t, X_states, args);

    % Augment output states
    m = X_vec(3);
    k = X_vec(4);
    c = X_vec(5);
    states = [X_next(:);
              m;
              k;
              c];
end

function states = springMassDamperParamEst_ODE(t, X_vec, dt, func, args)

    options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
    t_target = t + dt;

    % Extract states to integrate
    X_states = X_vec(1:2, :);
    % Update parameters of system
    args{2} = X_vec(3);
    args{3} = X_vec(4);
    args{4} = X_vec(5);
    % Apply integration method
    [~, X] = ode23tb(@(t, X) func(t, X, args{:}),[t, t_target] , X_states, options);
    X_next = X(end, :);

    % Augment output states
    m = X_vec(3);
    k = X_vec(4);
    c = X_vec(5);
    states = [X_next(:);
              m;
              k;
              c];
end