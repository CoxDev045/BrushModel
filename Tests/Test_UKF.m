set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear;close all;clc
% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 500;%0.37;     % Spring stiffness (N/m)
c = 5;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7.5;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
% F = 1 * sin(0.5 * pi * t_output_points) .* (t_output_points <= 5) + randn(size(t_output_points)) * 0.1;
% % % Pulse parameters
% % pulse_height = 1000; % Amplitude of each pulse
% % pulse_width = 0.05; % Duration of each pulse
% % pulse_separation = 2.5; % Time between the start of pulses
% % 
% % % Generate the pulse train
% % F = zeros(size(t_output_points));
% % for t = 0:pulse_separation:t_final - pulse_separation
% %     % Find the indices for the current pulse
% %     pulse_indices = (t_output_points >= t) & (t_output_points < t + pulse_width);
% % 
% %     % Apply the pulse height to the corresponding indices
% %     F(pulse_indices) = pulse_height;
% % end

% Sinusoidal sweep parameters
F_amplitude = 100;    % Maximum force amplitude
f_start = 0;      % Starting frequency in Hz
f_end = 10.0;        % Ending frequency in Hz

% Generate the sinusoidal sweep
F = F_amplitude * sawtooth(2 * pi * (f_start + (f_end - f_start) * (t_output_points / t_final) / 2) .* t_output_points);

F_ext = @(t) interp1(t_output_points, F, t, "linear");


x1 = zeros(length(t_output_points), 3);
v1 = zeros(size(x1));
a1 = zeros(size(x1));
e1 = zeros(size(x1));
normP = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), size(x1, 2));

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.0;    % Initial displacement of mass
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
my_dynamics = @(t, X) springMassDamperDynamicsExpanded(t, X, args{:});
tic;
[~, X_sol] = ode23tb(my_dynamics,t_output_points , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)

X_sol = X_sol + 0.1 * randn(size(X_sol));

options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
numStates = 2; % Displacement and velocity
numStatesMeasure = 1; % Displacement and velocity
args_ukf = {F_ext, m, k, c};

% Initialise error with error at initial time
X_vec = [x1(1, 1); v1(1, 1)];
X_pred_init = springMassDamper_step(t_init, X_vec, dt, @springMassDamperDynamics, args_ukf);
X_sol_init = X_sol(1, :).';
P_init = (X_sol_init - X_pred_init) * (X_sol_init - X_pred_init).' + eye(numStates) * dt; % add identity matrix to stabilise P_init

UKF_opts = struct('ProcessNoiseCov', eye(numStates) * dt, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', P_init, ...
                  'dt', 1/fs_output, ...
                  'NumParticles', numStates);

EKF_opts = struct('ThreshScal', sqrt(eps) * ones(numStates, 1),...
                  'facF', [], ...
                  'facH', [], ...
                  'ProcessNoiseCov', ones(numStates) * dt, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', P_init, ...
                  'dt', 1/fs_output, ...
                  'NumParticles', numStates);
A = [0,1;-k/m, -c/m];
B = [0;1/m];
Phi = expm(A * dt);
Gamma = A \ ( ( Phi - eye(size(Phi)) ) * B) ;

KF_opts = struct('StateTransition', Phi, ...
                  'ControlMatrix', Gamma,...
                  'Forcing', F_ext,...
                 'MeasurementModel', eye(numStatesMeasure, numStates), ... 
                 'ProcessNoiseCov', ones(numStates) * dt, ...
                 'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                 'ErrorCov', P_init, ...
                 'dt', 1/fs_output, ...
                 'NumParticles', numStates);

P_history = zeros(numStates, numStates, length(t_output_points), size(x1, 2));
P_history(:, :, 1, 1) = P_init;
P_history(:, :, 1, 2) = P_init;
P_history(:, :, 1, 3) = P_init;

args_ukf = {F_ext, m, k, c};
my_dynamics_ukf = @(t, X) springMassDamper_step(t, X, dt, @springMassDamperDynamics, args_ukf);

for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);
    %%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1);];
    if i <= 5000
        Y_vec = X_sol(i, 1).';
    else
        Y_vec = [];
    end
    tic;
    [X_next, ~, P, UKF_opts] = unscentedKalmanFilter(my_dynamics_ukf, t, X_vec, Y_vec, UKF_opts);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    e1(i+1, 1) = rms(X_sol(i, :).' - X_next);
    normP(i+1, 1) = norm(P);
    P_history(:, :, i+1, 1) = P;

    %%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2);];
    if i <= 5000
        Y_vec = X_sol(i, 1).';
    else
        Y_vec = [];
    end
    tic;
    [X_next, ~, P, EKF_opts] = extendedKalmanFilter(my_dynamics_ukf, t, X_vec, Y_vec, EKF_opts);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    e1(i+1, 2) = rms(X_sol(i, :).' - X_next);
    normP(i+1, 2) = norm(P);
    P_history(:, :, i+1, 2) = P;

    %%%%%%%%%%%%%% KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 3); v1(i, 3);];
    if i <= 5000
        Y_vec = X_sol(i, 1).';
    else
        Y_vec = [];
    end
    tic;
    [X_next, ~, P, KF_opts] = linearKalmanFilter(t, X_vec, Y_vec, KF_opts);
    time_to_solve(i, 3) = toc;
    x1(i+1, 3) = X_next(1);
    v1(i+1, 3) = X_next(2);
    e1(i+1, 3) = rms(X_sol(i, :).' - X_next);
    normP(i+1, 3) = norm(P);
    P_history(:, :, i+1, 3) = P;
end
%%

% Extract the standard deviations from the covariance matrix
std_x = squeeze(sqrt(P_history(1, 1, :, 1))); % Squeeze converts the 1x1xN matrix to a 1xN vector
std_v = squeeze(sqrt(P_history(2, 2, :, 1)));

std_x_ekf = squeeze(sqrt(P_history(1, 1, :, 2))); % Squeeze converts the 1x1xN matrix to a 1xN vector
std_x_kf = squeeze(sqrt(P_history(1, 1, :, 3)));
% Set the confidence interval multiplier (1 for 95%)
conf_mult = 1;

% Calculate the upper and lower bounds for the confidence interval
upper_bound_x = x1(:, 1) + conf_mult * std_x;
lower_bound_x = x1(:, 1) - conf_mult * std_x;

upper_bound_v = v1(:, 1) + conf_mult * std_v;
lower_bound_v = v1(:, 1) - conf_mult * std_v;

upper_bound_x_ekf = min(10, x1(:, 2) + conf_mult * std_x_ekf);
lower_bound_x_ekf = max(-10, x1(:, 2) - conf_mult * std_x_ekf);

upper_bound_x_kf = min(10, x1(:, 3) + conf_mult * std_x_kf);
lower_bound_x_kf = max(-10, x1(:, 3) - conf_mult * std_x_kf);

% Construct Matrices used for colouring in the confidence bound
fill_T = [t_output_points, fliplr(t_output_points)];
fill_X = [upper_bound_x.', fliplr(lower_bound_x.')];
fill_V = [upper_bound_v.', fliplr(lower_bound_v.')];

fill_X_ekf = [upper_bound_x_ekf.', fliplr(lower_bound_x_ekf.')];
fill_X_kf = [upper_bound_x_kf.', fliplr(lower_bound_x_kf.')];

figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1(:, 1))
hold on
plot(t_output_points, upper_bound_x, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Unscented Kalman Filter')
legend('Predicted', '', '', 'Error Band')
ylim([-1, 1])

nexttile
plot(t_output_points, x1(:, 2))
hold on
plot(t_output_points, upper_bound_x_ekf, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x_ekf, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X_ekf, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Extended Kalman Filter')
ylim([-1, 1])
legend('Predicted', '', '', 'Error Band')

nexttile
plot(t_output_points, x1(:, 3))
hold on
plot(t_output_points, upper_bound_x_kf, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x_kf, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X_kf, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Linear Kalman Filter')
ylim([-1, 1])
legend('Predicted', '', '', 'Error Band')


%%
lgd = {'UKF', 'EKF', 'KF'};
figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1)
hold on
plot(t_output_points, upper_bound_x, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Displacement [m]')
% ylim([-0.5, 0.5])

nexttile
plot(t_output_points, v1)
hold on
plot(t_output_points, upper_bound_v, 'r--', 'LineWidth', 1);
plot(t_output_points, lower_bound_v, 'r--', 'LineWidth', 1);
fill(fill_T, fill_V, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
% ylim([-2, 2])

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
% ------------------------------------------------------------
%        Parameter Estimation through state augmentation
% ------------------------------------------------------------
clear;close all;clc
% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 500;%0.37;     % Spring stiffness (N/m)
c = 5;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7.5;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
% F = 1 * sin(0.5 * pi * t_output_points) .* (t_output_points <= 5) + randn(size(t_output_points)) * 0.1;
% % % Pulse parameters
% % pulse_height = 1000; % Amplitude of each pulse
% % pulse_width = 0.05; % Duration of each pulse
% % pulse_separation = 2.5; % Time between the start of pulses
% % 
% % % Generate the pulse train
% % F = zeros(size(t_output_points));
% % for t = 0:pulse_separation:t_final - pulse_separation
% %     % Find the indices for the current pulse
% %     pulse_indices = (t_output_points >= t) & (t_output_points < t + pulse_width);
% % 
% %     % Apply the pulse height to the corresponding indices
% %     F(pulse_indices) = pulse_height;
% % end

% Sinusoidal sweep parameters
F_amplitude = 100;    % Maximum force amplitude
f_start = 0;      % Starting frequency in Hz
f_end = 10.0;        % Ending frequency in Hz

% Generate the sinusoidal sweep
F = F_amplitude * sawtooth(2 * pi * (f_start + (f_end - f_start) * (t_output_points / t_final) / 2) .* t_output_points);

F_ext = @(t) interp1(t_output_points, F, t, "linear");


x1 = zeros(length(t_output_points), 3);
v1 = zeros(size(x1));
e1 = zeros(size(x1));
normP = zeros(size(x1));
M = zeros(size(x1));
K = zeros(size(x1));
C = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), size(x1, 2));

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.0;    % Initial displacement of mass
v1(1, :) = 0.0;    % Initial velocity of mass

M(1, :) = 1;
K(1, :) = 500;
C(1, :) = 5;


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
[~, X_sol] = ode23tb(my_dynamics,t_output_points , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)

X_sol_noisy = X_sol + 0.1 * randn(size(X_sol));
%%
numStates = 5; % Displacement and velocity
numStatesMeasure = 2; % Displacement and velocity
args_ukf = {F_ext};

% Initialise error with error at initial time
X_vec = [x1(1, 1); v1(1, 1); M(1, 1); K(1, 1); C(1, 1); ];
X_pred_init = springMassDamperParamEst(t_init, X_vec, dt, @springMassDamperDynamics, args_ukf);

X_sol_init = [X_sol_noisy(1, :).'; M(1, 1); K(1, 1); C(1, 1);];

P_squared = (X_sol_init - X_pred_init) * (X_sol_init - X_pred_init).' + eye(numStates) * dt; % add identity matrix to stabilise P_init
P_init = chol(P_squared); % Implementation of SR-UKF

UKF_opts = struct('ProcessNoiseCov', ones(numStates) * dt, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', P_init, ...
                  'dt', 1/fs_output, ...
                  'NumParticles', numStates);

EKF_opts = struct('ThreshScal', sqrt(eps) * ones(numStates, 1),...
                  'facF', [], ...
                  'facH', [], ...
                  'ProcessNoiseCov', ones(numStates) * dt, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', P_init, ...
                  'dt', 1/fs_output, ...
                  'NumParticles', numStates);
A = [0,1,0,0,0;
    -k/m, -c/m,0,0,0;
    0,0,1,0,0;
    0,0,0,1,0;
    0,0,0,0,1];
B = [0; 1/m; 0; 0; 0];
Phi = expm(A * dt);
Gamma = A \ ( ( Phi - eye(size(Phi)) ) * B) ;

KF_opts = struct('StateTransition', Phi, ...
                  'ControlMatrix', Gamma,...
                  'Forcing', F_ext,...
                 'MeasurementModel', eye(numStatesMeasure, numStates), ... 
                 'ProcessNoiseCov', ones(numStates) * dt, ...
                 'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                 'ErrorCov', P_init, ...
                 'dt', 1/fs_output, ...
                 'NumParticles', numStates);

P_history = zeros(numStates, numStates, length(t_output_points), size(x1, 2));
P_history(:, :, 1, 1) = P_init;
P_history(:, :, 1, 2) = P_init;
P_history(:, :, 1, 3) = P_init;

args_ukf = {F_ext, m, k, c};
my_dynamics_ukf = @(t, X) springMassDamperParamEst(t, X, dt, @springMassDamperDynamics, args_ukf);
%%
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);
    if i <= 5000
        Y_vec = X_sol_noisy(i, 1).';
    else
        Y_vec = [];
    end
    %%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1); M(i, 1); K(i, 1); C(i, 1)];
    tic;
    [X_next, ~, P, UKF_opts] = unscentedKalmanFilter(my_dynamics_ukf, t, X_vec, Y_vec, UKF_opts);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    e1(i+1, 1) = rms(X_sol_noisy(i, :).' - X_next(1:2));
    normP(i+1, 1) = norm(P);
    M(i+1, 1) = X_next(3);
    K(i+1, 1) = X_next(4);
    C(i+1, 1) = X_next(5);
    P_history(:, :, i+1, 1) = P;

    %%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2); M(i, 2); K(i, 2); C(i, 2)];
    tic;
    [X_next, ~, P, EKF_opts] = extendedKalmanFilter(my_dynamics_ukf, t, X_vec, Y_vec, EKF_opts);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    e1(i+1, 2) = rms(X_sol_noisy(i, :).' - X_next(1:2));
    normP(i+1, 2) = norm(P);
    M(i+1, 2) = X_next(3);
    K(i+1, 2) = X_next(4);
    C(i+1, 2) = X_next(5);
    P_history(:, :, i+1, 2) = P;

    %%%%%%%%%%%%%% KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 3); v1(i, 3); M(i, 3); K(i, 3); C(i, 3)];
    A = [0,1,0,0,0;
        -K(i, 3)/M(i, 3), -C(i, 3)/M(i, 3),0,0,0;
        0,0,1,0,0;
        0,0,0,1,0;
        0,0,0,0,1];
    B = [0; 1/M(i, 3); 0; 0; 0];
    Phi = expm(A * dt);
    
    if cond(A) <= 1e12
        Gamma = A \ ( ( Phi - eye(size(Phi)) ) * B) ;
    else
        Gamma = pinv(A) * ( ( Phi - eye(size(Phi)) ) * B);
    end
    KF_opts.ControlMatrix = Gamma;
    KF_opts.StateTransition = Phi;

    tic;
    [X_next, ~, P, KF_opts] = linearKalmanFilter(t, X_vec, Y_vec, KF_opts);
    time_to_solve(i, 3) = toc;
    x1(i+1, 3) = X_next(1);
    v1(i+1, 3) = X_next(2);
    e1(i+1, 3) = rms(X_sol_noisy(i, :).' - X_next(1:2));
    normP(i+1, 3) = norm(P);
    M(i+1, 3) = X_next(3);
    K(i+1, 3) = X_next(4);
    C(i+1, 3) = X_next(5);
    P_history(:, :, i+1, 3) = P;
end
%%
figure
T = tiledlayout('horizontal');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
hold on
plot(X_sol(1:5000, 1), X_sol(1:5000, 2), 'k.')
plot(x1(1:5000, 1), v1(1:5000, 1), '.')
grid on
xlabel('Displacement [m]')
ylabel('Velocity [m/s]')
legend('Measured','UKF')

nexttile
hold on
plot(X_sol(1:5000, 1), X_sol(1:5000, 2), 'k.')
plot(x1(1:5000, 2), v1(1:5000, 2), '.')
grid on
xlabel('Displacement [m]')
ylabel('Velocity [m/s]')
legend('Measured','EKF')

nexttile
hold on
plot(X_sol(1:5000, 1), X_sol(1:5000, 2), 'k.')
plot(x1(1:5000, 3), v1(1:5000, 3), '.')
grid on
xlabel('Displacement [m]')
ylabel('Velocity [m/s]')
legend('Measured','KF')

%%
% Extract the standard deviations from the covariance matrix
std_x = squeeze(sqrt(P_history(1, 1, :, 1))); % Squeeze converts the 1x1xN matrix to a 1xN vector
std_v = squeeze(sqrt(P_history(2, 2, :, 1)));

std_x_ekf = squeeze(sqrt(P_history(1, 1, :, 2))); % Squeeze converts the 1x1xN matrix to a 1xN vector
std_x_kf = squeeze(sqrt(P_history(1, 1, :, 3)));
% Set the confidence interval multiplier (1 for 95%)
conf_mult = 1;

% Calculate the upper and lower bounds for the confidence interval
upper_bound_x = x1(:, 1) + conf_mult * std_x;
lower_bound_x = x1(:, 1) - conf_mult * std_x;

upper_bound_v = v1(:, 1) + conf_mult * std_v;
lower_bound_v = v1(:, 1) - conf_mult * std_v;

upper_bound_x_ekf = min(10, x1(:, 2) + conf_mult * std_x_ekf);
lower_bound_x_ekf = max(-10, x1(:, 2) - conf_mult * std_x_ekf);

upper_bound_x_kf = min(10, x1(:, 3) + conf_mult * std_x_kf);
lower_bound_x_kf = max(-10, x1(:, 3) - conf_mult * std_x_kf);

% Construct Matrices used for colouring in the confidence bound
fill_T = [t_output_points, fliplr(t_output_points)];
fill_X = [upper_bound_x.', fliplr(lower_bound_x.')];
fill_V = [upper_bound_v.', fliplr(lower_bound_v.')];

fill_X_ekf = [upper_bound_x_ekf.', fliplr(lower_bound_x_ekf.')];
fill_X_kf = [upper_bound_x_kf.', fliplr(lower_bound_x_kf.')];

figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1(:, 1))
hold on
plot(t_output_points, upper_bound_x, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Unscented Kalman Filter')
legend('Predicted', '', '', 'Error Band')
ylim([-1, 1])

nexttile
plot(t_output_points, x1(:, 2))
hold on
plot(t_output_points, upper_bound_x_ekf, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x_ekf, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X_ekf, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Extended Kalman Filter')
ylim([-1, 1])
legend('Predicted', '', '', 'Error Band')

nexttile
plot(t_output_points, x1(:, 3))
hold on
plot(t_output_points, upper_bound_x_kf, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x_kf, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X_kf, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Linear Kalman Filter')
ylim([-1, 1])
legend('Predicted', '', '', 'Error Band')


%%
lgd = {'UKF', 'EKF', 'KF'};
figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
plot(t_output_points, x1)
hold on
plot(t_output_points, upper_bound_x, 'k--', 'LineWidth', 1);
plot(t_output_points, lower_bound_x, 'k--', 'LineWidth', 1);
fill(fill_T, fill_X, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Displacement [m]')
% ylim([-0.5, 0.5])

nexttile
plot(t_output_points, v1)
hold on
plot(t_output_points, upper_bound_v, 'r--', 'LineWidth', 1);
plot(t_output_points, lower_bound_v, 'r--', 'LineWidth', 1);
fill(fill_T, fill_V, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on
legend(lgd)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
% ylim([-2, 2])

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

figure
T = tiledlayout('vertical');
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile; plot(t_output_points, M); legend('UKF', 'EKF', 'KF');
nexttile; plot(t_output_points, K); legend('UKF', 'EKF', 'KF');
nexttile; plot(t_output_points, C); legend('UKF', 'EKF', 'KF')

%%
clear;close all;clc
% --- 1. Define System Parameters ---
m = 1;%7.64e-10;      % Mass (kg)
k = 500;%0.37;     % Spring stiffness (N/m)
c = 5;%1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7.5;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);

% Define "experimental" forcing function to represent real world data
% F = 1 * sin(0.5 * pi * t_output_points) .* (t_output_points <= 5) + randn(size(t_output_points)) * 0.1;
% % % Pulse parameters
% % pulse_height = 1000; % Amplitude of each pulse
% % pulse_width = 0.05; % Duration of each pulse
% % pulse_separation = 2.5; % Time between the start of pulses
% % 
% % % Generate the pulse train
% % F = zeros(size(t_output_points));
% % for t = 0:pulse_separation:t_final - pulse_separation
% %     % Find the indices for the current pulse
% %     pulse_indices = (t_output_points >= t) & (t_output_points < t + pulse_width);
% % 
% %     % Apply the pulse height to the corresponding indices
% %     F(pulse_indices) = pulse_height;
% % end

% Sinusoidal sweep parameters
F_amplitude = 100;    % Maximum force amplitude
f_start = 0;      % Starting frequency in Hz
f_end = 10.0;        % Ending frequency in Hz

% Generate the sinusoidal sweep
F = F_amplitude * sawtooth(2 * pi * (f_start + (f_end - f_start) * (t_output_points / t_final) / 2) .* t_output_points);

F_ext = @(t) interp1(t_output_points, F, t, "linear");


x1 = zeros(length(t_output_points), 1);
v1 = zeros(size(x1));
a1 = zeros(size(x1));
e1 = zeros(size(x1));
normP = zeros(size(x1));

time_to_solve = zeros(length(t_output_points), size(x1, 2));

% --- 3. Define Initial State ---
% The state vector X is [x; v] (displacement; velocity)
x1(1, :) = 0.0;    % Initial displacement of mass
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
my_dynamics = @(t, X) springMassDamperDynamicsExpanded(t, X, args{:});
tic;
[~, X_sol] = ode23tb(my_dynamics,t_output_points , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)

X_sol = X_sol + 0.1 * randn(size(X_sol));
%%
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off");
numStates = 2; % Displacement and velocity
numStatesMeasure = 1; % Displacement and velocity
args_ukf = {F_ext, m, k, c};

% Initialise error with error at initial time
X_vec = [x1(1, 1); v1(1, 1)];
X_pred_init = springMassDamper_step(t_init, X_vec, dt, @springMassDamperDynamics, args_ukf);
X_sol_init = X_sol(1, :).';
P_init = (X_sol_init - X_pred_init) * (X_sol_init - X_pred_init).' + eye(numStates) * dt; % add identity matrix to stabilise P_init

UKF_opts = struct('ProcessNoiseCov', eye(numStates) * dt, ...
                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                  'ErrorCov', P_init, ...
                  'dt', 1/fs_output, ...
                  'NumParticles', numStates);

NumParams = 3;
P_theta = 1e-6 * eye(NumParams);
UKF_theta_opts = struct('ProcessNoiseCov', eye(NumParams) * dt, ...
                        'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
                        'ErrorCov', P_theta, ...
                        'dt', 1/fs_output, ...
                        'NumParticles', NumParams);

% % EKF_opts = struct('ThreshScal', sqrt(eps) * ones(numStates, 1),...
% %                   'facF', [], ...
% %                   'facH', [], ...
% %                   'ProcessNoiseCov', ones(numStates) * dt, ...
% %                   'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
% %                   'ErrorCov', P_init, ...
% %                   'dt', 1/fs_output, ...
% %                   'NumParticles', numStates);
% % A = [0,1;-k/m, -c/m];
% % B = [0;1/m];
% % Phi = expm(A * dt);
% % Gamma = A \ ( ( Phi - eye(size(Phi)) ) * B) ;
% % 
% % KF_opts = struct('StateTransition', Phi, ...
% %                   'ControlMatrix', Gamma,...
% %                   'Forcing', F_ext,...
% %                  'MeasurementModel', eye(numStatesMeasure, numStates), ... 
% %                  'ProcessNoiseCov', ones(numStates) * dt, ...
% %                  'MeasurementNoiseCov', eye(numStatesMeasure) * 0.1, ...
% %                  'ErrorCov', P_init, ...
% %                  'dt', 1/fs_output, ...
% %                  'NumParticles', numStates);

P_history = zeros(numStates, numStates, length(t_output_points), size(x1, 2));
for i = 1:length(size(x1, 2))
    P_history(:, :, 1, i) = P_init;
end

args_ukf = {F_ext, m, k, c};
my_dynamics_ukf = @(t, X) springMassDamper_step(t, X, dt, @springMassDamperDynamics, args_ukf);

theta_vec = [m; k; c]; 
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);
    %%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1);];
    if i <= 5000
        Y_vec = X_sol(i, 1).';
    else
        Y_vec = [];
    end
    % ----------- State Estimation ---------------------
    tic;
    [X_next, ~, P, UKF_opts] = unscentedKalmanFilter(my_dynamics_ukf, t, X_vec, Y_vec, UKF_opts);
    V_pred = UKF_opts.V_pred;
    % ----------- Parameter Estimation -----------------
        % Use the innovation from the state filter as the measurement
    if ~isempty(V_pred)
        Y_theta_vec = V_pred;
    
        % The measurement for the parameter filter is the innovation
        % The dynamics are assumed to be a random walk
        [theta_next, ~, P_theta, UKF_theta_opts] = unscentedKalmanFilter(@parameterDynamics, t, theta_vec, Y_theta_vec, UKF_theta_opts);

        % Update the parameter state for the next time step
        theta_vec = theta_next;
        % Save the updated parameter covariance
        P_theta_history(:, :, i+1) = P_theta;
    else
        % No innovation means no measurement for the parameter filter
        % So only predict the parameters forward (which is just a random walk)
        [theta_next, ~, P_theta, UKF_theta_opts] = unscentedKalmanFilter(@parameterDynamics, t(i), theta_vec, [], UKF_theta_opts);
        
        theta_vec = theta_next;
        P_theta_history(:, :, i+1) = P_theta;
    end
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    e1(i+1, 1) = rms(X_sol(i, :).' - X_next);
    normP(i+1, 1) = norm(P);
    P_history(:, :, i+1, 1) = P;

end
%%
function dX = springMassDamperDynamics(t, X, Forcing, m, k, c)

    dx = X(2, :);
    dvx = (Forcing(t) - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end

function theta_next = parameterDynamics(t, theta)
    theta_next = theta; % Parameters assumed constant
end

function dX = springMassDamperDynamicsExpanded(t, X, Forcing, m, k, c)
    
    if t <= 5
        F = Forcing(t);
    else
        F = 1 * sin(10 * 2 * pi * t);
    end  
    dx = X(2, :);
    dvx = (F - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end

function states = springMassDamper_step(t, X_vec, dt, func, args)
    % Extract states to integrate
    X_states = X_vec(:);
    X_next = evaluateTRBDF2(func, dt, t, X_states, args);
    % Augment output states
    states = X_next(:);
end

function states = springMassDamperParamEst(t, X_vec, dt, func, args)

    % t_target = t + dt;
    % hRK_current = dt;

    % Extract states to integrate
    X_states = X_vec(1:2, :);
    % Update parameters of system
    args{2} = X_vec(3);
    args{3} = X_vec(4);
    args{4} = X_vec(5);

    % while t < t_target
    %     hRK_current = min(hRK_current,  t_target - t);
    %     % Apply integration method
    %     [X_next, h_next] = adaptive_ODE(func, dt, t, X_states, args);
    %     % Update current time based on the step taken
    %     t = t + hRK_current;
    %     % Update time step
    %     hRK_current = h_next;
    %     % Update solution
    %     X_states = X_next;
    % end

    X_next = evaluateTRBDF2(func, dt, t, X_states, args);
    % Augment output states
    [m, k, c] = deal(args{2:end});
    states = [X_next(:);
              m;
              k;
              c];
end