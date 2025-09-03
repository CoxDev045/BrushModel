clear;close all;clc

R = 0:.002:1;
TH = 2*pi*(0:.002:1); 
X = R'*cos(TH);
Y = R'*sin(TH);
t_val = 1;

Z = rosenbrock(t_val, [X(:),Y(:)]);
Z = reshape(Z,size(X));

figure
surf(X, Y, Z)
shading flat

options.JacTol = 1e-3;
options.ThreshScal = sqrt(eps) * [1;1];
options.MaxIters = 1000;
options.FunTol = 1e-6;

x0 = rand(2, 1);
objectiveFunc = @(t, x) rosenbrock(t, x);

[x_next, Fty_next, exitflag, output] = solveTrustRegionDogLeg(objectiveFunc, t_val, x0, options);

%%
clear;close all;clc
% --- 1. SETUP THE PROBLEM ---
% Define system parameters
m = 10;
c = 2;
k = 100;
% Create a forcing function (time-dependent)
forcing_func = @(t) 5 * sin(2 * pi * t);
args = {m, c, k, forcing_func};

% Define a wrapper function for the solver
% The wrapper must be in the format F(t, y, args), as required by numjac
% The solveTrustRegionDogLeg function is designed to minimize F(y), not F(t, y)
% so the wrapper must not include t
my_dynamics_solver = @(t, y) springMassDamperDynamics(t, y, args{:});

% Set up solver options
options.JacTol = 1e-3;
options.ThreshScal = sqrt(eps) * [1;1];
options.MaxIters = 1000;
options.FunTol = 1e-6;

% --- 2. DEFINE TIME POINTS & INITIAL GUESS ---
t_points = 0:0.1:10;
y_solution = zeros(length(t_points), 2);

% Initial guess for the first time point (t=0)
y0 = [0; 0]; 

% --- 3. LOOP THROUGH TIME POINTS ---
for i = 1:length(t_points)
    t_val = t_points(i);
    
    % Call the solver to find the root at the current time
    [y_root, Fval, exitflag, output] = solveTrustRegionDogLeg(my_dynamics_solver,t_val, y0, options);
    
    % Store the result
    y_solution(i, :) = y_root';
    J_solution(:, :, i) = output.Jacobian;
    
    % Update the initial guess for the next time step
    % This is a good strategy as the solution from the current time
    % is likely close to the solution at the next time.
    y0 = y_root;
end

% --- 4. VISUALIZE THE RESULTS ---
figure;
plot(t_points, y_solution(:, 1));
xlabel('Time');
ylabel('Equilibrium Position (x)');
title('Time-Varying Equilibrium Position of a SDOF System');
grid on;

% Add a plot for the forcing function for comparison
hold on;
plot(t_points, arrayfun(forcing_func, t_points) / k, 'r--'); % x_eq = F(t)/k for this system
legend('Calculated Equilibrium Position', 'Analytical Equilibrium Position (F(t)/k)');
hold off;

%%
clear;close all;clc
% --- 1. Define System Parameters ---
m = 7.64e-10;      % Mass (kg)
k = 0.37;     % Spring stiffness (N/m)
c = 1.40e-4;    % Damping coefficient (Ns/m)

% --- 2. Define Simulation Time Span and Output Points ---
t_init = 0;             % Start time (s)
t_final = 7;           % End time (s)
fs_output = 1000;       % Desired output sampling frequency (Hz)
                        % This determines how many points ode45 returns.
dt = 1/fs_output;
t_output_points = linspace(t_init, t_final, t_final * fs_output);
% Define "experimental" forcing function to represent real world data
F = 0.1 * sin(2 * pi * t_output_points) .* (t_output_points <= 5) + randn(size(t_output_points)) * 0.001;
forcing_func = @(t) interp1(t_output_points, F, t, "linear", "extrap");

args = {m, c, k, forcing_func};

x1 = zeros(length(t_output_points), 3);
v1 = zeros(size(x1));
time_to_solve = zeros(size(x1));

x1(1, :) = 0.1;    % Initial displacement of mass
v1(1, :) = 0.0;    % Initial velocity of mass

options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","off"); 
for i = 1:length(t_output_points)-1
    t = t_output_points(i);
    t_target = t_output_points(i+1);

    %%%%%%%%%%%%%% TR-BDF2 (Own Solver) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 1); v1(i, 1);];
    tic;
    X_next = evaluateTRBDF2(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 1) = toc;
    x1(i+1, 1) = X_next(1);
    v1(i+1, 1) = X_next(2);
    % e1(i+1) = rms(X_sol(i, :).' - X_next);


    %%%%%%%%%%%%%% TR-BDF2 (fsolve) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 2); v1(i, 2);];
    tic;
    X_next = trbdf2(@springMassDamperDynamics, dt, t, X_vec, args);
    time_to_solve(i, 2) = toc;
    x1(i+1, 2) = X_next(1);
    v1(i+1, 2) = X_next(2);
    % e1(i+1) = rms(X_sol(i, :).' - X_next);

    %%%%%%%%% ODE23tb at each time step %%%%%%%%%%%%%%%%%%%%%%%%
    X_vec = [x1(i, 3); v1(i, 3);];
    tic;
    [~, X_next] = ode23s(@(t, X) springMassDamperDynamics(t, X, args{:}),[t, t_target] , X_vec, options);
    time_to_solve(i, 3) = toc;
    x1(i+1, 3) = X_next(end, 1);
    v1(i+1, 3) = X_next(end, 2);
    % % e1(i+1, 10) = rms( X_sol(i, :) - X_next(end, :) );
    % e1(i+1, 1) = rms(x1(i+1, 1) - x1(i+1, 2) );
    % e1(i+1, 2) = rms(v1(i+1, 1) - v1(i+1, 2) );
    % e1(i+1, 2) = rms(v1(i+1, 1) - v1(i+1, 2) );
    % 
    
end


lgd = {'TR-BDF2 (Own Solver)', 'TR-BDF2 (fsolve)', 'ode23tb'};
figure
T = tiledlayout('vertical');

nexttile
plot(t_output_points, x1)
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
legend(lgd)

nexttile
plot(t_output_points, v1)
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend(lgd)

nexttile
semilogy(t_output_points, time_to_solve)
grid on
xlabel('Time [s]')
ylabel('log_{10}(CPU time) [s]')
legend(lgd)

%%
function y = rosenbrock(t, x, varargin)
[~, M] = size(x);
if M ~= 2
    % Input is column vector
    xp = x(:).';
else
    xp = x;
end

% y = 100 * (xp(:, 2) - xp(:, 1).^2).^2 + (1 - xp(:, 1)).^2;
y = xp(:, 1).^2 + xp(:, 2).^2;
end

function dX = springMassDamperDynamics(t, X, m, k, c, Forcing)

    dx = X(2, :);
    dvx = (Forcing(t) - k .* X(1, :) - c .* X(2, :)) ./m;
    dX = [dx; dvx];
end