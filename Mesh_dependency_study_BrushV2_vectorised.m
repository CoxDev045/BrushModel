clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
% if isempty(gcp('nocreate'))
%     pool = parpool('Processes');
% else
%     delete(gcp("nocreate"))
%     pool = parpool('Processes');
% end

%%

clearvars; clc;
close all;

total_time = tic;
pressure_start = tic;
fprintf("Starting simulation. \n")

% Contact Width Parameters
a1 = 1.03;
a2 = 0.44;
a3 = 0;

% Vertical Load
Fz = 100;

% Contact Width (y)
a = a1 * Fz^a2 + a3;
% Contact Length (x)
b = 9;

contact_area = a * b;

n_values = [3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80];
SR = -0.3:0.005:0.3;

re = 30;
alpha = deg2rad(0);
fs = 4.71e4; % Hz
dt = 1 / fs;
omega_z = 0;

t_initial = 0;
t_final = 1;
time_span = linspace(t_initial, t_final, t_final * fs);

omega = [linspace(0, 10, length(time_span)/2), repmat(10, 1, length(time_span)/2)];

% Initialise Pressure function coefficients
nx = 2;
ny = 1;
n_params = [nx, ny];

lambda_x = 1;
lambda_y = 1;
lambda = [lambda_x, lambda_y];

xe = 0; % mm
ye = 0; % mm
COP = [xe, ye];

One_D = false;

for sr_idx = 1:numel(SR)
    v0 = omega * re / (SR(sr_idx) + 1);
    
    for k = 1:numel(n_values)
        numBrushes = n_values(k);
        contact_area_per_brush = contact_area / (numBrushes^2);
        
        x_vals = linspace(-a, a, numBrushes);
        y_vals = linspace(-b, b, numBrushes);
        [X, Y] = meshgrid(x_vals, y_vals);
        
        % Compute 2D pressure distribution
        [~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
        P_grid = real(P_grid);
        fprintf('Pressure Grid Successfully calculated in %.6f \n', toc(pressure_start));
        
        % Preallocate simulation results
        sim_solution = zeros(numBrushes, numBrushes, length(time_span), 9);
        shift_amount_cumulative = cumsum(omega * dt * re);
        shift_amount_cumulative = floor(shift_amount_cumulative);
        
        for i = 1:length(time_span)
            % Shift pressure distribution
            sim_solution(:, :, i, 1) = circshift(P_grid, [-shift_amount_cumulative(i), 0]);
        end
        
        brushArray = BrushVec(X(:), Y(:), P_grid(:), zeros(1, numBrushes^2));
        
        start_brush_sim = tic;
        fprintf('Starting brush model simulation! \n');
        
        for i = 1:length(time_span)
            %%%%%%%%%%%%%% Update Pressure %%%%%%%%%%%%%%%%%%%%
            tempPress = sim_solution(:, :, i, 1);
            brushArray.press = tempPress(:);

            %%%%%%%%%%%%%% Update Pressure Dependent Properties %%%%%%%%%%%%%%%%
            brushArray.ky = brushArray.ky_0 + brushArray.ky_l .* brushArray.press;
            brushArray.kx = brushArray.phi .* brushArray.ky;

            %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
            brushArray = brushArray.update_brush(omega(i), omega_z, re, v0(i), alpha, dt);

            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution(:, :, i, 4) = reshape(brushArray.delta_x, numBrushes, numBrushes);
            sim_solution(:, :, i, 5) = reshape(brushArray.delta_y, numBrushes, numBrushes);
            sim_solution(:, :, i, 6) = reshape(brushArray.tauX, numBrushes, numBrushes);
            sim_solution(:, :, i, 7) = reshape(brushArray.tauY, numBrushes, numBrushes);
            sim_solution(:, :, i, 8) = reshape(brushArray.slide, numBrushes, numBrushes);
            sim_solution(:, :, i, 9) = reshape(brushArray.mu, numBrushes, numBrushes);

            %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
            sim_solution(:, :, i, 2) = hypot(sim_solution(:, :, i, 4), sim_solution(:, :, i, 5));
            sim_solution(:, :, i, 3) = hypot(sim_solution(:, :, i, 6), sim_solution(:, :, i, 7));
        end
        
        fprintf('Finished Simulation in %.2fs \n', toc(start_brush_sim));
        
        % Compute force integration
        dx = x_vals(2) - x_vals(1);
        dy = y_vals(2) - y_vals(1);
        dA = dx * dy;
        
        pressure_mask = sim_solution(:, :, :, 1) > 1e-3;
        total_valid_points = sum(sum(pressure_mask, 1), 2);
        
        forceX = squeeze(trapz(trapz(sim_solution(:, :, :, 6) .* pressure_mask, 1), 2)) * dA;
        forceY = squeeze(trapz(trapz(sim_solution(:, :, :, 7) .* pressure_mask, 1), 2)) * dA;
        forceTotal = squeeze(trapz(trapz(sim_solution(:, :, :, 3) .* pressure_mask, 1), 2)) * dA;
        avg_mu = squeeze(trapz(trapz(sim_solution(:, :, :, 9) .* pressure_mask, 1), 2) ./ max(total_valid_points));
        
        % Save results
        filename = sprintf('Slip_dependency_study/Slip_dependency_study_n%d_SR%.4f.mat', n_values(k), SR(sr_idx));
        save(filename, 'SR', 'time_span', 'forceX', 'forceY', 'forceTotal', 'avg_mu', 'omega', 'v0');
        fprintf('Saved results for SR = %.2f to %s\n', SR(sr_idx), filename);
    end
end

fprintf('Total Simulation time: %.2fs \n', toc(total_time));



%%
close all; clc;

fprintf("Starting simulation. \n")
tic;
% Test with a 10x10 grid of brushes
a1 = 1.03;
a2 = 0.44;
a3 = 0;

Fz = 100;

a = a1 * Fz^a2 + a3;
b = 9;
contact_area = a * b;

re = 200; % mm
alpha = deg2rad(0);
max_iters = 100;
tol = 1e-4;
fs = 1000; % Hz
dt = 1 / fs;
omega_z = 0;

t_initial = 0;
t_final = 1;
time_span = linspace(t_initial, t_final, t_final * fs);

omega = 0.01 * ones(1,length(time_span));
omega(length(time_span)/4+1:length(time_span)/2) = linspace(0.01, 1, length(time_span)/4);
% omega(length(time_span)/2+1:2*length(time_span)/2) = linspace(0.1, 1, length(time_span)/2);
omega(length(time_span)/2+1:length(time_span)) = 1;

% SR = (omega * re - v0) / v0;
SR = -0.3:0.005:0.3;

% Grid Sizes
n = [3, 5, 10, 20, 50];%, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000].';

% Initialize results struct array with empty fields for each SR value
results(numel(SR)) = struct('SR_val', [], 'n_values', [], 'time_span', [], ...
                            'time', [], 'force', [], 'avg_mu', [], 'omega', [], 'v0', []);


% Initialise Pressure function coefficients
nx = 4;
ny = 2;
n_params = [nx, ny];
lambda_x = 1;
lambda_y = 1;
lambda = [lambda_x, lambda_y];
xe = 0; % mm
ye = 0; % mm
COP = [xe, ye];
One_D = false;

for sr_idx = 1:numel(SR)
    SR_val = SR(sr_idx);
    v0 = omega * re / (SR_val + 1);

    % Temporary storage for each grid size n at this SR value
    time_data = cell(1, numel(n));
    force_data = cell(1, numel(n));
    avg_mu_data = cell(1, numel(n));

    for k = 1:numel(n)
        contact_area_per_brush = contact_area / (n(k) * n(k));
        x_vals = linspace(-a, a, n(k));
        y_vals = linspace(-b, b, n(k));
        
        [X, Y] = meshgrid(x_vals, y_vals);
        [~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
        P_grid = real(P_grid);
        
        tau_x_grid = zeros(n(k), n(k), length(time_span));
        mu_grid = zeros(n(k), n(k), length(time_span));

        time = zeros(length(time_span), 1);
        force = zeros(length(time_span), 1);
        avg_mu = zeros(length(time_span), 1);
        
        tic;  % Start timing
        for q = 1:length(time_span)
            tic;  % Start timing
            for i = 1:n(k)
                for j = 1:n(k)
                    brush = Brush(x_vals(i), y_vals(j), P_grid(i, j), Fz);
                    
                    if abs(P_grid(i, j)) > 1e-3
                        while brush.x >= -a
                            brush = brush.update_brush(omega(q), omega_z, re, v0(q), alpha, dt, tol, max_iters);
                        end
                    end
                    
                    tau_x_grid(i, j, q) = brush.tauX;
                    mu_grid(i, j, q) = brush.mu;
                end
            end
            time(q) = toc;  % Stop timing 
            avg_mu(q) = mean(mu_grid(:, :, q), 'all');
            force(q) = trapz(trapz(tau_x_grid(:, :, q))) * contact_area_per_brush;
        end
        fprintf('Finished Simulation in %.2fs \n', toc);

        time_data{k} = time;
        force_data{k} = force;
        avg_mu_data{k} = avg_mu;
    end

    % Store results for this SR value
    results(sr_idx).SR_val = SR_val;
    results(sr_idx).n_values = n;
    results(sr_idx).time_span = time_span;
    results(sr_idx).time = time_data;
    results(sr_idx).force = force_data;
    results(sr_idx).avg_mu = avg_mu_data;
    results(sr_idx).omega = omega;
    results(sr_idx).v0 = v0;
end
toc
% Save results to file once all simulations complete
save('Slip_dependency_study.mat', 'results');

