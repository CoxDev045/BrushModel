clearvars; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal')
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

numBrushes = 10;

re = 36;
alpha = deg2rad(0);
fs = 1e6; % Hz
dt = 1 / fs;
omega_z = 0;

t_final = 2;
SR = 0.1;
rpm = 100;

t_initial = 0;
num_steps = t_final * fs;
time_span = linspace(t_initial, t_final, num_steps);

% Angular velocity: ramp up to max RPM for first half, then constant
omega_max = rpm * 2 * pi / 60; % Convert RPM to rad/s
half_steps = floor(num_steps/2);

omega_profile = [linspace(0, omega_max, half_steps), ...
            omega_max * ones(1, num_steps - half_steps)]';

omega = repmat(omega_profile, 1, numBrushes^2);

% Calculate v0 based on slip ratio
v0 = omega * re ./ (SR + 1);

% % [time_span, omega, v0] = setupSimulationParams(t_final, fs, numBrushes, re, SR, rpm);

x_vals = linspace(-a, a, numBrushes);
y_vals = linspace(-b, b, numBrushes);

% Calculate area element for integration
dx = x_vals(2) - x_vals(1);
dy = y_vals(2) - y_vals(1);
dA = dx * dy;

[X, Y] = meshgrid(x_vals, y_vals);

% Initialise PRessure function coefficients
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

% Calculate 2D pressure distribution
[~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
P_grid = real(P_grid);

% Clear unnecessary variables
clear x_vals y_vals dx dy a1 a2 a3 COP One_D contact_area lambda lambda_x lambda_y n_params nx ny xe ye;

sim_solution = zeros(numBrushes * numBrushes, length(time_span), 9); 

% Calculate the amount to shift per timestep
shift_amount = (omega * dt * re); 
shift_amount = cumsum(shift_amount);
% % shift_amount_cumulative = cumsum(shift_amount);

% figure
% hold on
% plot(time_span, shift_amount_cumulative)

shift_amount = floor(shift_amount);

% plot(time_span, shift_amount_cumulative)
% grid on
% legend('Continuous', 'Discrete')
% ylabel('Units Shifted')
% xlabel('Time [s]')
for i = 1:length(time_span)
    shifted_press = circshift(P_grid, [-shift_amount(i), 0]);
    % Shift the pressure distribution downwards and wrap it
    sim_solution(:, i, 1) = shifted_press(:);
end
 
fprintf('Pressure Grid Successfully calculated in %.6f \n', toc(pressure_start));

% Pre-allocate array of empty objects
brushArray = BrushVec_improved(X(:), Y(:), P_grid(:));

% % simParams.t = time_span;
% % simParams.omega = omega;
% % simParams.v0 = v0;

clear pressure_start P_grid shifted_press shift_amount;
% % clear time_span omega v0;
% % whos

%
start_brush_sim = tic;
fprintf('Starting brush model simulation! \n')

for i = 1:length(time_span)
    %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
    brushArray = brushArray.update_brush(sim_solution(:, i, 1), omega(i, :), omega_z, re, v0(i, :), alpha, dt);
    %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sim_solution( :, i, 4) = brushArray.delta_x;
    sim_solution( :, i, 5) = brushArray.delta_y;
    sim_solution( :, i, 6) = brushArray.tauX;
    sim_solution( :, i, 7) = brushArray.tauY;
    sim_solution( :, i, 8) = brushArray.slide;
    sim_solution( :, i, 9) = brushArray.mu;

    %%%%%%%%%%%%%% Wrap x-values if edge is reached %%%%%%%%%%%%%%%%%%%%
    brushArray.x(brushArray.x < -a) = a;
    %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
    sim_solution( :, i, 2) = hypot(sim_solution( :, i, 4), sim_solution( :, i, 5));
    sim_solution( :, i, 3) = hypot(sim_solution( :, i, 6), sim_solution( :, i, 7));
    %%%%%%%%%%%%%% Increment counter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(i, 10000) == 0
        fprintf("%d iterations completed! %d remaining.\n", i, length(time_span) - i);
    end
end

fprintf('Brush Sim ended successfully! in %.6f \n', toc(start_brush_sim))
fprintf('Total Simulation time: %.2fs \n', toc(total_time));
%

sim_solution = fillmissing(sim_solution, 'constant', 0);

pressure_mask = sim_solution( :, :, 1) > 1e-3;
total_valid_points = max(trapz(pressure_mask( :, :), 1)); % Count valid points
%
% Intergrate stress to get force
forceX = squeeze(trapz(sim_solution( :, :, 6) .* pressure_mask)) * dA;

forceY = squeeze(trapz(sim_solution( :, :, 7) .* pressure_mask)) * dA;

forceTotal = squeeze(trapz(sim_solution( :, :, 3) .* pressure_mask)) * dA;

avg_mu = squeeze(trapz(sim_solution( :, :, 9) .* pressure_mask)./ total_valid_points); % Avoid division by zero


colorbar_names = cell(9,1);
colorbar_names{1} = "Pressure Distribution [MPa]";
colorbar_names{2} = "Total Displacement [mm]";
colorbar_names{3} = "Total Shear Stress [MPa]";
colorbar_names{4} = "Relative X Displacement [mm]";
colorbar_names{5} = "Relative Y Displacement [mm]";
colorbar_names{6} = "\tau_x (X Shear Stress) [MPa]";
colorbar_names{7} = "\tau_y (Y Shear Stress) [MPa]";
colorbar_names{8} = "Sliding Condition [Boolean]";
colorbar_names{9} = "Friction Coefficient";


% Moving Average
winSize = 50;
b_val = (1 / winSize) * ones(1, winSize);
a_val = 1;

forceTotal_filt = filtfilt(b_val, a_val, forceTotal);
forceX_filt = filtfilt(b_val, a_val, forceX);
forceY_filt = filtfilt(b_val, a_val, forceY);

figure
subplot(311)
plot(time_span, (forceX));
hold on
plot(time_span, forceX_filt)
grid on;
title("Longitudinal Force [N]")
legend('Simulated', 'MA Filtered')

subplot(312)
plot(time_span, forceY);
hold on
plot(time_span, forceY_filt)
grid on;
title("Lateral Force [N]")
legend('Simulated', 'MA Filtered')

subplot(313)
plot(time_span, forceTotal);
hold on
plot(time_span, hypot(forceX, forceY));
plot(time_span, (avg_mu * Fz));
grid on;
title("Total Force [N]")
legend("Total Force (from stresses) [N]", "Total Force (from forces) [N]","Avg Friction times Vertical Force [N]")
ylim([-10, 100])

[Y_mag2, half_freqz] = One_sided_fft(forceTotal, fs);
N = 2 * length(Y_mag2);

PSD = (1 /(fs * N)) * real(Y_mag2  .* conj(Y_mag2 ));

figure
loglog(half_freqz, PSD)
grid on
xlabel('Frequency [kHz]')
ylabel('Magnitude [N]')
%%

% Reshape solution array to grid before animation
% load('Animation_ready_rolling_20x20_10slip_verlet_solver_6s.mat')
sim_solution = reshape(sim_solution, numBrushes, numBrushes, length(time_span), 9);
% % % sim_solution_ds = downsample(sim_solution, 50);
% % % sim_solution_ds = reshape(ans, numBrushes, numBrushes, length(time_span) / 50, 9);

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

x_vals = linspace(-a, a, numBrushes);
y_vals = linspace(-b, b, numBrushes);

fprintf('Simulation Solution is ready for animation! \n');

% Initialize Video
video_filename = sprintf('BrushSim_%dHz_%drpm_slip%d.mp4', fs, rpm, 100 * 0.1);
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 60;  
open(v);

% Initialize Figure
if ishandle(1)
    close(1);
end
fh = figure(1);
fh.WindowState = 'maximized';
h = gobjects(1, 9);
im = gobjects(1, 9);

% Predefined Color Limits
clim_values = [
    0, 0.5;
    0, 1.0;
    0, 0.2;
    -0.5, 0.5;
    -0.5, 0.5;
    0, 0.2;
    0, 0.2;
    0, 1;
    0, 1.3
];

% Initialize subplots
for j = 1:9
    h(j) = subplot(3, 3, j);
    im(j) = imagesc(y_vals, x_vals, sim_solution(:, :, 1, j));
    c = colorbar;
    c.Label.String = colorbar_names{j};
    title(colorbar_names{j});
    xlabel('Lateral y-direction [mm]');
    ylabel('Longitudinal x-direction [mm]');
    set(h(j), 'CLim', clim_values(j, :));
end

pause(1);

plot_ind = [1:10:1e4, 1e4:100:length(time_span)];

% Animation Loop
for t = plot_ind
    for j = 1:9
        set(im(j), 'CData', sim_solution(:, :, t, j));
    end
    
    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Reduce lag
    if mod(t, 10) == 0
        pause(0.001);
    end
end

% Finalize Video
close(v);
fprintf("Animation Successfully saved as Video!\n");


%%


% Plot an initial frame and store the handle to the image object
im = imagesc(y_vals, x_vals, sim_solution(:, :, 1, 1));  % Initial frame
c = colorbar;
c.Label.String = "Pressure [MPa]";
title("Pressure Distribution");
xlabel('Lateral y-direction [mm]');
ylabel('Longitudinal x-direction [mm]');

% Loop to update the image data
for t = 1:length(time_span)

    pause(0.01);  % Pause to control animation speed
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'CData', sim_solution(:, :, t, 1));  % Update the data
end

%%
function [time_span, omega, v0] = setupSimulationParams(t_final, fs, numBrushes, re, SR, rpm)
    % Sets up simulation parameters for brush dynamics
    %
    % Parameters:
    %   t_final    - Final time (seconds)
    %   fs         - Sampling frequency (Hz)
    %   numBrushes - Number of brushes in one dimension
    %   re         - Effective radius (m)
    %   SR         - Slip ratio
    %   rpm        - Maximum rotational speed (RPM)
    
    t_initial = 0;
    num_steps = t_final * fs;
    time_span = linspace(t_initial, t_final, num_steps);
    
    % Angular velocity: ramp up to max RPM for first half, then constant
    omega_max = rpm * 2 * pi / 60; % Convert RPM to rad/s
    half_steps = floor(num_steps/2);
    
    omega_profile = [linspace(0, omega_max, half_steps), ...
                    omega_max * ones(1, num_steps - half_steps)]';
    
    omega = repmat(omega_profile, 1, numBrushes^2);
    
    % Calculate v0 based on slip ratio
    v0 = omega * re ./ (SR + 1);
end