clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
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


numBrushes = 20;

re = 30;
alpha = deg2rad(0);
fs = 2.5e5; % Hz
dt = 1 / fs;
omega_z = 0;

t_initial = 0;
t_final = 2;
time_span = linspace(t_initial, t_final, t_final * fs);

omega = [linspace(0, 10, length(time_span) / 2), repmat(10, 1, length(time_span) / 2)];

SR = 0.1;
v0 = omega * re / (SR + 1);

Fz = 100;
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

x_vals = linspace(-a, a, numBrushes);
y_vals = linspace(-b, b, numBrushes);
dx = x_vals(2) - x_vals(1);
dy = y_vals(2) - y_vals(1);
dA = dx * dy;

[X, Y] = meshgrid(x_vals, y_vals);

%
[~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
P_grid = real(P_grid);

% Ensure output directory exists
outputDir = 'BrushSim_2_500_000Hz_20x20_solution';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

segmentLength = 5e4;
numSegments = length(time_span) / segmentLength;
%
for i = 1:numSegments
    % Compute start and end indices for the segment
    startIdx = (i - 1) * segmentLength + 1;
    
    % Ensure last segment includes remaining time points
    if i == numSegments
        endIdx = length(time_span);
    else
        endIdx = i * segmentLength;
    end

    % Extract time segment
    time_vec = time_span(startIdx:endIdx);
    omega_vec = omega(startIdx:endIdx);
    v0_vec = v0(startIdx:endIdx);

    sim_solution = simulateBrushModel_mex(numBrushes, P_grid(:), time_vec(:), dt, omega_vec(:), v0_vec(:), re, alpha, omega_z, X(:), Y(:));
    
    % Save results for this segment
    filename = fullfile(outputDir, sprintf('sim_solution_part%d.mat', i));
    save(filename, 'sim_solution', '-v7.3');
    fprintf('Saved results for t = %.2f to t = %.2f\n', time_vec(1), time_vec(end));

    % Clear variable to free up memory
    clear sim_solution;


end
toc(total_time)

%%
post_process = tic;
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

toc(post_process)
%%
% % 
% % fprintf('Simulation Solution is ready for animation! \n');
% % 
% % % Save simulation as video
% % close all;
% % fprintf("Starting animation ...\n")
% % 
% % % Create figure
% % figure;
% % % Pre-create the subplots and store the handles
% % h = gobjects(1, 9);  % Preallocate for subplot handles
% % 
% % % File name for the video
% % video_filename = 'BrushV2_sim_8_47_100Hz_5slip_80x80_verlet_solver_RollingPres.mp4';
% % 
% % % Create a VideoWriter object
% % v = VideoWriter(video_filename, 'MPEG-4');  % Use 'MPEG-4' for .mp4 format
% % v.FrameRate = 120;  % Adjust the frame rate as needed
% % open(v);
% % 
% % % Initialize subplots
% % for j = 1:9
% %     h(j) = subplot(3, 3, j);
% %     % Plot an initial frame and store the handle to the image object
% %     im(j) = imagesc(y_vals, x_vals, sim_solution(:, :, 1, j));  % Initial frame
% %     c = colorbar;
% %     c.Label.String = colorbar_names{j};
% %     title(colorbar_names{j});
% %     xlabel('Lateral y-direction [mm]');
% %     ylabel('Longitudinal x-direction [mm]');
% % 
% %     set(h(j), 'CLim', [0, 0.5]);
% % 
% %     if j == 2
% %         set(h(j), 'CLim', [0, 2.0]);  % Use 'set' to change the color limit
% %         % title('Displacement Magnitude [mm]')
% %     elseif j == 3
% %         set(h(j), 'CLim', [0, 1.0]);
% %         % title('Shear Stress Magnitude [MPa]')
% %     elseif j == 4
% %         set(h(j), 'CLim', [-11, 0.1]);
% %         % title('X Displacement [mm]')
% %     elseif j == 5
% %         set(h(j), 'CLim', [-0.5, 0.5]);
% %         % title('Y Displacement [mm]')
% %     elseif j == 6
% %         set(h(j), 'CLim', [-2, 0.2]);
% %         % title('Longitudinal Shear Stress [MPa]')
% %     elseif j == 7
% %         set(h(j), 'CLim', [-0.2, 0.2]);
% %         % title('Lateral Shear Stress [MPa]')
% %     elseif j == 8
% %         set(h(j), 'CLim', [0, 1]);
% %         % title('Sliding Check [boolean]')
% %     elseif j == 9
% %         set(h(j), 'CLim', [0, 1.3]);
% %         % title('Friction Coefficient')
% %     end
% % 
% % end
% % 
% % pause(1);
% % 
% % % Loop to update the image data
% % for t = 1:length(time_span)
% % 
% %     pause(0.01);  % Pause to control animation speed
% % 
% %     for j = 1:9
% %         % Update the data in the existing image plot rather than redrawing it
% %         set(im(j), 'CData', sim_solution(:, :, t, j));  % Update the data
% %     end
% % 
% %     % Capture the plot as a frame
% %     frame = getframe(gcf);
% % 
% %     % Write the frame to the video file
% %     writeVideo(v, frame);
% % end
% % 
% % % Close the video file
% % close(v);
% % 
% % fprintf("Animation Successfully saved as Video!\n")
% % %%
% % 
% % 
% % % Plot an initial frame and store the handle to the image object
% % im = imagesc(y_vals, x_vals, sim_solution(:, :, 1, 1));  % Initial frame
% % c = colorbar;
% % c.Label.String = "Pressure [MPa]";
% % title("Pressure Distribution");
% % xlabel('Lateral y-direction [mm]');
% % ylabel('Longitudinal x-direction [mm]');
% % 
% % % Loop to update the image data
% % for t = 1:length(time_span)
% % 
% %     pause(0.01);  % Pause to control animation speed
% %     % Update the data in the existing image plot rather than redrawing it
% %     set(im, 'CData', sim_solution(:, :, t, 1));  % Update the data
% % end