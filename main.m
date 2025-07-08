clear; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%

clearvars; clc;
close all;



total_time = tic;
pressure_start = tic; 
fprintf("Starting simulation. \n")
fs_save = 1e3; % Hz
fs_sim = 1e3; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 10;
isRolling = true;

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

fprintf("Initialised brush model simulation in %.2fs! \n", toc(pressure_start))

%%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%
% % Press = reshape(P_grid.shifted(:, 1), numElems, numElems);
sliding = false;
K = min(size(model_input.omega), [], 'all');
% for i = 1:K
%     sim_solution{i} = simulateBrushModel_V2(model_input);
% end
sim_solution = simulateBrushModel_V2_mex(model_input);

%%


%%%%%%%%%% Sliding tyre %%%%%%%%%%%%%%%%%%%
%   P_grid must be constant, since the brushes will not move through contact patch
%   Supply v0 instead of omega
% Press = P_grid_subsampled; %reshape(P_grid.shifted(:, 1), numElems, numElems);
% K = min(size(v0), [], 'all');
% for i = 1:K
%     sim_solution{i} = simulateBrushModel_V2(numElems, LenTime_sim, LenTime_save, ...
%                                                         Press, dt_sim, dt_save, v0(:, i), ...
%                                                         re, alpha, omega_z, X, Y);
% end

fprintf("Finished simulation in %.2fs! \n", toc(total_time))
whos sim_solution omega v0
 
shift_amount_cumulative = (cumsum(v0 * dt_sim));
shift_amount = (gradient(floor(shift_amount_cumulative)) > 0);

post_process = tic;
for i = 1:K
    working_data = sim_solution{i};

    P_threshold = 0.02;
    pressure_mask = max(working_data(:, :, 1), P_threshold) > P_threshold;
    total_valid_points = max(sum(pressure_mask( :, :), 1)); % Count valid points

    % Intergrate stress to get force
    forceX(i, :) = squeeze(trapz(working_data( :, :, 6))) * dA;
    forceX_medfilt(i, :) = medfilt1(forceX(i, :));
    
    forceY(i, :) = squeeze(trapz(working_data( :, :, 5))) * dA;
    forceY_medfilt(i, :) = medfilt1(forceY(i, :));
    
    % % %%%%%%%%%%%%%%%%%% Save forces for regression purposes %%%%%%%%%%%%%%%%%%%%
    % % X_data = cat(2, forceX(:), forceY(:));
    % % save('BrushModel_data.mat', 'v0', 'SR', 't_save', 'X_data', '-v7.3')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    forceTotal(i, :) = squeeze(trapz(working_data( :, :, 3))) * dA;
    % forceTotal = dot(working_data( :, :, 7), P_grid.save, 1) * dA;
    forceTotal_medfilt(i, :) = medfilt1(forceTotal(i, :));

    avg_mu(i, :) = squeeze(sum(working_data( :, :, 8) .* pressure_mask)./ total_valid_points); % Avoid division by zero
end


colorbar_names = cell(9,1);
colorbar_names{1} = "Pressure Distribution [MPa]";
colorbar_names{2} = "Total Displacement [mm]";
colorbar_names{3} = "Total Shear Stress [MPa]";
colorbar_names{4} = "Relative X Displacement [mm]";
colorbar_names{5} = "Relative Y Displacement [mm]";
colorbar_names{6} = "\tau_x (X Shear Stress) [MPa]";
colorbar_names{7} = "\tau_y (Y Shear Stress) [MPa]";
% colorbar_names{9} = "Sliding Condition [Boolean]";
colorbar_names{8} = "Friction Coefficient";
colorbar_names{9} = 'Sliding Velocity [mm/s]';


% Moving Average
winSize = 10;
b_val = (1 / winSize) * ones(1, winSize);
a_val = 1;

for i = 1:K
    forceTotal_filt(i, :) = filtfilt(b_val, a_val, double(forceTotal_medfilt(i, :)));
    forceX_filt(i, :) = filtfilt(b_val, a_val, double(forceX_medfilt(i, :)));
    forceY_filt(i, :) = filtfilt(b_val, a_val, double(forceY_medfilt(i, :)));
end

dt_ratio = uint32(dt_ratio);

if sliding
    vel = max(abs(v0));
else
    vel = max(abs(re * omega));
end

for i = 1:K
    % lgd{i} = sprintf("Linear Velocity of wheel v_{max} = %.1f m/s", max(re * omega(:, i), [], "all"));
    lgd{i} = sprintf("Velocity of Road element v_{max} = %.1f m/s", vel(i));
    lgd2{i} = sprintf("v_{max} = %.1f m/s", vel(i)); 
end

figure
% plot(t_save, 40 * shift_amount(1:dt_ratio:end))
subplot(121)
hold on
% plot(t_save, omega * re)
plot(t_save, v0)
grid on
title('Input Velocities')
legend(lgd)
xlabel('Time[s]');ylabel('Velocity [m/s]')
% ylim([-2 * Fz, 2 * Fz])
hold off

subplot(122)
plot(t_save, forceTotal)
grid on
title('Total Force')
xlabel('Time[s]');ylabel('Force [N]')
ylim([0, 2 * Fz])
legend(lgd2)

for i = 1:K
    [Y_mag2(i, :), half_freqz] = One_sided_fft(forceTotal(i, :), fs_sim);
    N = 2 * length(Y_mag2(i, :));

    PSD(i, :) = (1 /(fs_sim * N)) * real(Y_mag2(i, :)  .* conj(Y_mag2(i, :) ));
end
figure
loglog(half_freqz, Y_mag2)

grid on
title('PSD of Force and Velocity')
xlabel('Frequency [Hz]');ylabel('Magnitude [N^2/Hz]')
legend(lgd2)

for i = 1:K
    lgd2{i} = sprintf("Simulated: v_{max} = %.1f m/s", vel(i)); 
end

figure
subplot(131)
plot(t_save, (forceX));
hold on
grid on;
title("Longitudinal Force [N]")
legend(lgd2)
hold off
xlabel('Time [s]');ylabel('Force [N]')
ylim([0, 2 * Fz])

subplot(132)
plot(t_save, forceY);
grid on;
title("Lateral Force [N]")
legend(lgd2)
hold off
xlabel('Time [s]');ylabel('Force [N]')

for i = 1:K
    lgd2{i} = sprintf("Total Force: v_{max} = %.1f m/s",vel(i) ); 
    lgd2{i+K} = sprintf("Avg mu X F_z: v_{max} = %.1f m/s", vel(i) );
end

subplot(133)
plot(t_save, forceTotal);
hold on
plot(t_save, (avg_mu * Fz), '--');
title('Total Force')
legend(lgd2)
grid on
ylim([0, 2 * Fz])
hold off
xlabel('Time [s]');ylabel('Force [N]')

toc(post_process)

%%
plot_ind = 1:10:LenTime_save;
figure
for i = 1:3
    working_data = sim_solution{i};
    
    disp = squeeze(working_data(1:200, plot_ind, 4)).';
    force = squeeze(working_data(1:200, plot_ind, 6)).';
    slide_vel = squeeze(working_data(1:200, plot_ind, 9)).';
    
    subplot(2,3,i)
    plot(disp, force, '.')
    grid on
    ylim([-0.2, 0.2])
    title(sprintf('Phase Plot v_{max} = %.2f'),max(abs(v0(:, i))))
    xlabel('Rel. Disp. [mm]')
    ylabel('Long. Stress [MPa]')
    
    
    subplot(2,3,i+3)
    plot(slide_vel, force, '.')
    grid on
    ylim([-0.2, 0.2])
    title(sprintf('Phase Plot v_{max} = %.2f'), max(abs(v0(:, i))))
    xlabel('Sliding Vel. [mm/s]')
    ylabel('Long. Stress [MPa]')
end

% % figure(7)
% % comet(disp(:, 10), force(:, 10))
% % grid on
% % ylim([-0.2, 0.2])

%%

plot_solution = 0;

working_data = sim_solution{3};

if size(plot_solution, 1) ~= numBrushes 
    plot_solution = reshape(working_data, numBrushes, numBrushes, LenTime_save, 9);
end

fprintf('Simulation Solution is ready for animation! \n');

% Initialize Figure
if ishandle(7)
    close(7);
    if exist(v)
        close(v);
    end
end

% Initialize Video
% video_filename = sprintf('NoisyPress_10s_%dN_BrushSim_%dHz_%drpm_slip%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, rpm, SR, omega_z, alpha);
video_filename = sprintf('ISTVS_TM700_Sliding_100s_%dN_BrushSim_%dHz_max_omega%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, max(abs(v0(:, 3))), omega_z, alpha);
Video_path = fullfile(strcat("Animations/",video_filename));
v = VideoWriter(Video_path, 'MPEG-4');
v.FrameRate = 60;  
open(v);

fh = figure(7);
h = gobjects(1, 9);
im = gobjects(1, 9);

clim_values = [
    0, 0.52;
    0, 0.05;
    0, 0.6;
    -0.05, 0.05;
    -0.05, 0.05;
    -0.6, 0.6;
    -0.1, 0.1;
    0, 2.5;
    -10, 10;
];

% Initialize subplots
for j = 1:9
    h(j) = subplot(3, 3, j);
    im(j) = imagesc(x_vals, y_vals, plot_solution(:, :, 1, j));
    c = colorbar;
    c.Label.String = colorbar_names{j};
    title(colorbar_names{j});
    ylabel('Lateral y-direction [mm]');
    xlabel('Longitudinal x-direction [mm]');
    set(h(j), 'CLim', clim_values(j, :));
end

pause(1);

plot_ind = 1:100:LenTime_save;

% Animation Loop
for t = plot_ind
    for j = 1:9
        set(im(j), 'CData', plot_solution(:, :, t, j));
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
plot_solution = 0;

working_data = sim_solution{3};

if size(plot_solution, 1) ~= numBrushes 
    plot_solution = cat(3, working_data(:, :, 4), working_data(:, :, 6));
    plot_solution = reshape(plot_solution, numBrushes, numBrushes, LenTime_save, 2);
end

fprintf('Simulation Solution is ready for animation! \n');

% Initialize Figure
if ishandle(8)
    close(8);
    close(v);
end

% Initialize Video
video_filename = sprintf('ISTVS_Sliding_100s_SurfacePlot_%dN_BrushSim_%dHz_max_omega_%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, max(abs(v0(:, 3))), omega_z, alpha);
Video_path = fullfile(strcat("Animations/",video_filename));
v = VideoWriter(Video_path, 'MPEG-4');
v.FrameRate = 60;  
open(v);

fh = figure(8);
h = gobjects(1, 2);
im = gobjects(1, 2);

% Predefined Color Limits
clim_values = [
    -0.01, 0.05;
    -0.6, 0.6;
];

names = cell(1,2);
names{1} = 'Relative Displacement [mm]';
names{2} = 'Longitudinal Shear Stress [MPa]';

% Initialize subplots
for j = 1:2
    h(j) = subplot(1, 2, j);
    im(j) = surf(X, Y, plot_solution(:, :, 1, j));
    shading interp
    colormap jet
    c = colorbar;
    c.Label.String = names{j};
    title(names{j});
    ylabel('Lateral y-direction [mm]');
    xlabel('Longitudinal x-direction [mm]');
    set(h(j), 'CLim', clim_values(j, :));
    set(h(j), 'Zlim', clim_values(j, :));
end

pause(1);

plot_ind = 1:100:LenTime_save;

% Animation Loop
for t = plot_ind
    for j = 1:2
        set(im(j), 'ZData', plot_solution(:, :, t, j));
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
if size(sim_solution, 1) ~= numBrushes 
    pressure_mask = P_grid_save > 1e-3;
    P_grid_save = reshape(P_grid_save .* pressure_mask, numBrushes, numBrushes, LenTime_save);
end

% Initialize Video
video_filename = 'Shifted_Pressure';
Video_path = fullfile(strcat("Animations/",video_filename));
v = VideoWriter(Video_path, 'MPEG-4');
v.FrameRate = 60;  
open(v);

% Plot an initial frame and store the handle to the image object
im = imagesc(y_vals, x_vals, P_grid_save(:, :, 1));  % Initial frame
c = colorbar;
c.Label.String = "Pressure [MPa]";
title("Pressure Distribution");
ylabel('Lateral y-direction [mm]');
xlabel('Longitudinal x-direction [mm]');

% Loop to update the image data
for t = 1:50:length(time_span)
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'CData', P_grid_save(:, :, t));  % Update the data

    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Reduce lag
    if mod(t, 10) == 0
        pause(0.001);
    end
end

close(v);

%%
fprintf('Viewing Boolean array! \n');

if size(bool_array, 1) ~= numBrushes 
    bool_array = reshape(bool_array, numBrushes, numBrushes, LenTime_save, 2);
    colorbar_names{1} = "Sliding Condition";
    colorbar_names{2} = "Passing Condition";
end

% Initialize Figure
if ishandle(4)
    close(4);
end

fh = figure(4);
% fh.WindowState = 'maximized';
h = gobjects(1, 2);
im = gobjects(1, 2);

% Predefined Color Limits
clim_values = [
    0, 1;
    0, 1;
];

% Initialize subplots
for j = 1:2
    h(j) = subplot(2, 1, j);
    im(j) = imagesc(y_vals, x_vals, bool_array(:, :, 1, j));
    c = colorbar;
    c.Label.String = colorbar_names{j};
    title(colorbar_names{j});
    ylabel('Lateral y-direction [mm]');
    xlabel('Longitudinal x-direction [mm]');
    set(h(j), 'CLim', clim_values(j, :));
end

pause(1);

plot_ind = 1:1:LenTime_save;

% Animation Loop
for t = plot_ind
    for j = 1:2
        set(im(j), 'CData', bool_array(:, :, t, j));
    end
    if ~any(bool_array(:, :, t, 1))
        pause(10);
    else
        pause(0.001);
    end
    % end
end

%%

if size(sim_solution, 1) ~= numBrushes 
    X_grid_save = reshape(squeeze(sim_solution(:, :, 8)), numBrushes, numBrushes, LenTime_save);
end

% Plot an initial frame and store the handle to the image object
im = imagesc(y_vals, x_vals, X_grid_save(:, :, 1));  % Initial frame
c = colorbar;
c.Label.String = "X values [mm]";
title("Distribution of X_values over time");
ylabel('Lateral y-direction [mm]');
xlabel('Longitudinal x-direction [mm]');


plot_ind = 1:1:LenTime_save;

% Loop to update the image data
for t = plot_ind
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'CData', X_grid_save(:, :, t));  % Update the data

    % Reduce lag
    if mod(t, 10) == 0
        pause(1);
    end
end
