clear; close all; clc

run('main.m')
v0 = model_input.v0;
omega = model_input.omega;
dt_sim = model_input.dt_sim;
Fz = 560 * 9.81;

fprintf("Finished simulation in %.2fs! \n", toc(total_time))
whos sim_solution omega v0


shift_amount_cumulative = (cumsum(v0 * dt_sim));
shift_amount = (gradient(floor(shift_amount_cumulative)) > 0);

post_process = tic;
forceX              = zeros(K, model_input.LenTime_save);
forceY              = zeros(K, model_input.LenTime_save);
forceTotal          = zeros(K, model_input.LenTime_save);
forceX_medfilt      = zeros(K, model_input.LenTime_save);
forceY_medfilt      = zeros(K, model_input.LenTime_save);
forceTotal_medfilt  = zeros(K, model_input.LenTime_save);
avg_mu              = zeros(K, model_input.LenTime_save);

forceX_filt         = zeros(K, model_input.LenTime_save);
forceY_filt         = zeros(K, model_input.LenTime_save);
forceTotal_filt     = zeros(K, model_input.LenTime_save);

Y_mag2              = zeros(round(model_input.LenTime_save / 2), K);
PSD                 = zeros(round(model_input.LenTime_save / 2), K);

lgd = cell(1, K);
lgd2 = cell(1, K);


for i = 1:K
    working_data = sim_solution{i};

    P_threshold = 0.02;
    pressure_mask = max(working_data.PressGrid, P_threshold) > P_threshold;
    total_valid_points = max(sum(pressure_mask( :, :), 1)); % Count valid points

    % Intergrate stress to get force
    forceX(i, :) = squeeze(sum(working_data.tauX)) * model_input.dA;
    forceX_medfilt(i, :) = medfilt1(forceX(i, :));
    
    forceY(i, :) = squeeze(trapz(working_data.tauY)) * model_input.dA;
    forceY_medfilt(i, :) = medfilt1(forceY(i, :));
    
    % % %%%%%%%%%%%%%%%%%% Save forces for regression purposes %%%%%%%%%%%%%%%%%%%%
    % % X_data = cat(2, forceX(:), forceY(:));
    % % save('BrushModel_data.mat', 'v0', 'SR', 't_save', 'X_data', '-v7.3')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    forceTotal(i, :) = squeeze(sum(working_data.TotalStress)) * model_input.dA;
    % forceTotal = dot(working_data( :, :, 7), P_grid.save, 1) * dA;
    forceTotal_medfilt(i, :) = medfilt1(forceTotal(i, :));

    avg_mu(i, :) = squeeze(sum(working_data.mu .* pressure_mask)./ total_valid_points); % Avoid division by zero
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


% % Moving Average
% winSize = 10;
% b_val = (1 / winSize) * ones(1, winSize);
% a_val = 1;
% 
% for i = 1:K
%     forceTotal_filt(i, :) = filtfilt(b_val, a_val, double(forceTotal_medfilt(i, :)));
%     forceX_filt(i, :) = filtfilt(b_val, a_val, double(forceX_medfilt(i, :)));
%     forceY_filt(i, :) = filtfilt(b_val, a_val, double(forceY_medfilt(i, :)));
% end

dt_ratio = uint32(model_input.dt_ratio);

if ~isRolling
    vel = max(abs(v0));
else
    vel = max(abs(model_input.re * omega));
end

for i = 1:K
    lgd{i} = sprintf("Linear Velocity of wheel v_{max} = %.1f m/s", max(model_input.re * omega(:, i), [], "all"));
    lgd{i+1} = sprintf("Velocity of Road element v_{max} = %.1f m/s", vel(i));
    lgd2{i} = sprintf("v_{max} = %.1f m/s", vel(i)); 
end

T = tiledlayout('horizontal');
T.Padding = "tight";
T.TileSpacing = "tight";

nexttile
hold on
plot(t_save, omega(1:dt_ratio:end) * model_input.re)
plot(t_save, v0(1:dt_ratio:end))
hold off
grid on
title('Input Velocities')
legend(lgd)
xlabel('Time[s]');ylabel('Velocity [m/s]')

nexttile
plot(t_save, forceTotal)
grid on
title('Total Force')
xlabel('Time[s]');ylabel('Force [N]')
ylim([0, 2 * Fz])
legend(lgd2)

for i = 1:K
    [Y_mag2(:, i), half_freqz] = One_sided_fft(forceTotal(i, :), fs_sim);
    N = 2 * length(Y_mag2(:, i));

    PSD(:, i) = (1 /(fs_sim * N)) * real(Y_mag2(:, i)  .* conj(Y_mag2(:, i) ));
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
T = tiledlayout('horizontal');
T.Padding = "tight";
T.TileSpacing = "tight";

nexttile
plot(t_save, (forceX));
hold on
grid on;
title("Longitudinal Force [N]")
legend(lgd2)
hold off
xlabel('Time [s]');ylabel('Force [N]')
ylim([0, 2 * Fz])

nexttile
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

nexttile
plot(t_save, forceTotal);
hold on
plot(t_save, -(avg_mu * Fz), '--');
title('Total Force')
legend(lgd2)
grid on
ylim([0, 2 * Fz])
hold off
xlabel('Time [s]');ylabel('Force [N]')

%%%%%%%%%%%%%%%%%% Plot Force vs Slip %%%%%%%%%%%%%%%%%%%%%%%
% Identify time range where slip ranges from 0% to 100%
ind = (t_save >= 11) .* (t_save <= 101);
% Remove all the indices where ind is less than 1
ind = ind > 0;

% Plot force vs slip using indices calculated
if isRolling
    SR_for_plot = model_input.SR(1:dt_ratio:end);
    figure
    plot(SR_for_plot(ind), forceX(ind))
    hold on
    plot(SR_for_plot(ind), avg_mu(ind) .* Fz)
    grid on
    xlabel('Longitudinal Slip')
    ylabel('Force Simulated [N]')
    title('Force vs Slip graph generated from Brush Model')
    legend('Simulated Force [N]', '\mu_{avg} \times F_z [N]', Location='best')
else
    v0_for_plot = model_input.v0(1:dt_ratio:end);
    figure
    plot(v0_for_plot(ind), abs(forceX(ind)))
    hold on
    plot(v0_for_plot(ind), avg_mu(ind) .* Fz)
    grid on
    xlabel('Longitudinal Slip Velocity [mm/s]')
    ylabel('Force Simulated [N]')
    title('Force vs Slip graph generated from Brush Model')
    legend('Simulated Force [N]', '\mu_{avg} \times F_z [N]', Location='best')
end
toc(post_process)


plot_ind = 1:10:model_input.LenTime_save;
figure
T = tiledlayout('horizontal');
T.Padding = "tight";
T.TileSpacing = "tight";
for i = 1:K
    working_data = sim_solution{i};
    
    disp = squeeze(working_data.delta_x(1:200, plot_ind)).';
    force = squeeze(working_data.tauX(1:200, plot_ind)).';
    slide_vel = squeeze(working_data.vs(1:200, plot_ind)).';
    
    % subplot(2,3,i)
    nexttile
    plot(disp, force, '.')
    grid on
    % ylim([-0.2, 0.2])
    title(sprintf('Phase Plot v_{max} = %.2f', max( abs( v0(:, i) ) ) ) )
    xlabel('Rel. Disp. [mm]')
    ylabel('Long. Stress [MPa]')
    
    
    % subplot(2,3,i+3)
    nexttile
    plot(slide_vel, force, '.')
    grid on
    % ylim([-0.2, 0.2])
    title( sprintf('Phase Plot v_{max} = %.2f', max( abs( v0(:, i) ) ) ) )
    xlabel('Sliding Vel. [mm/s]')
    ylabel('Long. Stress [MPa]')
end

% % figure(7)
% % comet(disp(:, 10), force(:, 10))
% % grid on
% % ylim([-0.2, 0.2])

%

plot_solution = struct();

working_data = sim_solution{1};
savedFields = fieldnames(working_data);

if size(plot_solution, 1) ~= numBrushes
    for i = 1:length(savedFields)
        plot_solution.(savedFields{i}) = reshape(working_data.(savedFields{i}), numBrushes, numBrushes, model_input.LenTime_save);
    end
end

fprintf('Simulation Solution is ready for animation! \n');

% Initialize Figure
if ishandle(7)
    close(7);
    % if exist(v)
    %     close(v);
    % end
end

% Initialize Video
% video_filename = sprintf('NoisyPress_10s_%dN_BrushSim_%dHz_%drpm_slip%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, rpm, SR, omega_z, alpha);
% % video_filename = sprintf('TM700_Sliding_100s_%dN_BrushSim_%dHz_max_omega%.2f_omegaZ%.2f_alpha%.2f.mp4', ...
% %                          Fz, model_input.fs_sim, max(abs(v0(:))), model_input.omega_z, model_input.alpha);
% % Video_path = fullfile(strcat("C:\Users\coxde\OneDrive\Masters\BrushV2\Animations/",video_filename));
% % v = VideoWriter(Video_path, 'MPEG-4');
% % v.FrameRate = 60;  
% % open(v);

figure(7);
h = gobjects(1, 9);
im = gobjects(1, 9);

clim_values = [
    0, 1.0;
    0, 0.15;
    0, 0.8;
    -0.15, 0.15;
    -0.15, 0.15;
    -0.8, 0.8;
    -0.1, 0.1;
    0, 2.5;
    -10, 10;
];

% Initialize subplots
for j = 1:length(savedFields)
    h(j) = subplot(3, 3, j);
    im(j) = imagesc(model_input.X(1, :), model_input.Y(:, 1), plot_solution.(savedFields{j})(:, :, 1));
    c = colorbar;
    c.Label.String = savedFields{j};
    title(savedFields{j});
    ylabel('Lateral y-direction [mm]');
    xlabel('Longitudinal x-direction [mm]');
    set(h(j), 'CLim', clim_values(j, :));
end

pause(1);

plot_ind = 1:50:model_input.LenTime_save;

% Animation Loop
for t = plot_ind
    for j = 1:9
        set(im(j), 'CData', plot_solution.(savedFields{j})(:, :, t));
    end

    % Capture frame
    frame = getframe(gcf);
    % % writeVideo(v, frame);

    % Reduce lag
    if mod(t, 10) == 0
        pause(0.001);
    end
end

% Finalize Video
% close(v);
fprintf("Animation Successfully saved as Video!\n");

%%
% % plot_solution = 0;
% % 
% % working_data = sim_solution{3};
% % 
% % if size(plot_solution, 1) ~= numBrushes 
% %     plot_solution = cat(3, working_data(:, :, 4), working_data(:, :, 6));
% %     plot_solution = reshape(plot_solution, numBrushes, numBrushes, LenTime_save, 2);
% % end
% % 
% % fprintf('Simulation Solution is ready for animation! \n');
% % 
% % % Initialize Figure
% % if ishandle(8)
% %     close(8);
% %     close(v);
% % end
% % 
% % % Initialize Video
% % video_filename = sprintf('ISTVS_Sliding_100s_SurfacePlot_%dN_BrushSim_%dHz_max_omega_%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, max(abs(v0(:, 3))), omega_z, alpha);
% % Video_path = fullfile(strcat("Animations/",video_filename));
% % v = VideoWriter(Video_path, 'MPEG-4');
% % v.FrameRate = 60;  
% % open(v);
% % 
% % figure(8);
% % h = gobjects(1, 2);
% % im = gobjects(1, 2);
% % 
% % % Predefined Color Limits
% % clim_values = [
% %     -0.01, 0.05;
% %     -0.6, 0.6;
% % ];
% % 
% % names = cell(1,2);
% % names{1} = 'Relative Displacement [mm]';
% % names{2} = 'Longitudinal Shear Stress [MPa]';
% % 
% % % Initialize subplots
% % for j = 1:2
% %     h(j) = subplot(1, 2, j);
% %     im(j) = surf(X, Y, plot_solution(:, :, 1, j));
% %     shading interp
% %     colormap jet
% %     c = colorbar;
% %     c.Label.String = names{j};
% %     title(names{j});
% %     ylabel('Lateral y-direction [mm]');
% %     xlabel('Longitudinal x-direction [mm]');
% %     set(h(j), 'CLim', clim_values(j, :));
% %     set(h(j), 'Zlim', clim_values(j, :));
% % end
% % 
% % pause(1);
% % 
% % plot_ind = 1:100:LenTime_save;
% % 
% % % Animation Loop
% % for t = plot_ind
% %     for j = 1:2
% %         set(im(j), 'ZData', plot_solution(:, :, t, j));
% %     end
% % 
% %     % Capture frame
% %     frame = getframe(gcf);
% %     writeVideo(v, frame);
% % 
% %     % Reduce lag
% %     if mod(t, 10) == 0
% %         pause(0.001);
% %     end
% % end
% % 
% % % Finalize Video
% % close(v);
% % fprintf("Animation Successfully saved as Video!\n");
% % 
% % %%
% % if size(sim_solution, 1) ~= numBrushes 
% %     pressure_mask = P_grid_save > 1e-3;
% %     P_grid_save = reshape(P_grid_save .* pressure_mask, numBrushes, numBrushes, LenTime_save);
% % end
% % 
% % % Initialize Video
% % video_filename = 'Shifted_Pressure';
% % Video_path = fullfile(strcat("Animations/",video_filename));
% % v = VideoWriter(Video_path, 'MPEG-4');
% % v.FrameRate = 60;  
% % open(v);
% % 
% % % Plot an initial frame and store the handle to the image object
% % im = imagesc(y_vals, x_vals, P_grid_save(:, :, 1));  % Initial frame
% % c = colorbar;
% % c.Label.String = "Pressure [MPa]";
% % title("Pressure Distribution");
% % ylabel('Lateral y-direction [mm]');
% % xlabel('Longitudinal x-direction [mm]');
% % 
% % % Loop to update the image data
% % for t = 1:50:length(time_span)
% %     % Update the data in the existing image plot rather than redrawing it
% %     set(im, 'CData', P_grid_save(:, :, t));  % Update the data
% % 
% %     % Capture frame
% %     frame = getframe(gcf);
% %     writeVideo(v, frame);
% % 
% %     % Reduce lag
% %     if mod(t, 10) == 0
% %         pause(0.001);
% %     end
% % end
% % 
% % close(v);

%%
% % fprintf('Viewing Boolean array! \n');
% % 
% % if size(bool_array, 1) ~= numBrushes 
% %     bool_array = reshape(bool_array, numBrushes, numBrushes, LenTime_save, 2);
% %     colorbar_names{1} = "Sliding Condition";
% %     colorbar_names{2} = "Passing Condition";
% % end
% % 
% % % Initialize Figure
% % if ishandle(4)
% %     close(4);
% % end
% % 
% % fh = figure(4);
% % % fh.WindowState = 'maximized';
% % h = gobjects(1, 2);
% % im = gobjects(1, 2);
% % 
% % % Predefined Color Limits
% % clim_values = [
% %     0, 1;
% %     0, 1;
% % ];
% % 
% % % Initialize subplots
% % for j = 1:2
% %     h(j) = subplot(2, 1, j);
% %     im(j) = imagesc(y_vals, x_vals, bool_array(:, :, 1, j));
% %     c = colorbar;
% %     c.Label.String = colorbar_names{j};
% %     title(colorbar_names{j});
% %     ylabel('Lateral y-direction [mm]');
% %     xlabel('Longitudinal x-direction [mm]');
% %     set(h(j), 'CLim', clim_values(j, :));
% % end
% % 
% % pause(1);
% % 
% % plot_ind = 1:1:LenTime_save;
% % 
% % % Animation Loop
% % for t = plot_ind
% %     for j = 1:2
% %         set(im(j), 'CData', bool_array(:, :, t, j));
% %     end
% %     if ~any(bool_array(:, :, t, 1))
% %         pause(10);
% %     else
% %         pause(0.001);
% %     end
% %     % end
% % end
% % 
% % %%
% % 
% % if size(sim_solution, 1) ~= numBrushes 
% %     X_grid_save = reshape(squeeze(sim_solution(:, :, 8)), numBrushes, numBrushes, LenTime_save);
% % end
% % 
% % % Plot an initial frame and store the handle to the image object
% % im = imagesc(y_vals, x_vals, X_grid_save(:, :, 1));  % Initial frame
% % c = colorbar;
% % c.Label.String = "X values [mm]";
% % title("Distribution of X_values over time");
% % ylabel('Lateral y-direction [mm]');
% % xlabel('Longitudinal x-direction [mm]');
% % 
% % 
% % plot_ind = 1:1:LenTime_save;
% % 
% % % Loop to update the image data
% % for t = plot_ind
% %     % Update the data in the existing image plot rather than redrawing it
% %     set(im, 'CData', X_grid_save(:, :, t));  % Update the data
% % 
% %     % Reduce lag
% %     if mod(t, 10) == 0
% %         pause(1);
% %     end
% % end
