clearvars; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%

clearvars; 
close all; clc;

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

% a = 212.5000 / 2; % mm
% b = 234.3750 / 2; % mm

contact_area = a * b;

% n = [3, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000].';
numBrushes = 85;
contact_area_per_brush = contact_area / (numBrushes * numBrushes);

% Parameters for displacement computation

% omega = 1;
% v0 = 0;

re = 30;
alpha = deg2rad(0);
fs = 4.71e4; % Hz
dt = 1 / fs;
omega_z = 0;

t_initial = 0;
t_final = 1;
time_span = linspace(t_initial, t_final, t_final * fs);

omega = linspace(0, t_final, length(time_span));

% SR = (omega * re - v0) / v0;
SR = 0.05;
v0 = omega * re / (SR + 1); %zeros(size(omega));

% Initialize time and force arrays
time = zeros(length(time_span), 1);
% force = zeros(length(time_span), 1);

x_vals = linspace(-a, a, numBrushes);
y_vals = linspace(-b, b, numBrushes);

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

sim_solution = zeros(numBrushes, numBrushes, length(time_span), 9);
% sim_solution(:, :, :, 1) = repmat(P_grid, [1, 1, length(time_span)]);

% Calculate the amount to shift per timestep
shift_amount = (omega * dt * re);  % Adjust re for scale if needed
shift_amount_cumulative = cumsum(shift_amount);


figure
hold on
% plot(time_span, shift_amount)
plot(time_span, shift_amount_cumulative)

shift_amount_cumulative = floor(shift_amount_cumulative);

plot(time_span, shift_amount_cumulative)
grid on
legend('Continuous', 'Discrete')
ylabel('Units Shifted')
xlabel('Time [s]')

for i = 1:length(time_span)
    % Shift the pressure distribution downwards and wrap it
    sim_solution(:, :, i, 1) = circshift(P_grid, [-shift_amount_cumulative(i), 0]);
end

% % % sim_solution(:, :, :, 1) = arrayfun(@(i) circshift(P_grid, [-shift_amount_cumulative(i), 0]), 1:length(time_span), 'UniformOutput', false);
% % % sim_solution(:, :, :, 1) = cat(3, sim_solution{:, :, :, 1});

    
fprintf('Pressure Grid Successfully calculated in %.6f \n', toc(pressure_start));
clear pressure_start;

% Pre-allocate array of empty objects
vs = zeros(1, numBrushes * numBrushes);
brushArray(1,1) = BrushVec(X(:), Y(:), P_grid(:), vs(:));
%
% % for i = 1:numBrushes
% %     for j = 1:numBrushes
% %         brushArray(i, j).x = x_vals(i); 
% %         brushArray(i, j).y = y_vals(j);
% %         brushArray(i, j).press = sim_solution(i, j, 1, 1);
% %         brushArray(i, j).ky = brushArray(i, j).ky_0 + brushArray(i, j).ky_l * brushArray(i, j).press;
% %         brushArray(i, j).kx = brushArray(i, j).phi + brushArray(i, j).ky;
% %         brushArray(i, j).vs = 0;
% %     end
% % end

% [brushArray.x] = (); 
% [brushArray.y] = ();
% [brushArray.press] = (P_grid(:));
% [brushArray.ky] = deal([brushArray.ky_0] + [brushArray.ky_l] .* [brushArray.press]);
% [brushArray.kx] = deal([brushArray.phi] + [brushArray.ky]);
% [brushArray.vs] = deal(0);

%%
start_brush_sim = tic;
fprintf('Starting brush model simulation! \n')

% Initialize each object with specific values
for i = 1:length(time_span)
    for j = 1:numBrushes
        for k = 1:numBrushes
            brush = brushArray(j, j);
            brush.press = sim_solution(j, k, i, 1);
            brush.ky = brush.ky_0 + brush.ky_l * brush.press;
            brushkx = brush.phi + brush.ky;
            brush = brush.update_brush(omega(i), omega_z, re, v0(i), alpha, dt);

            % Store results from the brush object
            % sim_solution(j, k, i, 2) = "Total Displacement [mm]";
            % sim_solution(j, k, i, 3) = "Total Shear Stress [MPa]";
            % if isnan(brushArray(j, k).vs)
            %     fprintf('For brush located at (%.4f, %.4f), vs was calculated to be NAN at %.4fs. \n', x_vals(i), y_vals(j), time_span(i));
            % end
            if ~isnan(brushArray(j, k).delta_x)
                sim_solution(j, k, i, 4) = brushArray(j, k).delta_x;
            else
                brushArray(j, k).delta_x = 0;
                sim_solution(j, k, i, 4) = brushArray(j, k).delta_x;
            end
            if ~isnan(brushArray(j, k).delta_y)
                sim_solution(j, k, i, 5) = brushArray(j, k).delta_y;
            else
                brushArray(j, k).delta_y = 0;
                sim_solution(j, k, i, 5) = brushArray(j, k).delta_y;
            end
           if ~isnan(brushArray(j, k).tauX)
                sim_solution(j, k, i, 6) = brushArray(j, k).tauX;
            else
                brushArray(j, k).tauX = 0;
                sim_solution(j, k, i, 6) = brushArray(j, k).tauX;
            end
           if ~isnan(brushArray(j, k).tauY)
                sim_solution(j, k, i, 7) = brushArray(j, k).tauY;
            else
                brushArray(j, k).tauY = 0;
                sim_solution(j, k, i, 7) = brushArray(j, k).tauY;
            end

            sim_solution(j, k, i, 8) = brushArray(j, k).slide;

            if ~isnan(brushArray(j, k).mu)
                sim_solution(j, k, i, 9) = brushArray(j, k).mu;
            else
                brushArray(j, k).mu = 0;
                sim_solution(j, k, i, 9) = brushArray(j, k).mu;
            end

            if brushArray(j, k).x < -a
                brushArray(j, k).x = a;
            end
        end
    end
    % Compute total displacement over grid
    sim_solution(:, :, i, 2) = hypot(sim_solution(:, :, i, 4), sim_solution(:, :, i, 5));
    % Compute total Stress over grid
    sim_solution(:, :, i, 3) = hypot(sim_solution(:, :, i, 6), sim_solution(:, :, i, 7));
    % % % % Intergrate stress to get force
    % % % force(i) = trapz(trapz(sim_solution(:, :, i, 6))) * contact_area_per_brush;
    if mod(i, 10000) == 0
        fprintf("%d iterations completed! %d remaining.\n", i, length(time_span) - i);
    end
end

% Intergrate stress to get force
force = squeeze(trapz(trapz(sim_solution(:, :, :, 6), 1), 2)) * contact_area_per_brush;

fprintf('Brush Sim ended successfully! in %.6f \n', toc(start_brush_sim))
fprintf('Total Simulation time: %.2fs \n', toc(total_time));


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


plot(time_span, (force));
grid on;
title('Friction Force Calculated from shear stress');
xlabel('Time[s]');ylabel('Force [N]')

%%

fprintf('Simulation Solution is ready for animation! \n');

% Save simulation as video
close all;
fprintf("Starting animation ...\n")

% Create figure
figure;
% Pre-create the subplots and store the handles
h = gobjects(1, 9);  % Preallocate for subplot handles

% File name for the video
video_filename = 'BrushV2_sim_5_75_000Hz_5slip_10x10_Verlet_solver_RollingPres.mp4';

% Create a VideoWriter object
v = VideoWriter(video_filename, 'MPEG-4');  % Use 'MPEG-4' for .mp4 format
v.FrameRate = 120;  % Adjust the frame rate as needed
open(v);

% Initialize subplots
for j = 1:9
    h(j) = subplot(3, 3, j);
    % Plot an initial frame and store the handle to the image object
    im(j) = imagesc(y_vals, x_vals, sim_solution(:, :, 1, j));  % Initial frame
    c = colorbar;
    c.Label.String = colorbar_names{j};
    title(colorbar_names{j});
    xlabel('Lateral y-direction [mm]');
    ylabel('Longitudinal x-direction [mm]');
    
    % % if j == 2
    % %     set(h(j), 'CLim', [0, 3]);  % Use 'set' to change the color limit
    % %     % title('Displacement Magnitude [mm]')
    % % elseif j == 3
    % %     set(h(j), 'CLim', [0, 0.6]);
    % %     % title('Shear Stress Magnitude [MPa]')
    % % elseif j == 4
    % %     set(h(j), 'CLim', [0, 3]);
    % %     % title('X Displacement [mm]')
    % % elseif j == 5
    % %     set(h(j), 'CLim', [0, 0.1]);
    % %     % title('Y Displacement [mm]')
    % % elseif j == 6
    % %     set(h(j), 'CLim', [0, 0.6]);
    % %     % title('Longitudinal Shear Stress [MPa]')
    % % elseif j == 7
    % %     set(h(j), 'CLim', [0, 1]);
    % %     % title('Lateral Shear Stress [MPa]')
    % % elseif j == 8
    % %     % set(h(j), 'CLim', [0, 1]);
    % %     % title('Sliding Check [boolean]')
    % % elseif j == 9
    % %     % set(h(j), 'CLim', [0, 1]);
    % %     % title('Friction Coefficient')
    % % end

    set(h(j), 'CLim', [0, 0.5]);

    if j == 2
        set(h(j), 'CLim', [0, 2.0]);  % Use 'set' to change the color limit
        % title('Displacement Magnitude [mm]')
    elseif j == 3
        set(h(j), 'CLim', [0, 1.0]);
        % title('Shear Stress Magnitude [MPa]')
    elseif j == 4
        set(h(j), 'CLim', [-11, 0.1]);
        % title('X Displacement [mm]')
    elseif j == 5
        set(h(j), 'CLim', [-0.5, 0.5]);
        % title('Y Displacement [mm]')
    elseif j == 6
        set(h(j), 'CLim', [-2, 0.2]);
        % title('Longitudinal Shear Stress [MPa]')
    elseif j == 7
        set(h(j), 'CLim', [-0.2, 0.2]);
        % title('Lateral Shear Stress [MPa]')
    elseif j == 8
        set(h(j), 'CLim', [0, 1]);
        % title('Sliding Check [boolean]')
    elseif j == 9
        set(h(j), 'CLim', [0, 1.3]);
        % title('Friction Coefficient')
    end

end

pause(1);

% Loop to update the image data
for t = 1:length(time_span)

    pause(0.01);  % Pause to control animation speed
    
    for j = 1:9
        % Update the data in the existing image plot rather than redrawing it
        set(im(j), 'CData', sim_solution(:, :, t, j));  % Update the data
        % title(h(j), ['Time Step: ', num2str(t)]);

        % % Adjust color limits based on time step and subplot index
        % if t <= 400
        %     if j == 2
        %         set(h(j), 'CLim', [0, 10]);  % Use 'set' to change the color limit
        %     elseif j == 3
        %         set(h(j), 'CLim', [0, 0.8]);
        %     elseif j == 4
        %         set(h(j), 'CLim', [0, 10]);
        %     elseif j == 5
        %         set(h(j), 'CLim', [0, 0.1]);
        %     elseif j == 6
        %         set(h(j), 'CLim', [-0.6, 0]);
        %     elseif j == 7
        %         set(h(j), 'CLim', [0, 1]);
        %     end
        % else
        %     if j == 2
        %         set(h(j), 'CLim', [0, 8]);
        %     elseif j == 3
        %         set(h(j), 'CLim', [0, 0.8]);
        %     elseif j == 4
        %         set(h(j), 'CLim', [0, 8]);
        %     elseif j == 5
        %         set(h(j), 'CLim', [0, 0.1]);
        %     elseif j == 6
        %         set(h(j), 'CLim', [-0.8, 0.8]);
        %     elseif j == 7
        %         set(h(j), 'CLim', [0, 1]);
        %     end
        % end
    end
    
    % Capture the plot as a frame
    frame = getframe(gcf);
    
    % Write the frame to the video file
    writeVideo(v, frame);
end

% Close the video file
close(v);

fprintf("Animation Successfully saved as Video!\n")
%%


% Plot an initial frame and store the handle to the image object
im = imagesc(y_vals, x_vals, P_grid_moving(:, :, 1));  % Initial frame
c = colorbar;
c.Label.String = "Pressure [MPa]";
title("Pressure Distribution");
xlabel('Lateral y-direction [mm]');
ylabel('Longitudinal x-direction [mm]');

% Loop to update the image data
for t = 1:length(time_span)

    pause(0.01);  % Pause to control animation speed
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'CData', P_grid_moving(:, :, t));  % Update the data
end