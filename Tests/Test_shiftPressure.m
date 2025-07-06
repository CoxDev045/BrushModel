clearvars; clc;
close all;

total_time = tic;
pressure_start = tic; 
fprintf("Starting simulation. \n")

% Contact Width Parameters
a1 = 1.03;
a2 = 0.44;
a3 = 195 * 0.35;

% Vertical Load
Fz = 3000;

% Contact Width (y)
a = a1 * Fz^a2 + a3;
% Contact Length (x)
b = 195 / 2;


numBrushes = 20;

re = 15 * 25.4 + 0.55 * 195; % 195 / 55 R15
alpha = deg2rad(0);

fs_sim = 1e5; % Hz
dt_sim = 1 / fs_sim;

fs_save = 1e3; % Hz
dt_save = 1 / fs_save;

omega_z = 0;

t_initial = 0;
t_final = 3;
t_save = linspace(t_initial, t_final, t_final * fs_save + 1);

rpm = 275;
omega_max = rpm * 2 * pi / 60; % Convert RPM to rad/s

omega = zeros(fs_sim * t_final+1, 1);%[linspace(0, omega_max, fs_sim * t_final/2), repmat(omega_max, 1,  fs_sim * t_final/2+1)];

SR = 0.1;
% v0 = omega * re / (SR + 1);

% Initialise PRessure function coefficients
nx = 1;
ny = 1.5;
n_params = [nx, ny];

lambda_x = 1.25;
lambda_y = 1;
lambda = [lambda_x, lambda_y];

xe = -1; % mm
ye = -1; % mm
COP = [xe, ye];

One_D = false;

x_vals = linspace(-b, b, numBrushes);
y_vals = linspace(-a, a, numBrushes);
dx = x_vals(2) - x_vals(1);
dy = y_vals(2) - y_vals(1);
dA = dx * dy;

[X, Y] = meshgrid(x_vals, y_vals);


[~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
P_grid = real(P_grid);

contourf(X, Y, P_grid .* (P_grid > 0.02), 25)
colorbar
axis equal

% % figure
% % subplot(211)
% % imagesc(P_grid)
% % colorbar
% % title('Uniform')
% % 
% % P_grid = P_grid .* perlin(size(P_grid));
% % 
% % subplot(212)
% % imagesc(P_grid)
% % colorbar
% % title('With noise added')
%%
numBrushes = uint16(numBrushes);
numElems = uint16(numBrushes^2);
LenTime_sim = single(fs_sim * t_final + 1);
LenTime_save = single(fs_save * t_final + 1);
dt_ratio = dt_save / dt_sim;

[P_grid_shifted, P_grid_save] = shiftPressure(numBrushes, numElems, LenTime_sim, dt_ratio, P_grid, X, Y, omega(:), dt_sim, re);
 
fprintf("Generated pressure distribution in %.2fs! \n", toc(pressure_start))

% whos

clear P_grid
P_grid.shifted = P_grid_shifted;
P_grid.save = P_grid_save;

% P_grid.shifted = (P_grid.shifted) + 0.05 * (0.5 - rand(size(P_grid.shifted)));
% P_grid.save = (P_grid.save) + 0.05 * (0.5 - rand(size(P_grid.save)));

save("Sliding_Pressure_fs_sim_1e5_fs_save_1e3_len_3s.mat", 'P_grid', '-v7.3')
%%
if size(P_grid.save, 1) ~= numBrushes 
    pressure_mask = P_grid.save > 0.02;
    P_grid_save = reshape(P_grid.save .* pressure_mask, numBrushes, numBrushes, LenTime_save);
end

% Plot an initial frame and store the handle to the image object
im = imagesc(x_vals, y_vals, P_grid_save(:, :, 1));  % Initial frame
c = colorbar;
c.Label.String = "Pressure [MPa]";
title("Pressure Distribution");
ylabel('Lateral y-direction [mm]');
xlabel('Longitudinal x-direction [mm]');
axis equal

% Loop to update the image data
for t = 1:LenTime_save
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'CData', P_grid_save(:, :, t));  % Update the data

    % Reduce lag
    if mod(t, 10) == 0
        pause(0.1);
    end
end