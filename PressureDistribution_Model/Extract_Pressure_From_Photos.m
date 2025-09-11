set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear; close all;clc

PhotoDir = 'C:\Users\coxde\OneDrive\Masters\StaticTyreTestRig';

Files = dir(PhotoDir);
Files = Files(~[Files.isdir]);
FileNames = {Files.name};
FileNames = FileNames(1:end-3);


sigma = 5;
N = 100;

FinalImages = zeros(N, N, length(FileNames));

for i = 1:length(FileNames)
    fullName = fullfile(PhotoDir, FileNames{i});
    im = imread(fullName);
    
    % Convert to grayscale
    im_gray = rgb2gray(im);
    
    % Crop image to desired size
    im_cropped = imcrop(im_gray);

    imshow(im_cropped)

    % Apply gaussian blur before downsampling
    im_blurred = imgaussfilt(im_cropped, sigma);

    % Downsample to size NxN
    im_downsampled = imresize(im_blurred, [N, N]);

    imshow(im_downsampled)

    % Normalise values to range [0, 1]
    im_normalised = 1 - im2double(im_downsampled);

    % Save final image
    FinalImages(:, :, i) = im_normalised;
end

%%
clear;close all;clc
load("PressImages.mat") 

X_lengths = [70, 85, 119.6, 134.7, 146.5, 166, 196, 198, 205];
Y_lengths = [68.9, 68.9, 104.9, 115.7, 122.3, 148, 164, 168, 175];
Z_loads = [1500, 2000, 2500, 3000, 3500, 4000, 5000, 5500, 6000];

A_mat = [Z_loads(:), Z_loads(:).^2, log10(Z_loads(:))];

% Ax = b
X_vec_x = A_mat \ X_lengths(:);
X_vec_y = A_mat \ Y_lengths(:);

figure
T = tiledlayout('flow');
T.TileSpacing = 'tight';
T.Padding = "compact";

for i = 1:9
    nexttile
    imshow(FinalImages(:, :, i))
end

figure
hold on
scatter(Z_loads, X_lengths)
plot(Z_loads, A_mat * X_vec_x)
scatter(Z_loads, Y_lengths)
plot(Z_loads, A_mat * X_vec_y)
grid on
legend('Measured X', 'X Model', 'Measured Y', 'Y Model')
xlabel('Vertical Load [N]')
ylabel('Length [mm]')
%%
X_min = -1 * X_lengths./2;
X_max = X_lengths./2;
Y_min = -1 * Y_lengths./2;
Y_max = Y_lengths./2;

X_grids = cell(1, length(X_lengths));
Y_grids = cell(size(X_grids));

for i = 1:length(X_lengths)
    x_vals = linspace(X_min(i), X_max(i), 100);
    y_vals = linspace(Y_min(i), Y_max(i), 100);
    [X, Y] = meshgrid(x_vals, y_vals);

    X_grids{i} = X;
    Y_grids{i} = Y;
end

% Find the overall min and max of all X and Y coordinates
all_X_values = cell2mat(X_grids');
all_Y_values = cell2mat(Y_grids');
min_X_global = min(all_X_values(:));
max_X_global = max(all_X_values(:));
min_Y_global = min(all_Y_values(:));
max_Y_global = max(all_Y_values(:));

% Define the resolution of your master canvas
num_points_master = 500; % Choose a high enough resolution to avoid data loss
x_master_vec = linspace(min_X_global, max_X_global, num_points_master);
y_master_vec = linspace(min_Y_global, max_Y_global, num_points_master);

% Define your desired final size
N = 60;
x_plot_vec = linspace(min_X_global, max_X_global, N);
y_plot_vec = linspace(min_Y_global, max_Y_global, N);

[X_master, Y_master] = meshgrid(x_master_vec, y_master_vec);

[X_plot, Y_plot] = meshgrid(x_plot_vec, y_plot_vec);

% Create a cell array to store the padded pressure fields
padded_pressure_slices = cell(1, size(FinalImages,3));

% Loop through each unique pressure field
for i = 1:size(FinalImages, 3)
    current_P = FinalImages(:, :, i);
    current_X = X_grids{i};
    current_Y = Y_grids{i};
    
    % Interpolate the current pressure field onto the master canvas
    % Outside the current field's domain, the result will be NaN
    P_temp = interp2(current_X, current_Y, current_P, X_master, Y_master, 'linear');
    
    % % Replace the NaNs with zeros to complete the padding
    % P_temp(isnan(P_temp)) = 0;
    
    % Store the padded result
    padded_pressure_slices{i} = P_temp;
end

% Create a cell array to store the final NxN pressure fields
final_pressure_slices = cell(size(padded_pressure_slices));

% Loop through each padded pressure field
for i = 1:length(padded_pressure_slices)
    current_P_padded = padded_pressure_slices{i};
    
    % Resample the padded field to your final NxN size
    final_pressure_slices{i} = imresize(current_P_padded, [N, N], 'lanczos2');
end

% Now final_pressure_slices contains all your pressure fields in a
% uniform NxN format, with correct physical scaling.
figure;
ax = gca;
hold(ax, 'on');

% Loop through each padded pressure field and plot it at its load level
for i = 1:length(final_pressure_slices)
    current_P_padded = final_pressure_slices{i};
    current_load = Z_loads(i);
    
    % Create a Z matrix for plotting, so all points are at the correct load level
    Z_plot = current_load * ones(size(X_plot));
    
    % Plot the filled contour slice on the 3D axes
    % We use 'surf' with a specified color data to get a filled surface plot
    surf(ax, X_plot, Y_plot, Z_plot, current_P_padded, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.8);
end

% --- 4. Final Plot Settings ---
colormap(ax, 'jet');
% cb = colorbar(ax, 'Location', 'westoutside');
% ylabel(cb, 'Pressure (units)');
xlabel(ax, 'X Position (m)');
ylabel(ax, 'Y Position (m)');
zlabel(ax, 'Vertical Load (N)');
title(ax, 'Pressure Distribution Growth with Increasing Load');
view(ax, 3); % Set a 3D view
grid(ax, 'on');
axis(ax, 'tight');
hold(ax, 'off');