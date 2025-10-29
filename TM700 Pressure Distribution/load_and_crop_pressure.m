clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')

%%
% Load data
FileName = "TM700Fz560Tr100r2_RawData.csv";
measured_array = Measured_Pressure.load_raw_data(FileName);

% Crop Pressure grid to remove any small and unwanted value that can
% influence the subsampling
% Define how many rows the calibration stick takes up
row = 3;
% Define how many columns the calibration stick takes up
col = 2;
% Define how long the calibration stick is
width = 32;
% Define a threshold for the values to be removed
threshold = 15;
% Crop pressure grid
P_grid_cropped = Measured_Pressure.crop_pressure_grid(measured_array, row, col, width, threshold);

figure
subplot(1,2,1)
contourf(measured_array)
xlabel('Pixels')
ylabel('Pixels')
title('Measured Pressure Field From TM700 with 100% Tread')

subplot(1,2,2)
contourf(P_grid_cropped)
xlabel('Pixels')
ylabel('Pixels')
title('Cropped Pressure Field From TM700 with 100% Tread')
%%
% Calibrate pressure field so that the measured values can be interpreted
Calibration_stick_length = 100; %mm
Calibration_stick_pixels = width; % How many pixels corresponds to the length
Fz = 560 * 9.81; % Vertical load in N
% Here we'll use the same threshold as before
pressure_array = Measured_Pressure.calibrate_pressure_field(P_grid_cropped, Calibration_stick_length, Calibration_stick_pixels, Fz, threshold);

figure
contourf(pressure_array)
xlabel('Pixels')
ylabel('Pixels')
title('Measured Pressure Field From TM700 with 0% Tread')
c = colorbar;
c.Label.String = 'Pressure [MPa]';

%%
% Subsample grid so that it can be used with the brush model
target_size = [10, 10]; % The grid size for the brush model
P_grid_subsampled = Measured_Pressure.subsample_pressure_grid(pressure_array, target_size, Calibration_stick_length, Calibration_stick_pixels);

figure
subplot(1,2,1)
contourf(pressure_array)
xlabel('Pixels')
ylabel('Pixels')
title('Measured Pressure Field From TM700 with 100% Tread')
c = colorbar;
c.Label.String = 'Pressure [MPa]';

subplot(1,2,2)
contourf(P_grid_subsampled)
xlabel('Pixels')
ylabel('Pixels')
title('Subsampled Pressure Field From TM700 with 100% Tread')
c = colorbar;
c.Label.String = 'Pressure [MPa]';

clearvars -except P_grid_subsampled

pressure_file_name = 'TM700Fz560Tr100r2_SubSampled';
save(pressure_file_name, 'P_grid_subsampled');