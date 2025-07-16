clear; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%

clearvars; clc;
close all;

total_time = tic;
pressure_start = tic; 
fprintf("Starting simulation. \n")
fs_save = 1e2; % Hz
fs_sim = 1e3; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 120;
isRolling = false;
t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

fprintf("Initialised brush model simulation in %.2fs! \n", toc(pressure_start))

K = min(size(model_input.omega), [], 'all');
sim_solution = cell(1, K);

%%%%%%%% Simulate brush model %%%%%%%%%%%%%%%%%%%%%
for i = 1:K
    sim_solution{i} = simulateBrushModel_V2(model_input);
end

load gong.mat
sound(y, Fs)