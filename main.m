clear; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%

clearvars; clc;
close all;



total_time = tic;
pressure_start = tic; 
fprintf("Starting simulation. \n")
fs_save = 1e3; % Hz
fs_sim = 1e5; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 100;
isRolling = true;
t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

fprintf("Initialised brush model simulation in %.2fs! \n", toc(pressure_start))
%
%%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%
% % Press = reshape(P_grid.shifted(:, 1), numElems, numElems);
sliding = false;
K = min(size(model_input.omega), [], 'all');
sim_solution = cell(1, K);
for i = 1:K
    sim_solution{i} = simulateBrushModel_V2_mex(model_input);
end
% tic;
% sim_solution = simulateBrushModel_V2_mex(model_input);
% toc



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

 