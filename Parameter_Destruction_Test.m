set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear; close all; clc

numBrushes  = 20;
isRolling   = true;
fs_sim      = 1e3;
fs_save     = 1e3;
t_initial   = 0;
t_final     = 40;

[model_input, sim_solution] = main(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );
v0 = model_input.v0;
omega = model_input.omega;
dt_sim = model_input.dt_sim;
Fz = model_input.Fz;
K = min(size(model_input.omega), [], 'all');