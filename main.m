function [model_input, sim_solution] = main(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final)
    clearvars; clc;
    close all;
    
    total_time = tic; 
    fprintf("Starting simulation. \n")
    
    model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);
    
    fprintf("Initialised brush model simulation in %.2fs! \n", toc(pressure_start))
    
    K = min(size(model_input.omega), [], 'all');
    sim_solution = cell(1, K);
    
    %%%%%%%% Simulate brush model %%%%%%%%%%%%%%%%%%%%%
    for i = 1:K
        sim_solution{i} = simulateBrushModel_V2(model_input);
    end
    toc(total_time)
    load gong.mat
    sound(y, Fs)
end

