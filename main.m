function [model_input, sim_solution] = main(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final)
     
    arguments (Input)
        numBrushes  (1, 1) single
        isRolling   (1, 1) logical
        fs_sim      (1, 1) single
        fs_save     (1, 1) single
        t_initial   (1, 1) single
        t_final     (1, 1) single
    end

    arguments (Output)
        model_input     (1, 1) struct
        sim_solution    (1, :) cell
    end
    
    total_time = tic; 
    fprintf("Starting simulation. \n")
    
    model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);
    
    fprintf("Initialised brush model simulation in %.2fs! \n", toc(total_time))
    % compileMex(model_input)
    
    K = min(size(model_input.omega), [], 'all');
    sim_solution = cell(1, K);
    
    %%%%%%%% Simulate brush model %%%%%%%%%%%%%%%%%%%%%
    for i = 1:K
        sim_solution{i} = simulateBrushModel_V2(model_input);
    end
    fprintf("Finished simulation in %.2fs! \n", toc(total_time))
    load chirp.mat%gong.mat
    sound(y, 1 * Fs)
end

