function [sim_solution, bool_array] = simulateBrushModel(max_numElems, max_LenTime, P_grid, time_span, omega, v0, dt, re, alpha, omega_z, X, Y) %#codegen -args

    arguments (Input)
        max_numElems (1, 1) uint16 {coder.Constant} % **Must be a constant**
        max_LenTime (1, 1) single {coder.Constant} % **Must be a constant**
        P_grid      (:, :) single
        time_span   (:, 1) single
        omega       (:, 1) single
        v0          (:, 1) single
        dt          (1, 1) single
        re          (1, 1) single
        alpha       (1, 1) single
        omega_z     (1, 1) single
        X           (:, 1) single
        Y           (:, 1) single
    end
    
    % Define variable-size arrays with user-defined max sizes
    coder.varsize('P_grid', [max_numElems, max_LenTime], [1, 1]);
    coder.varsize('time_span', [max_LenTime, 1], [1, 0]);
    coder.varsize('omega', [max_LenTime, 1], [1, 0]);
    coder.varsize('sim_solution', [max_numElems, max_LenTime, 7], [1, 1, 0]);
    coder.varsize('slide_array', [max_numElems, max_LenTime], [1, 1]);

    % Ensure input arrays do not exceed the set max sizes
    assert(size(P_grid, 1) <= max_numElems);
    assert(size(P_grid, 2) <= max_LenTime);
    assert(size(time_span, 1) <= max_LenTime);
    assert(size(omega, 1) <= max_LenTime);


    % Allocate memory based on dynamic sizes
    sim_solution = zeros(max_numElems, max_LenTime, 7, 'single');
    bool_array = false(max_numElems, max_LenTime, 2);

    progress_steps = round(linspace(1, max_LenTime, 10)); % 10 checkpoints
    % sim_solution(:, :, 1) = shiftPressure(P_grid, X, Y, omega, time_span, dt, re);

    brushArray = BrushVec_CPP(X, Y, P_grid(:, 1), max_numElems);

    for i = 1:max_LenTime
        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        brushArray = brushArray.update_brush(P_grid(:, i), omega(i), omega_z, re, v0(i), alpha, dt);
        %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sim_solution(:, i, 3:7) = cat(3, brushArray.delta_x, brushArray.delta_y, ...
                                      brushArray.tauX, brushArray.tauY, ...
                                      brushArray.mu);
        bool_array(:, i, 1:2) = cat(3, brushArray.slide, brushArray.passed);
    
        %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
        sim_solution(:, i, 1) = hypot(sim_solution(:, i, 3), sim_solution(:, i, 4));
        sim_solution(:, i, 2) = hypot(sim_solution(:, i, 5), sim_solution(:, i, 6));
        %%%%%%%%%%%%%% Increment counter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update progress at defined steps
        if any(i == progress_steps)
            fprintf('%d%% completed (%d/%d steps)\n', round(100 * i / length(time_span)), i, length(time_span));
        end
    end

end
