function varargout = simulateBrushModel_V2(max_numElems, max_LenSimTime, max_LenSaveTime, ...
                                                            P_grid, dt_sim, dt_save,...
                                                            omega, SR, ... % Uncomment for rolling
                                                            re, alpha, omega_z, X, Y) %#codegen -args
                                                            % v0, ... % Uncomment for sliding 
                                                            
    arguments (Input)
        max_numElems (1, 1) uint16  % **Must be a constant**
        max_LenSimTime (1, 1) double % **Must be a constant**
        max_LenSaveTime (1, 1) double  % **Must be a constant**
        P_grid      (:, :) double
        dt_sim      (1, 1) double
        dt_save     (1, 1) double
        % Change to omega for rolling
        omega          (:, :) double
        % Comment out for sliding
        SR          (1, 1) double
        re          (1, 1) double
        alpha       (1, 1) double
        omega_z     (1, 1) double
        X           (:, :) double
        Y           (:, :) double
    end
    

    % Ensure input arrays do not exceed the set max sizes
    assert(size(P_grid, 1) <= max_numElems);
    assert(size(P_grid, 2) <= max_numElems);
    % assert(size(omega, 1) <= max_LenSimTime);
    assert(size(omega, 1) <= max_LenSimTime); % ** Change to v0 for sliding **
    assert(size(X, 1) <= max_numElems);
    assert(size(X, 2) == max_numElems);
    assert(size(Y, 1) <= max_numElems);
    assert(size(Y, 2) == max_numElems);

    % Define ratio of sim time to sample time
    dt_ratio = int32(dt_save / dt_sim);

    % Allocate memory based on dynamic sizes
    sim_solution = zeros(max_numElems * max_numElems, max_LenSaveTime, 9);
    if nargout > 1
        bool_array = false(max_numElems * max_numElems, max_LenSaveTime, 2);
    end

    progress_steps = round(linspace(1, max_LenSimTime, 10)); % 10 checkpoints
    
    % Initialise grid of brushes
    brushArray = BrushVec_CPP(X(:), Y(:), P_grid(:), max_numElems * max_numElems);
   
    % % % For sliding, assume omega is zero
    % % omega = 0;

    %%%%%%%%%%%%%%%% Calculations for Pressure distribution %%%%%%%%%%%%%%%
    maxX = max(X, [], 'all');
    maxY = max(Y, [], 'all');
    minX = min(X, [], 'all');
    minY = min(Y, [], 'all');
    
    % % numBrushes = sqrt(max_numElems);
    % % PressX = reshape(X, numBrushes, numBrushes);
    % % PressY = reshape(X, numBrushes, numBrushes);

    % The amount the pressure distribution will have shifted due to rolling
    shift_amount = cumsum(omega * dt_sim * re);
    % % shift_amount = 0; % For sliding the shift due to rolling is zero 

    % counter for saving results
    j = 1;
    for i = int32(1):int32(max_LenSimTime)

        tempPress = shiftPressure(X, Y, P_grid, shift_amount(i), maxX, minX, brushArray.p_0); % Remove index at shift_amount for sliding

        % Since v0 is an constant multiple of omega, just calculate v0
        % instead of saving the values to an array
        v0 = omega(i) * re / (SR + 1);

        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        brushArray = brushArray.update_brush(tempPress, omega(i), omega_z, re, v0, alpha, dt_sim); % Index at v0 for sliding,
                                                                                                   % Index at omega for rolling 

        if mod(i, dt_ratio) == 0
            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution(:, j, 1) = tempPress;
            sim_solution(:, j, 4:9) = cat(3, brushArray.delta_x, brushArray.delta_y, ...
                                          brushArray.tauX, brushArray.tauY, ...
                                          brushArray.mu, brushArray.vs);
            if nargout > 1
                bool_array(:, j, :) = cat(2, brushArray.slide, brushArray.passed);
            end
            %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
            sim_solution(:, j, 2) = hypot(sim_solution(:, j, 4), sim_solution(:, j, 5));
            sim_solution(:, j, 3) = hypot(sim_solution(:, j, 6), sim_solution(:, j, 7));

            %%%%%%%%%%%%%%% Increment counter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            j = j + 1;
        end
        %%%%%%%%%%%%%% Update progress at defined steps %%%%%%%%%%%%%%%%%%%
        if any(i == progress_steps)
            fprintf('%d%% completed (%d/%d steps)\n', round(100 * i / int32(max_LenSimTime)), i, int32(max_LenSimTime));
        end
    end

    varargout{1} = sim_solution;
    if nargout > 1
        varargout{2} = bool_array;
    end

end

function tempPress = shiftPressure(X, Y, P, shift_amount, max_val, min_val, Pmin)

    % Corrects for the tendency of the mod function to shift the leading edge to the centre of the contact patch
    X_shift_corrector = mod(X, max_val - abs(min_val));

    % Calculate shift amount due to tyre rolling
    X_shifted = mod(X_shift_corrector + shift_amount, max_val - (min_val)) + min_val;
    
    % Linearly interpolate pressure along grid
    tempPress = interp2(X, Y, P, X_shifted, Y, 'linear');

    % Mask Pressure grid to remove negative or small values (done to smooth out final outputs)
    tempPress = max(tempPress(:), Pmin);
end