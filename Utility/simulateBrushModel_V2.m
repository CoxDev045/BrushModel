function varargout = simulateBrushModel_V2(model_input) %#codegen -args
                                                             
                                                            
    arguments (Input)
        model_input    (1, 1) struct
    end

    % Define the required fields for your model_input struct
    requiredFields = {'numElems',...
                      'LenTime_sim' ,...
                      'LenTime_save',...
                      'dt_sim'      ,...
                      'dt_save'     ,...
                      'dt_ratio'    ,...
                      'Press'       ,...
                      'omega'       ,...
                      'SR'          ,...
                      're'          ,...
                      'alpha'       ,...
                      'omega_z'     ,...
                      'X'           ,...
                      'Y'           ,...
                      'v0'};

    for i = 1:length(requiredFields)
        currentField = requiredFields{i};
        if ~isfield(model_input, currentField)
            % For codegen, using 'assert' with a simpler message is often preferred
            % for fundamental checks that indicate a design/input problem.

            % This will cause an assertion failure in the generated C code.
            coder.ceval('printf', 'ERROR: Missing required field %s.\n', currentField); % For console output
            assert(false, ['Missing required field: ', currentField, ' in model_input struct.']); 
            % The assert message itself will be a fixed string at compile time.
        end
    end

    % If all checks pass, continue with your simulation logic
    coder.ceval('printf', 'All required fields are present. Proceeding with simulation...\n');

    % Define ratio of sim time to sample time
    dt_ratio = int32(dt_save / dt_sim);

    % Allocate memory based on dynamic sizes
    % sim_solution = zeros(max_numElems * max_numElems, max_LenSaveTime, 9);
    sim_solution = struct('PressGrid', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'TotalDisplacement', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'TotalStress', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'delta_x', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'delta_y', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'tauX', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'tauY', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'mu', zeros(max_numElems * max_numElems, max_LenSaveTime),...
                          'vs', zeros(max_numElems * max_numElems, max_LenSaveTime) ...
                          );
    sim_sol_fieldnames = fieldnames(sim_solution);

    if nargout > 1  
        % bool_array = false(max_numElems * max_numElems, max_LenSaveTime, 2);
        bool_array = struct('', false(max_numElems * max_numElems, max_LenSaveTime),...
                            '', false(max_numElems * max_numElems, max_LenSaveTime) ...
                          );
    end

    progress_steps = round(linspace(1, max_LenSimTime, 10)); % 10 checkpoints
    
    % Initialise grid of brushes
    brushArray = BrushVec_CPP(model_input.X(:), model_input.Y(:), ...
                              model_input.P_grid(:), model_input.numElems^2);
   
    % % % For sliding, assume omega is zero
    % % omega = 0;

    %%%%%%%%%%%%%%%% Calculations for Pressure distribution %%%%%%%%%%%%%%%
    maxX = max(X, [], 'all');
    % maxY = max(Y, [], 'all');
    minX = min(X, [], 'all');
    % minY = min(Y, [], 'all');

    % The amount the pressure distribution will have shifted due to rolling
    shift_amount = cumsum(model_input.omega * model_input.dt_sim * model_input.re);
    % % shift_amount = 0; % For sliding the shift due to rolling is zero 

    % counter for saving results
    j = 1;
    for i = int32(1):int32(model_input.max_LenSimTime)

        tempPress = shiftPressure(model_input.X, model_input.Y, ...
                                  model_input.P_grid, ...
                                  shift_amount(i), ...  % Remove index at shift_amount for sliding
                                  maxX, minX, brushArray.p_0); 

        % Since v0 is an constant multiple of omega, just calculate v0
        % instead of saving the values to an array
        v0 = omega(i) * re / (SR + 1);

        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        brushArray = brushArray.update_brush(tempPress, model_input.omega(i), ...  % Add Index at omega for rolling 
                                             model_input.omega_z, model_input.re, ...
                                             model_input.v0, ... % Add index at v0 for sliding,
                                             model_input.alpha, ...
                                             model_input.dt_sim); 
                                                                                                  

        if mod(i, dt_ratio) == 0
            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution(:, j, 1) = tempPress;
            
            for k = 4:length(sim_sol_fieldnames)
                sim_solution.(sim_sol_fieldnames{k}) = brushArray.(sim_sol_fieldnames{k});
            end

            if nargout > 1
                bool_array(:, j, :) = cat(2, brushArray.slide, brushArray.passed);
            end
            %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
            sim_solution.TotalDisplacement = hypot(sim_solution.delta_x, sim_solution.delta_y);
            sim_solution.TotalStress = hypot(sim_solution.tauX, sim_solution.tauY);

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