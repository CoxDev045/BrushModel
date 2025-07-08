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
                      'P_grid'      ,...
                      'omega'       ,...
                      'SR'          ,...
                      're'          ,...
                      'alpha'       ,...
                      'omega_z'     ,...
                      'X'           ,...
                      'Y'           ,...
                      'v0'};

    validateInputs(requiredFields, model_input);
    
    [sim_solution, bool_array] = initialiseSolutionStructs(model_input, nargout);

    sim_sol_fieldnames = fieldnames(sim_solution);

    progress_steps = round(linspace(1, model_input.LenTime_save, 10)); % 10 checkpoints
    
    % Initialise grid of brushes
    brushArray = BrushVec_CPP(model_input.X(:), model_input.Y(:), ...
                              model_input.P_grid(:), ...
                              model_input.numElems^2);
   
    % % % For sliding, assume omega is zero
    % % omega = 0;

    %%%%%%%%%%%%%%%% Calculations for Pressure distribution %%%%%%%%%%%%%%%
    maxX = max(model_input.X, [], 'all');
    % maxY = max(Y, [], 'all');
    minX = min(model_input.X, [], 'all');
    % minY = min(Y, [], 'all');

    % The amount the pressure distribution will have shifted due to rolling
    shift_amount = cumsum(model_input.omega * model_input.dt_sim * model_input.re);
    % % shift_amount = 0; % For sliding the shift due to rolling is zero 

    for i = int32(1):int32(model_input.LenTime_save)

        tempPress = shiftPressure(model_input.X, model_input.Y, ...
                                  model_input.P_grid, ...
                                  shift_amount(i), ...  % Remove index at shift_amount for sliding
                                  maxX, minX, brushArray.p_0); 

        % % Since v0 is an constant multiple of omega, just calculate v0
        % % instead of saving the values to an array
        % v0 = model_input.omega(i, 1) * model_input.re / (model_input.SR + 1);

        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        brushArray = brushArray.update_brush(tempPress, model_input.omega(i, 1), ...  % Add Index at omega for rolling 
                                             model_input.omega_z, model_input.re, ...
                                             model_input.v0(i, 1), ... % Add index at v0 for sliding,
                                             model_input.alpha, ...
                                             model_input.dt_sim ...
                                             ); 
                                                                                                  

        if mod(i, model_input.dt_ratio) == 0
            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution.PressGrid = tempPress;
            
            for k = 4:length(sim_sol_fieldnames)
                sim_solution.(sim_sol_fieldnames{k}) = brushArray.(sim_sol_fieldnames{k});
            end

            if isstruct(bool_array)
                bool_array.isSliding = brushArray.slide ;
                bool_array.hasPassed = brushArray.passed;
            end
            %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
            sim_solution.TotalDisplacement = hypot(sim_solution.delta_x, sim_solution.delta_y);
            sim_solution.TotalStress = hypot(sim_solution.tauX, sim_solution.tauY);
        end
        %%%%%%%%%%%%%% Update progress at defined steps %%%%%%%%%%%%%%%%%%%
        if any(i == progress_steps)
            fprintf('%d%% completed (%d/%d steps)\n', round(100 * i / int32(model_input.LenTime_sim)), i, int32(model_input.LenTime_sim));
        end
    end

    varargout{1} = sim_solution;
    if nargout > 1
        varargout{2} = bool_array;
    end

end

function validateInputs(requiredFields, model_input)
    for i = 1:length(requiredFields)
        currentField = requiredFields{i};
        if ~isfield(model_input, currentField)
            % For codegen, using 'assert' with a simpler message is often preferred
            % for fundamental checks that indicate a design/input problem.

            % This will cause an assertion failure in the generated C code.
            if coder.target('MEX') || coder.target('Rtw')
                coder.ceval('printf', 'ERROR: Missing required field %s.\n', currentField); % For console output
            else
                fprintf('ERROR: Missing required field %s.\n', currentField);
                assert(false, ['Missing required field: ', currentField, ' in model_input struct.']); 
                % The assert message itself will be a fixed string at compile time.
            end
        end
    end

    if coder.target('MEX') || coder.target('Rtw')
        % If all checks pass, continue with your simulation logic
        coder.ceval('printf', 'All required fields are present. Proceeding with simulation...\n');
    else
        fprintf('All required fields are present. Proceeding with simulation...\n');
    end
end

function [sim_solution, bool_array] = initialiseSolutionStructs(model_input, numOutputs)
    % Allocate memory based on dynamic sizes
    % sim_solution = zeros(max_numElems * max_numElems, max_LenSaveTime, 9);
    sim_solution = struct('PressGrid',          zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'TotalDisplacement',  zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'TotalStress',        zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'delta_x',            zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'delta_y',            zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'tauX',               zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'tauY',               zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'mu',                 zeros(model_input.numElems^2, model_input.LenTime_save, 'single'),...
                          'vs',                 zeros(model_input.numElems^2, model_input.LenTime_save, 'single') ...
                          );
    if numOutputs > 1  
        % bool_array = false(max_numElems * max_numElems, max_LenSaveTime, 2);
        bool_array = struct('isSliding', false(model_input.numElems^2, model_input.LenTime_save),...
                            'hasPassed', false(model_input.numElems^2, model_input.LenTime_save) ...
                          );
    else
        bool_array = [];
    end
end