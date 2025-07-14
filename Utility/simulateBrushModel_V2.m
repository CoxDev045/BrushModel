function varargout = simulateBrushModel_V2(model_input) %#codegen -args
                                                             
                                                            
    arguments (Input)
        model_input    (1, 1) struct
    end
    % % Declare stdout and stderr as opaque C pointers
    % % These tell coder.ceval how to refer to the standard C stream pointers.
    % coder.cinclude('<stdio.h>'); % Include stdio.h for FILE* and fflush
    % stdout_ptr = coder.opaque('FILE*', 'stdout');
    
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
                      'v0'          ,...
                      're'          ,...
                      'alpha'       ,...
                      'omega_z'     ,...
                      'X'           ,...
                      'Y'           ,...
                      'dA'};

    % Validate inputs
    validateInputs(requiredFields, model_input);

    % Read variables from input struct to minimise accessing struct through
    % loop. (Minor speed increase over accessing struct at each iteration
    % in loop)
    omega   = model_input.omega(:, 1);  % Add Index at omega for rolling 
    omega_z = model_input.omega_z;
    re      = model_input.re;
    v0      = model_input.v0(:, 1); % Add index at v0 for sliding,
    alpha   = model_input.alpha;
    dt_sim  = model_input.dt_sim;
    % Simulate Measurement noise
    noiseVar = single(0.005);
    

    % Initialise solution structs
    [sim_solution, bool_array] = initialiseSolutionStructs(model_input, nargout);
    % Get field names from solution struct
    sim_sol_fieldnames = fieldnames(sim_solution);
    len_sol_fieldNames = int32(length(sim_sol_fieldnames));
    % Initialise steps for counter (default is 10 steps)
    progress_steps = single( round(linspace(1, model_input.LenTime_sim, 10)) ); % 10 checkpoints
    
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

    % counter for saving results
    j = single(1);
    for i = int32(1):int32(model_input.LenTime_sim)

        tempPress = shiftPressure(model_input.X, model_input.Y, ...
                                  model_input.P_grid, ...
                                  shift_amount(i), ...  % Remove index at shift_amount for sliding
                                  maxX, minX, brushArray.p_0); 

        % % Since v0 is an constant multiple of omega, just calculate v0
        % % instead of saving the values to an array
        % v0 = model_input.omega(i, 1) * model_input.re / (model_input.SR + 1);

        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        brushArray = brushArray.update_brush(tempPress, ... 
                                             omega(i), ...  % Add Index at omega for rolling
                                             omega_z, ...
                                             re, ...     
                                             v0(i), ...     % Add index at v0 for sliding,
                                             alpha, ...  
                                             dt_sim ... 
                                             ); 
                                                  
        % Check whether we should save the output or not
        shouldSave = mod(i, model_input.dt_ratio) == 0;
        %%%%%%%%%%%%%%% Sample Data %%%%%%%%%%%%%%%%%%%%%%%%
        if shouldSave
            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution.PressGrid(:, j) = tempPress;
            
            for k = int32(4):len_sol_fieldNames
                % Simulate Measurement noise
                measurementNoise = randn(1, 1, 'single') * noiseVar;
                sim_solution.(sim_sol_fieldnames{k})(:, j) = brushArray.(sim_sol_fieldnames{k}) + measurementNoise;
            end
            
            % Both the fields has to be non-empty for it to save
            if ~isempty(bool_array.isSliding) && ~isempty(bool_array.hasPassed)
                bool_array.isSliding(:, j) = brushArray.slide ;
                bool_array.hasPassed(:, j) = brushArray.passed;
            end
            %%%%%%%%%%%%%%% Increment counter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            j = single(j + 1);
        end
        %%%%%%%%%%%%%% Update progress at defined steps %%%%%%%%%%%%%%%%%%%
        if any(i == progress_steps)
            fprintf('%d%% completed (%d/%d steps)\n', ...
                        round(100 * i / int32(model_input.LenTime_sim)), i, int32(model_input.LenTime_sim));
        end
    end

    %%%%%%%%%%%%%% Calculate Magnitude of Displacement and Stresses %%%%
    sim_solution.TotalDisplacement = hypot(sim_solution.delta_x, sim_solution.delta_y);
    sim_solution.TotalStress = hypot(sim_solution.tauX, sim_solution.tauY);


    if nargout >= 1
        varargout{1} = sim_solution;
    end
    if nargout >= 2
        varargout{2} = bool_array;
    end

end

function validateInputs(requiredFields, model_input)
   
    % Initialise flag with false (no errors thus far)
    flag = false;
    % Initialise message with success message
    message = 'All required fields are present. Proceeding with simulation...\n';
    % Initialise missing field with empty string
    missingField = '';
    % Loop through all fields in requiredFields and check whether
    % model_input has them
    for i = 1:length(requiredFields)
        currentField = requiredFields{i};
        if ~isfield(model_input, currentField)
            % If field is missing, set flag to true
            flag = true;
            % Log the missing field
            missingField = currentField;

            % Generate Error message.
            message = sprintf('ERROR: Missing required field %s.\n', currentField);
        end
    end

    if flag
        % If there is a flag, print message and throw assertion error
        fprintf(message)
        assert(false, ['Missing required field: ', missingField, ' in model_input struct.']);
    else
        % If there is no flag, print success message
        fprintf(message);
    end    
end

function [SolStruct, boolStruct] = initialiseSolutionStructs(model_input, numOutputs)
    % Allocate memory based on dynamic sizes
    numBrushes = model_input.numElems^2;
    arraySize = model_input.LenTime_save;

    SolStruct = struct('PressGrid',          zeros(numBrushes, arraySize, 'single'),...
                       'TotalDisplacement',  zeros(numBrushes, arraySize, 'single'),...
                       'TotalStress',        zeros(numBrushes, arraySize, 'single'),...
                       'delta_x',            zeros(numBrushes, arraySize, 'single'),...
                       'delta_y',            zeros(numBrushes, arraySize, 'single'),...
                       'tauX',               zeros(numBrushes, arraySize, 'single'),...
                       'tauY',               zeros(numBrushes, arraySize, 'single'),...
                       'mu',                 zeros(numBrushes, arraySize, 'single'),...
                       'vs',                 zeros(numBrushes, arraySize, 'single') ...
                       );
    if numOutputs >= 2  
        boolStruct = struct('isSliding', false(numBrushes, arraySize),...
                            'hasPassed', false(numBrushes, arraySize) ...
                          );
    else
        boolStruct = struct('isSliding', [],...
                            'hasPassed', [] ...
                          );
    end
end