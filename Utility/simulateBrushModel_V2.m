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
                      'v0'          ,...
                      're'          ,...
                      'alpha'       ,...
                      'omega_z'     ,...
                      'X'           ,...
                      'Y'           ,...
                      'dA'          ,...
                      'isRolling'   ,...
                      'roadProfile' ,...     
                      'changeInRoadHeight'};

    % Validate inputs
    validateInputs(requiredFields, model_input);

    % Read variables from input struct to minimise accessing struct through
    % loop. (Minor speed increase over accessing struct at each iteration
    % in loop)
    omega   = model_input.omega;  % Add Index at omega for rolling 
    omega_z = model_input.omega_z;
    re      = model_input.re;
    v0      = model_input.v0; % Add index at v0 for sliding,
    alpha   = model_input.alpha;
    dt_sim  = model_input.dt_sim;
    isRolling = model_input.isRolling;
    StaticPressGrid = model_input.P_grid;
    StaticContactPatchShape = model_input.P_shape(:);

    numBrushes = single(model_input.numElems);
    Fz = model_input.Fz;

    RoadProfile = model_input.roadProfile;
    ChangeInRoadProfile = model_input.changeInRoadHeight;
    
    StaticX = model_input.X;
    StaticY = model_input.Y;
    dt_ratio = model_input.dt_ratio;
    LenTimeSim = int32(model_input.LenTime_sim);
    % numBrushes = model_input.numElems;
    % LenTimeSave = model_input.LenTime_save;
    noiseVar = 0.001;
    
    % Initialise solution structs
    [sim_solution, bool_array] = initialiseSolutionStructs(model_input, nargout);
    % Get field names from solution struct
    sim_sol_fieldnames = fieldnames(sim_solution);
    len_sol_fieldNames = int32(length(sim_sol_fieldnames));
    % Initialise steps for counter (default is 10 steps)
    progress_steps = single( round(linspace(1, model_input.LenTime_sim, 10)) ); % 10 checkpoints
    
    % Initialise grid of brushes
    brushArray = BrushVec_CPP(StaticX(:), StaticY(:), ...
                              StaticPressGrid(:), ...
                              numBrushes^2);

    %%%%%%%%%%%%%%%% Calculations for Pressure distribution %%%%%%%%%%%%%%%
    maxX = max(model_input.X, [], 'all');
    % maxY = max(Y, [], 'all');
    minX = min(model_input.X, [], 'all');
    % minY = min(Y, [], 'all');
    shift_amount_cumulative = single(0);
    
    Y_road = StaticY + abs(min(StaticY, [], "all"));
    contact_shape = [];
    % ----------------------------------------------------------
    L0      = single(217.214915);
    Fz_0    = single(5716.443731);
    q_L     = single(0.467604);
    W0      = single(225.436018);
    q_W     = single(0.477134);
    % ----------------------------------------------------------
    % counter for saving results
    j = single(1);
    for i = int32(1):int32(model_input.LenTime_sim)
        t_val = single(i-1) * dt_sim;
        X_shifted = reshape(brushArray.x, numBrushes, numBrushes);
        % Calculate contact patch dimensions
        a_current = 0.5 * L0 * (Fz(i) / Fz_0).^(q_L);
        b_current = 0.5 * W0 * (Fz(i) / Fz_0).^(q_W);
        % Calculate Pressure grid based on vertical force
        tempPress = calculatePressure(Fz(i), contact_shape, a_current, b_current, StaticX, StaticY);
        % Create mask for brushes that have vertical load applied
        % verticalMask = (StaticX(:).^2 / (a_current)^2 ) + (StaticY(:).^2 / (b_current)^2 ) <= 1;
        verticalMask = (abs(StaticX(:)) <= a_current) & (abs(StaticY(:)) <= b_current);
        
        tempPress = tempPress .* verticalMask .* StaticContactPatchShape;
        if isRolling
            % The amount the pressure distribution will have shifted due to rolling
            shift_amount = (omega(t_val) * dt_sim * re);
        else
            shift_amount = v0(t_val) * dt_sim;
        end

        % Update Road displacement
        shift_amount_cumulative = shift_amount_cumulative + shift_amount;
        X_road = X_shifted + shift_amount_cumulative + abs(min(X_shifted, [], "all"));
        

        z_disp = RoadProfile(X_road.', Y_road.').';
        z_dot = ChangeInRoadProfile(X_road.', Y_road.').';
        

        %%%%%%%%%%%%%% Use Update Properties and perform update step %%%%%%%
        
        brushArray = brushArray.update_brush(tempPress, ... 
                                             omega, ...  % Add Index at omega for rolling
                                             omega_z, ...
                                             re, ...     
                                             v0, ...     % Add index at v0 for sliding,
                                             alpha, ...  
                                             dt_sim, ... 
                                             t_val,...
                                             z_disp(:),...
                                             z_dot(:)); 
                                                  
        % Check whether we should save the output or not
        shouldSave = mod(i, dt_ratio) == 0;
        %%%%%%%%%%%%%%% Sample Data %%%%%%%%%%%%%%%%%%%%%%%%
        if shouldSave
            %%%%%%%%%%%%%% Save Simulation Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_solution.PressGrid(:, j) = tempPress;
            
            for k = int32(4):len_sol_fieldNames
                % Simulate Measurement noise
                % measurementNoise = randn(1, 1, 'single') * noiseVar;
                sim_solution.(sim_sol_fieldnames{k})(:, j) = brushArray.(sim_sol_fieldnames{k});% + measurementNoise;
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
                        round(100 * i / LenTimeSim), i, LenTimeSim);
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