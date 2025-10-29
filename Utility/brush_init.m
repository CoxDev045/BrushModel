function [model_input] = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final)
    
    % Initialise model input struct
    model_input = struct('numElems',        [],...
                         'LenTime_sim',     [], ...
                         'LenTime_save',    [], ...
                         'dt_sim',          [],...
                         'dt_save',         [],...
                         'dt_ratio',        [],...
                         'P_grid',          [],...
                         'omega',           [],...
                         'SR',              [],...
                         'v0',              [],...
                         're',              [],...
                         'alpha',           [],...
                         'omega_z',         [],...
                         'X',               [],...
                         'Y',               [],...
                         'dA',              [],...
                         'isRolling',       true, ...
                         'roadProfile',     [], ...
                         'changeInRoadHeight', [] ...
                         );
    
    % % % Contact Width Parameters
    % % a1 = 1.03;
    % % a2 = 0.44;
    % % a3 = 195 * 0.35;
    % % 
    % % % Vertical Load
    % % Fz = 560 * 9.81;
    
    % Contact Width (y)
    % a = a1 * Fz^a2 + a3;
    a = 215 / 2;
    % Contact Length (x)
    b = 200 / 2;
    model_input.A = a * b;

    model_input.re = single( 0.5 * 15 * 25.4 + 0.55 * 195 ); % 195 / 55 R15
    model_input.alpha = single( deg2rad(0) );
   
    model_input.omega_z = single(0);
    
    % t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );
    t_sim = single( linspace(t_initial, t_final, t_final * fs_sim + 1) );
    if isRolling
    
        edge0 = 0.005 * t_final;
        edge1 = 0.0875 * t_final;
        edge2 = 0.8375 * t_final;
        edge3 = 0.92 * t_final;
        
        wheel_linear_vel = 8; %* 0.01 / (model_input.re/100);% 8 m/s; 28.8 km/h         %1.522137451171875e+02;
        rel_vel_between_road_and_wheel = wheel_linear_vel * 1;
        % data_range_2 = 0.1 / 3; % 10m/s; 36 km/h
        % data_range_3 = 0.2 / 3; % 20 m/s; 72 km/h
        smoothStep = smootherstep(edge0, edge1, t_sim) .* (1 - smootherstep(edge2, edge3, t_sim));

        v0(:, 1) = rel_vel_between_road_and_wheel * smoothStep(:);
        v0 = single(v0);

        ramp = rampFunc(edge1, edge2, t_sim);

        model_input.SR = ramp(:);%; -0.2 * ones(size(v0))
        model_input.SR = single(model_input.SR);
        omega =  single( v0 ./ ((model_input.SR + 1) * model_input.re ) );

        % Parameterise v0 and omega as functions of time ( used for
        % integration schemes)
        % % model_input.v0 = @(t) interp1(t_sim, v0, t, "linear");
        % % model_input.omega = @(t) interp1(t_sim, omega, t, "linear");
        model_input.v0 = griddedInterpolant(t_sim, v0, 'linear', 'none');
        model_input.omega = griddedInterpolant(t_sim, omega, 'linear', 'none');
        
        
    else
        edge0 = 0.005 * t_final;
        edge1 = 0.0875 * t_final;
        edge2 = 0.8375 * t_final;
        edge3 = 0.92 * t_final;
        
        data_range = 8e-3;% 8 m/s
        % data_range_2 = 0.1 / 3; % 10m/s; 36 km/h
        % data_range_3 = 0.2 / 3; % 20 m/s; 72 km/h
        v0(:, 1) = data_range * smootherstep(edge0, edge1, t_sim) .* (1 - smootherstep(edge2, edge3, t_sim));
        v0 = single(v0);
        % model_input.v0(:, 2) =  model_input.v0(:, 1) * 5;
        % model_input.v0(:, 3) =  model_input.v0(:, 2) * 2;

        omega = zeros(size(v0), 'single');
        model_input.SR = ones(size(model_input.v0), 'single');
        model_input.isRolling = false;

        model_input.v0 = griddedInterpolant(t_sim, v0, 'linear', 'none');
        model_input.omega = griddedInterpolant(t_sim, omega, 'linear', 'none');
        
    end
    
    model_input.numElems = uint16(numBrushes);
    model_input.LenTime_sim = single(fs_sim * t_final + 1);
    model_input.dt_sim = single( 1 / fs_sim );
    model_input.LenTime_save = single(fs_save * t_final + 1);
    model_input.dt_save = single( 1 / fs_save );
    model_input.dt_ratio = int32(model_input.dt_save / model_input.dt_sim);
    
    % Load pressure distribution
    % dataPath = fullfile('TM700 Pressure Distribution/', 'TM700Fz560Tr100r2_SubSampled_20x20.mat');
    % P_grid = load(dataPath);
    % model_input.P_grid = single(P_grid.P_grid_subsampled);
    model_input.Fz = max(smootherstep(edge0, edge1, t_sim) .* (1 - smootherstep(edge2, edge3, t_sim)) * 560 * 9.81,  100);
    
    % Calculate contact patch dimensions
    % ----------------------------------------------------------
    L0      = single(217.214915);
    Fz_0    = single(5716.443731);
    q_L     = single(0.467604);
    W0      = single(225.436018);
    q_W     = single(0.477134);
    % ----------------------------------------------------------
    a_max = 0.5 * L0 * (max(model_input.Fz) / Fz_0).^(q_L);
    b_max = 0.5 * W0 * (max(model_input.Fz) / Fz_0).^(q_W);
    x_vals = linspace(-b_max, b_max, numBrushes);
    y_vals = linspace(-a_max, a_max, numBrushes);
    dx = x_vals(2) - x_vals(1);
    dy = y_vals(2) - y_vals(1);
    dA = abs(dx * dy);

    [X, Y] = meshgrid(x_vals, y_vals);

    model_input.X = X;
    model_input.Y = Y;
    
    model_input.dA = dA;

    P_grid = calculatePressure(min(model_input.Fz), [], a_max, b_max, X, Y);
    model_input.P_grid = reshape(P_grid, numBrushes, numBrushes);
    
    % Load contact patch shape
    P_shape = load('C:\Users\coxde\OneDrive\Masters\BrushV2\Brush Sims\TM700 Pressure Distribution\TM700_ContactPatchShape.mat');
    P_shape = imresize(double(P_shape.binary_Array), [numBrushes, numBrushes], "lanczos2");
    P_shape = P_shape >= 0.02;
    model_input.P_shape = P_shape;
    
    % --- User-defined parameters ---
    spatial_sampling_frequency_x = 100; % Number of samples per meter
    spatial_sampling_frequency_y = 100; % Number of samples per meter
    
    dist_x = sum(model_input.v0(t_sim)) * model_input.dt_sim; % Total physical length of the grid (e.g., meters)
    dist_y = 2 * a / 1000;   % Total physical width of the grid (e.g., meters)
    
    num_points_x = dist_x * spatial_sampling_frequency_x; % Number of points in X-direction
    num_points_y = dist_y * spatial_sampling_frequency_y; % Number of points in Y-direction
    
    roadProfile = generateRoadProfile_3D(1e-5, ...
                                         spatial_sampling_frequency_x, spatial_sampling_frequency_y, ...
                                         dist_x, dist_y).';
    
    %%%%%%%%%%%%%%%%%%%%%% Feed road into vertical model %%%%%%%%%%%%%%%%%%%
    [Xdist, Ydist] = ndgrid(linspace(0, dist_x * 1000, num_points_x + 2 * numBrushes),...
                              linspace(0, dist_y * 1000, num_points_y));
    
    roadProfile = padarray(roadProfile, double([numBrushes, 0]),0, 'both'); % Add an array of zeros in front of road
    changeInRoadHeight = gradient(roadProfile) * spatial_sampling_frequency_x;



    model_input.roadProfile = griddedInterpolant(single(Xdist), single(Ydist), single(roadProfile), 'linear', 'none');
    model_input.changeInRoadHeight = griddedInterpolant(single(Xdist), single(Ydist), single(changeInRoadHeight), 'linear', 'none');
    
    model_input.roadProfile = griddedInterpolant(single(Xdist), single(Ydist), single(roadProfile), 'linear', 'none');
    model_input.changeInRoadHeight = griddedInterpolant(single(Xdist), single(Ydist), single(changeInRoadHeight), 'linear', 'none');



end

function ramp = rampFunc(start, stop, time)

    start_time = start; % The time when the ramp would normally start its smooth transition
    stop_time = stop; % The time when the ramp should reach its maximum value (1) and then flatten out
    % This one would go from 0 to 1 between edge1+1 and stop_time
    original_sharp_ramp_modified = zeros(size(time));
    idx_original_linear = (time >= start_time+1) & (time <= stop_time+1);
    % Make rampe
    original_sharp_ramp_modified(idx_original_linear) = (time(idx_original_linear) - (start_time+1)) / ((stop_time+1) - (start_time+1));
    % Set value after ramp to 1
    original_sharp_ramp_modified(time > stop_time+1) = 1;
    % Clamp values to between 0 and 1
    ramp = max(0, min(1, original_sharp_ramp_modified));
   
end

