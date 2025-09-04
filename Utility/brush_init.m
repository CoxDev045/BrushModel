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
                         'isRolling',       true ...
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
    a = 220.8750 / 2;
    % Contact Length (x)
    b = 195 / 2;
    model_input.A = a * b;

    model_input.re = single( 0.5 * 15 * 25.4 + 0.55 * 195 ); % 195 / 55 R15
    model_input.alpha = single( deg2rad(0) );
   
    model_input.omega_z = single(0);
    
    % t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );
    t_sim = single( linspace(t_initial, t_final, t_final * fs_sim + 1) );
    if isRolling
    
        edge0 = 0.5;
        edge1 = 10.5;%25.5;
        edge2 = 100.5;%42;
        edge3 = 110.5;%67;
        
        wheel_linear_vel = 8; %* 0.01 / (model_input.re/100);% 8 m/s; 28.8 km/h         %1.522137451171875e+02;
        rel_vel_between_road_and_wheel = wheel_linear_vel * 1;
        % data_range_2 = 0.1 / 3; % 10m/s; 36 km/h
        % data_range_3 = 0.2 / 3; % 20 m/s; 72 km/h
        smoothStep = smootherstep(edge0, edge1, t_sim) .* (1 - smootherstep(edge2, edge3, t_sim));

        v0(:, 1) = rel_vel_between_road_and_wheel * smoothStep(:);
        v0 = single(v0);
        % model_input.omega(:, 2) =  model_input.omega(:, 1) * 5;
        % model_input.omega(:, 3) =  model_input.omega(:, 2) * 2; linspace(0, 1, length(model_input.omega)).' .* 
        
        ramp = rampFunc(edge1, edge2, t_sim);

        model_input.SR = ramp(:);
        model_input.SR = single(model_input.SR);
        omega =  single( v0 ./ ((model_input.SR + 1) * model_input.re ) );

        % Parameterise v0 and omega as functions of time ( used for
        % integration schemes)
        model_input.v0 = @(t) interp1(t_sim, v0, t, "linear");
        model_input.omega = @(t) interp1(t_sim, omega, t, "linear");
        
    else
        edge0 = 0.5;
        edge1 = 10.5;%25.5;
        edge2 = 100.5;%42;
        edge3 = 110.5;%67;
        
        data_range = 8;% 8 m/s
        % data_range_2 = 0.1 / 3; % 10m/s; 36 km/h
        % data_range_3 = 0.2 / 3; % 20 m/s; 72 km/h
        model_input.v0(:, 1) = data_range * smootherstep(edge0, edge1, t_sim) .* (1 - smootherstep(edge2, edge3, t_sim));
        model_input.v0 = single(model_input.v0);
        % model_input.v0(:, 2) =  model_input.v0(:, 1) * 5;
        % model_input.v0(:, 3) =  model_input.v0(:, 2) * 2;

        model_input.omega = zeros(size(model_input.v0), 'single');
        model_input.SR = ones(size(model_input.v0), 'single');
        model_input.isRolling = false;
    end
    
    x_vals = linspace(-b, b, numBrushes);
    y_vals = linspace(-a, a, numBrushes);
    dx = x_vals(2) - x_vals(1);
    dy = y_vals(2) - y_vals(1);
    dA = dx * dy;
    
    [X, Y] = meshgrid(x_vals, y_vals);
    
    model_input.X = single(X);
    model_input.Y = single(Y);
    model_input.dA = dA;

    model_input.numElems = uint16(numBrushes);
    model_input.LenTime_sim = single(fs_sim * t_final + 1);
    model_input.dt_sim = single( 1 / fs_sim );
    model_input.LenTime_save = single(fs_save * t_final + 1);
    model_input.dt_save = single( 1 / fs_save );
    model_input.dt_ratio = int32(model_input.dt_save / model_input.dt_sim);
    
    % load('Sliding_Pressure_fs_sim_1e5_fs_save_1e3_len_3s.mat')
    dataPath = fullfile('TM700 Pressure Distribution/', 'TM700Fz560Tr100r2_SubSampled_20x20.mat');
    P_grid = load(dataPath);
    model_input.P_grid = single(P_grid.P_grid_subsampled);


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

