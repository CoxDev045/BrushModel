function [model_input] = brush_init(numBrushes, fs_sim, fs_save, t_initial, t_final)
    
    % Initialise model input struct
    model_input = struct('numElems',[],...
                         'LenTime_sim',[], ...
                         'LenTime_save',[], ...
                         'dt_sim',[],...
                         'dt_save',[],...
                         'dt_ratio',[],...
                         'Press',[], ...
                         'omega',[], ...
                         'SR',[],...
                         're',[],...
                         'alpha',[],...
                         'omega_z',[],...
                         'X',[],...
                         'Y',[]);
    
    % Contact Width Parameters
    a1 = 1.03;
    a2 = 0.44;
    a3 = 195 * 0.35;
    
    % Vertical Load
    Fz = 560 * 9.81;%3000;
    
    % Contact Width (y)
    a = a1 * Fz^a2 + a3;
    % Contact Length (x)
    b = 195 / 2;
    
    model_input.re = 0.5 * 15 * 25.4 + 0.55 * 195; % 195 / 55 R15
    model_input.alpha = deg2rad(0);
   
    model_input.omega_z = 0;
    
    t_save = linspace(t_initial, t_final, t_final * fs_save + 1);
    
    edge0 = 0.5;
    edge1 = 10.5;%25.5;
    edge2 = 61.5;%42;
    edge3 = 81.5;%67;
    
    data_range = 0.01 / (model_input.re/100);% 1 m/s; 3.6 km/h         %1.522137451171875e+02;
    % data_range_2 = 0.1 / 3; % 10m/s; 36 km/h
    % data_range_3 = 0.2 / 3; % 20 m/s; 72 km/h
    model_input.omega(:, 1) = data_range * smootherstep(edge0, edge1, t_save) .* (1 - smootherstep(edge2, edge3, t_save));
    model_input.omega(:, 2) =  model_input.omega(:, 1) * 5;
    model_input.omega(:, 3) =  model_input.omega(:, 2) * 2;
    model_input.SR = 0.1;
    model_input.v0 =  model_input.omega * model_input.re / (model_input.SR + 1);
    
    x_vals = linspace(-b, b, numBrushes);
    y_vals = linspace(-a, a, numBrushes);
    dx = x_vals(2) - x_vals(1);
    dy = y_vals(2) - y_vals(1);
    dA = dx * dy;
    
    [X, Y] = meshgrid(x_vals, y_vals);
    
    model_input.X = X;
    model_input.Y = Y;

    model_input.numElems = uint16(numBrushes);
    model_input.LenTime_sim = single(fs_sim * t_final + 1);
    model_input.dt_sim = 1 / fs_sim;
    model_input.LenTime_save = single(fs_save * t_final + 1);
    model_input.dt_save = 1 / fs_save;
    model_input.dt_ratio = int32(model_input.dt_save / model_input.dt_sim);
    
    % load('Sliding_Pressure_fs_sim_1e5_fs_save_1e3_len_3s.mat')
    P_grid = load("C:\Users\coxde\OneDrive\Masters\BrushV2\TM700 Pressure Distribution\TM700Fz560Tr0_Buff_r1_SubSampled_20x20.mat");
    model_input.Press = P_grid.P_grid_subsampled;


end

