classdef ContactPatch %#codegen
    properties    
        contact_area            (1, 1) double = 1;
        numBrushes              (1, 1) double = 10;
        contact_area_per_brush  (1, 1) double = 1;
        re                      (1, 1) double = 30;
        alpha                   (1, 1) double = 0;
        fs                      (1, 1) double = 1e3; % Hz
        dt                      (1, 1) double = 1e-3;
        omega_z                 (1, 1) double = 0;

        % Pre-allocate array of empty objects
        brushArray (numBrushes, numBrushes) BrushV2 = BrushV2();
    end

    methods
        function obj = ContactPatch(numBrushes, contact_area, fs, re, omega_z)
            arguments
                numBrushes      (1, 1) double
                contact_area    (1, 1) double
                fs              (1, 1) double
                re              (1, 1) double
                omega_z         (1, 1) double
            end
            
            obj.numBrushes = numBrushes;
            obj.contact_area = contact_area;
            obj.contact_area_per_brush = contact_area / (numBrushes * numBrushes);
            obj.fs = fs;
            obj.re = re;
            obj.dt = 1/fs;
            obj.omega_z = omega_z;

        end

        function [time_span, omega, v0] = createSimInputs(t_initial, t_final, slip_ratio)
            time_span = linspace(t_initial, t_final, t_final * obj.fs);
            omega = linspace(0, t_final, length(time_span));
            v0 = omega * obj.re / (slip_ratio + 1);
        end
    end

    % Initialize time and force arrays
    time = zeros(length(time_span), 1);
    force = zeros(length(time_span), 1);

    x_vals = linspace(-a, a, numBrushes);
    y_vals = linspace(-b, b, numBrushes);

    [X, Y] = meshgrid(x_vals, y_vals);

    % Initialise PRessure function coefficients
    nx = 2;
    ny = 1;
    n_params = [nx, ny];
    
    lambda_x = 1;
    lambda_y = 1;
    lambda = [lambda_x, lambda_y];
    
    xe = 0; % mm
    ye = 0; % mm
    COP = [xe, ye];
    
    One_D = false;
    
    % Calculate 2D pressure distribution
    [~, P_grid] = ContactPressure(Fz, a, b, X, n_params, lambda, COP, One_D, Y);
    P_grid = real(P_grid);
    
    sim_solution = zeros(numBrushes, numBrushes, length(time_span), 9);
    % sim_solution(:, :, :, 1) = repmat(P_grid, [1, 1, length(time_span)]);
    
    % Calculate the amount to shift per timestep
    shift_amount = (omega * dt * re);  % Adjust re for scale if needed
    shift_amount_cumulative = cumsum(shift_amount);

    shift_amount_cumulative = floor(shift_amount_cumulative);


    for i = 1:length(time_span)
        % Shift the pressure distribution downwards and wrap it
        sim_solution(:, :, i, 1) = circshift(P_grid, [-shift_amount_cumulative(i), 0]);
    end


    %
    for i = 1:numBrushes
        for j = 1:numBrushes
            brushArray(i, j).x = x_vals(i); 
            brushArray(i, j).y = y_vals(j);
            brushArray(i, j).press = P_grid(i, j);
            brushArray(i, j).ky = brushArray(i, j).ky_0 + brushArray(i, j).ky_l * brushArray(i, j).press;
            brushArray(i, j).kx = brushArray(i, j).phi + brushArray(i, j).ky;
            brushArray(i, j).vs = 0;
        end
    end
    
    % Initialize each object with specific values
    for i = 1:length(time_span)
        for j = 1:numBrushes
            for k = 1:numBrushes
                brush = brushArray(j, j);
                brush.press = P_grid(j, k);
                brush.ky = brush.ky_0 + brush.ky_l * brush.press;
                brushkx = brush.phi + brush.ky;
                brush = brush.update_brush(omega(i), omega_z, re, v0(i), alpha, dt);
        
                % Store results from the brush object
                % sim_solution(j, k, i, 2) = "Total Displacement [mm]";
                % sim_solution(j, k, i, 3) = "Total Shear Stress [MPa]";
                % if isnan(brushArray(j, k).vs)
                %     fprintf('For brush located at (%.4f, %.4f), vs was calculated to be NAN at %.4fs. \n', x_vals(i), y_vals(j), time_span(i));
                % end
                if ~isnan(brushArray(j, k).delta_x)
                    sim_solution(j, k, i, 4) = brushArray(j, k).delta_x;
                else
                    brushArray(j, k).delta_x = 0;
                    sim_solution(j, k, i, 4) = brushArray(j, k).delta_x;
                end
                if ~isnan(brushArray(j, k).delta_y)
                    sim_solution(j, k, i, 5) = brushArray(j, k).delta_y;
                else
                    brushArray(j, k).delta_y = 0;
                    sim_solution(j, k, i, 5) = brushArray(j, k).delta_y;
                end
               if ~isnan(brushArray(j, k).tauX)
                    sim_solution(j, k, i, 6) = brushArray(j, k).tauX;
                else
                    brushArray(j, k).tauX = 0;
                    sim_solution(j, k, i, 6) = brushArray(j, k).tauX;
                end
               if ~isnan(brushArray(j, k).tauY)
                    sim_solution(j, k, i, 7) = brushArray(j, k).tauY;
                else
                    brushArray(j, k).tauY = 0;
                    sim_solution(j, k, i, 7) = brushArray(j, k).tauY;
                end
     
                sim_solution(j, k, i, 8) = brushArray(j, k).slide;
               
                if ~isnan(brushArray(j, k).mu)
                    sim_solution(j, k, i, 9) = brushArray(j, k).mu;
                else
                    brushArray(j, k).mu = 0;
                    sim_solution(j, k, i, 9) = brushArray(j, k).mu;
                end
    
                if brushArray(j, k).x < -a
                    brushArray(j, k).x = a;
                end
            end
        end
        % Compute total displacement over grid
        sim_solution(:, :, i, 2) = hypot(sim_solution(:, :, i, 4), sim_solution(:, :, i, 5));
        % Compute total Stress over grid
        sim_solution(:, :, i, 3) = hypot(sim_solution(:, :, i, 6), sim_solution(:, :, i, 7));
        % Intergrate stress to get force
        force(i) = trapz(trapz(sim_solution(:, :, i, 6))) * contact_area_per_brush;
        if mod(i, 10000) == 0
            fprintf("%d iterations completed! %d remaining.\n", i, length(time_span) - i);
        end
    end
end