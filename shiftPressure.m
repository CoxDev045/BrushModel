function [P_grid_shifted, P_grid_save] = shiftPressure(numBrushes, max_numElems, max_LenSimTime, dt_ratio, P_grid, X, Y, omega, dt_sim, re)

    arguments (Input)
        numBrushes      (1, 1) uint16  % **Must be a constant**
        max_numElems    (1, 1) uint16  % **Must be a constant**
        max_LenSimTime  (1, 1) single  % **Must be a constant**
        dt_ratio        (1, 1) single  % **Must be a constant**
        P_grid          (20, 20) single
        X               (20, 20) single
        Y               (20, 20) single
        omega           (:, 1) single
        dt_sim          (1, 1) single
        re              (1, 1) single
    end

    Press = P_grid;%reshape(P_grid, numBrushes, numBrushes);
    tempX = X;%reshape(X, numBrushes, numBrushes);
    tempY = Y;%reshape(Y, numBrushes, numBrushes);
    maxX = max(X, [], 'all');
    maxY = max(Y, [], 'all');
    minX = min(X, [], 'all');
    minY = min(Y, [], 'all');

    progress_steps = round(linspace(1, max_LenSimTime, 10)); % 10 checkpoints

    shift_amount = cumsum(omega * dt_sim * re);
    % % shift_mask = (gradient(floor(shift_amount)) > 0);

    P_grid_shifted = zeros(max_numElems, max_LenSimTime, 'single');
    P_grid_save = zeros(max_numElems, ceil(max_LenSimTime / dt_ratio), 'single');

    j = 1;
    
    for i = 1:max_LenSimTime
        X_shifted = mod(tempX + shift_amount(i), maxX - minX) + minX;
        % Y_shifted = mod(tempY - shift_amount_cumulative(i) - 1, maxY - minY) + minY;

        % % tempPress = circshift(tempPress, [0, shift_amount(i)]);

        tempPress = interp2(tempX, tempY, Press, X_shifted, tempY, 'linear', 0);

        % % tempPress = tempPress .* perlin(size(tempPress));
        if mod(i, dt_ratio) == 0
            P_grid_save(:, j) = tempPress(:);
            
            j = j + 1;

            % Shift pressure distribution
            P_grid_shifted(:, i) = tempPress(:);
        else
            % Shift pressure distribution
            P_grid_shifted(:, i) = tempPress(:);
        end

        %%%%%%%%%%%%%% Update progress at defined steps %%%%%%%%%%%%%%%%%%%
        if any(i == progress_steps)
            fprintf('%.2f%% completed (%.2f/%.2f steps)\n', round(100 * i / max_LenSimTime), i, max_LenSimTime);
        end
    end
end

% % function P_grid_shifted = shiftPressure(P_grid, X, Y, omega, time_span, dt, re)
% %     arguments (Input)
% %         P_grid (:, 1) single
% %         X (:, 1) single
% %         Y (:, 1) single
% %         omega (:, 1) single
% %         time_span (:, 1) single
% %         dt (1, 1) single
% %         re (1, 1) single
% %     end
% % 
% %     numElems = sqrt(size(P_grid, 1));
% %     tempPress = reshape(P_grid, numElems, numElems);
% %     tempX = reshape(X, numElems, numElems);
% %     tempY = reshape(Y, numElems, numElems);
% % 
% %     % Get grid properties
% %     maxX = max(X);
% %     maxY = max(Y);
% %     minX = min(X);
% %     minY = min(Y);
% % 
% %     % Calculate grid spacing
% %     dx = (maxX - minX) / (numElems - 1);
% %     dy = (maxY - minY) / (numElems - 1);
% % 
% %     % Calculate cumulative shift
% %     shift_amount_cumulative = cumsum(omega * dt * re);
% % 
% %     % Prepare output
% %     P_grid_shifted = zeros(size(P_grid, 1), length(time_span), 'single');
% % 
% %     % Original grid coordinates
% %     [Xq, Yq] = meshgrid(minX:dx:maxX, minY:dy:maxY);
% % 
% %     for i = 1:length(time_span)
% %         % Calculate shifted coordinates with proper wrapping
% %         % X_shifted = mod(tempX - shift_amount_cumulative(i), maxX - minX) + minX;
% %         Y_shifted = mod(tempY - shift_amount_cumulative(i), maxY - minY) + minY;
% % 
% %         % Interpolate the pressure at shifted coordinates
% %         P_shifted = interp2(Xq, Yq, tempPress, Xq, Y_shifted, 'spline', 0);
% % 
% %         % Store the shifted pressure
% %         P_grid_shifted(:, i) = P_shifted(:);
% %     end
% % end