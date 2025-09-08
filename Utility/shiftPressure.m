function tempPress = shiftPressure(X, Y, P, shift_amount, max_val, min_val)

    % Corrects for the tendency of the mod function to shift the leading edge to the centre of the contact patch
    % X_shift_corrector = mod(X, max_val - abs(min_val));

    % Calculate shift amount due to tyre rolling
    X_shifted = X + shift_amount;
    grid_width = max_val + abs(min_val); % a + |-a| = 2 * a 
    % X_shifted = mod(shifted_values - max_val, grid_width) + max_val;
    % Create a boolean mask for values that have exceeded the bounds
    is_too_high = X_shifted > max_val;
    is_too_low = X_shifted < min_val;
    
    % Apply the manual wrap using the masks
    X_shifted(is_too_high) = X_shifted(is_too_high) - grid_width;
    X_shifted(is_too_low) = X_shifted(is_too_low) + grid_width;
        
    % Linearly interpolate pressure along grid
    tempPress = interp2(X, Y, P, X_shifted, Y, "makima");

    % % % Define the width of the smooth transition zones
    % % ramp_dist = 0.2 * (grid_width);
    % % 
    % % % Define edges for smootherstep function
    % % % Ramp-up: from leading edge to just inside
    % % edge0 = min_val;        % Furthest edge of contact patch in the negative direction
    % % edge1 = min_val + ramp_dist; % 10% inboard of the trailing edge
    % % 
    % % % Ramp-down: from just inside to the trailing edge
    % % edge2 = max_val - ramp_dist; % 10% inboard of the leading edge
    % % edge3 = max_val;        % Furthest edge of contact patch in the positive direction
    % % 
    % % % Calculate the smooth ramp factor based on the new shifted x-values
    % % smooth_in = smootherstep(edge0, edge1, X_shifted);
    % % smooth_out = 1 - smootherstep(edge2, edge3, X_shifted);
    % % ramp_factor = smooth_in .* smooth_out;
    % % 
    % % % Apply the ramp factor to the interpolated pressure
    % % tempPress = tempPress .* ramp_factor;

    % Mask Pressure grid to remove negative or small values (done to smooth out final outputs)
    tempPress = max(tempPress(:), 0);
end