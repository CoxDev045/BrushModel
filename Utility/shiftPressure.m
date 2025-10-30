function tempPress = shiftPressure(XStatic, XCurrent, Y, P, shift_amount, max_val, min_val)

    % Calculate shift amount due to tyre rolling
    X_shifted = XCurrent + shift_amount;
    grid_width = max_val + abs(min_val); % a + |-a| = 2 * a 
    % X_shifted = mod(shifted_values - max_val, grid_width) + max_val;
    % Create a boolean mask for values that have exceeded the bounds
    is_too_high = X_shifted > max_val;
    is_too_low = X_shifted < min_val;
    
    % Apply the manual wrap using the masks
    X_shifted(is_too_high) = X_shifted(is_too_high) - grid_width;
    X_shifted(is_too_low) = X_shifted(is_too_low) + grid_width;
        
    % Linearly interpolate pressure along grid
    tempPress = interp2(XStatic, Y, P, X_shifted, Y, "linear");


    % Mask Pressure grid to remove negative or small values (done to smooth out final outputs)
    tempPress = max(tempPress(:), 0);
end