function tempPress = shiftPressure(X, Y, P, shift_amount, max_val, min_val, Pmin)

    % Corrects for the tendency of the mod function to shift the leading edge to the centre of the contact patch
    X_shift_corrector = mod(X, max_val - abs(min_val));

    % Calculate shift amount due to tyre rolling
    X_shifted = mod(X_shift_corrector + shift_amount, max_val - (min_val)) + min_val;
    
    % Linearly interpolate pressure along grid
    tempPress = interp2(X, Y, P, X_shifted, Y, 'linear');

    % Mask Pressure grid to remove negative or small values (done to smooth out final outputs)
    tempPress = max(tempPress(:), Pmin);
end