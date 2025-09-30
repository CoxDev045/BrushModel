function [P_grid] = calculatePressure(Fz, contact_shape, a, b, X, Y)
%CALCULATEPRESSURE Summary of this function goes here
    xe = [0, 0];
    lambda = [1, 1];
    n = [0.5, 0.5];

    [~, Pxy] = ContactPressure(Fz, a, b, X, n, lambda, xe, false, Y);
    P_grid = max(Pxy, 0);
    P_grid = P_grid(:);

    if ~isempty(contact_shape)
        P_grid = P_grid .* contact_shape;
    end

end

