function [X_next, updated_obj] = evaluateRK4_Brush(func, dt, t, X, slide_obj, args)

    [k1,~] = func(t, X, slide_obj, args);
    [k2,~] = func(t, X + 0.5 * dt * k1, slide_obj, args);
    [k3,~] = func(t, X + 0.5 * dt * k2, slide_obj, args);
    [k4,~] = func(t, X + dt * k3, slide_obj, args);
    % Calculate the next state vector
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    % Update slide_obj with the final, correct state
    % Should call the calculateStresses method one last time
    [~, updated_obj] = func(t + dt, X_next, slide_obj, args);
end