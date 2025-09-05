function [X_next, updated_obj] = evaluateRK4_Brush(func, dt, t, X, slide_obj, args)

    [k1,obj_k1] = func(t, X, slide_obj, args);
    [k2,obj_k2] = func(t + 0.5 *dt, X + 0.5 * dt * k1, obj_k1, args);
    [k3,obj_k3] = func(t + 0.5 * dt, X + 0.5 * dt * k2, obj_k2, args);
    [k4,obj_k4] = func(t + dt, X + dt * k3, obj_k3, args);
    % Calculate the next state vector
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    % Update slide_obj with the final, correct state
    % Should call the calculateStresses method one last time
    % % [~, updated_obj] = func(t + dt, X_next, slide_obj);
    updated_obj = obj_k4;
end