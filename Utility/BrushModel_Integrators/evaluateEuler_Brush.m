function [X_next, updated_obj] = evaluateEuler_Brush(func, dt, t, X_vec, obj, args)
    [k1, ~] = func(t, X_vec, obj, args);
    
    X_next = X_vec + dt * k1;

    [~, updated_obj] = func(t + dt, X_vec, obj, args);
end