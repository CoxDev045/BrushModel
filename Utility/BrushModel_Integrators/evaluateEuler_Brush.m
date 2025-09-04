function [X_next, updated_obj] = evaluateEuler_Brush(func, dt, t, X_vec, obj)
    [k1, ~] = func(t, X_vec, obj);
    
    X_next = X_vec + dt * k1;

    [~, updated_obj] = func(t + dt, X_vec, obj);
end