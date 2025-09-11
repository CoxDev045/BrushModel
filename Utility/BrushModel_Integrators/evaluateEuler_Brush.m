function X_next = evaluateEuler_Brush(func, dt, t, X_vec, obj, args)
    k1= func(t, X_vec, obj, args);
    
    X_next = X_vec + dt * k1;
end