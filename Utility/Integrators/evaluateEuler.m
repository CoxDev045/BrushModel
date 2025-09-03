function X_next = evaluateEuler(func, dt, t, X_vec)
    k1 = dt * func(t, X_vec);
    
    X_next = X_vec + k1;
end