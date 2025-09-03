function X_next = evaluateRK4(func, dt, t, X)

    k1 = dt * func(t, X);
    k2 = dt * func(t, X + 0.5 * k1);
    k3 = dt * func(t, X + 0.5 * k2);
    k4 = dt * func(t, X + k3);
    
    X_next = X + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end