function X_next = evaluateRKF5(func, dt, t, X_vec)

    k1 = dt * func(t, X_vec);
    k2 = dt * func(t, X_vec + 0.25 *  k1);
    k3 = dt * func(t, X_vec + 3/32 *  k1 + 9/32 *  k2);
    k4 = dt * func(t, X_vec + 1932/2197 *  k1 - 7200/2197 *  k2 + 7296/2197 *  k3);
    k5 = dt * func(t, X_vec + 439/216 *  k1 - 8 *  k2 + 3680/513 *  k3 - 845/4104 *  k4);
    k6 = dt * func(t, X_vec - 8/27 * k1 + 2 *  k2 - 3544/2565 *  k3 + 1859/4104 *  k4 - 11/40 *  k5);
    
    X_next =  X_vec + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4 - (9/50) * k5 + (2/55) * k6;
end