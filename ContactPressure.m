function [Px, Pxy] = ContactPressure(Fz, a, b, x, n, lambda, xe, One_D, y)
%%%%%%%%%%%%%%%%%%%%%%% Help %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Description %%%%%%%%%%%%%%%%%%%%%
% Calculates the contact pressure distribution using two shape parameters
% nl. n and lambda. Can be extended to 2D by supplying n and lambda for the
% y direction. For this one needs to have half the footprint width as a
% function of x and can be passed in through the 'w_x' parameter.
%%%%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%
%   Fz      := Vertical Force
%   a, b    := relative coordinate
%   A, B    := Regulation parameters ensuring boundary conditions are satisfied 
%   n       := shape parameter 1
%   lambda  := shape parameter 2
%   xe      := Distance from geometric centre to centre of pressure
%   w_x     := Half the footprint as a function of x
%   x, y    := Meshgrid of x and y values if 2D. Otherwise 1D array
%%%%%%%%%%%%%%%%%%%%%%% Authors %%%%%%%%%%%%%%%%%%%%

%%%%%%%% For 2D case %%%%%%%%%%%%
% Meng Zhang, Hans-Joachim Unrau, Martin Gießler, Frank Gauterin,
% A detailed tire tread friction model considering dynamic friction states,
% Tribology International,
% Volume 193,
% 2024,
% 109342,
% ISSN 0301-679X,
% https://doi.org/10.1016/j.triboint.2024.109342.

%%%%%%%%%%%%% For 1D case %%%%%%%%%%%%%%%%
% K. Guo & D. Lu (2007) UniTire: unified tire model for vehicle dynamic
% simulation, Vehicle System Dynamics, 45:S1, 79-99, DOI: 10.1080/00423110701816742

%%%%%%%%%%%%%%%%%%%%%%%% Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (One_D)

        Px = f_z(Fz, x, a, n, lambda, xe) / 100;
        Pxy = [];

    else

        Px = [];
        Pxy = f_z_y(Fz, x, y, a, b, n, lambda, xe);
    end
end

function Px = f_z(Fz, x, a, n, lambda, xe)
    Px = Fz / (2 * a) .* g(x/a, n, lambda, xe, a);
end

function g = g(x, n, lambda, xe, a)
    A = (2 * n + 1) .* (4 * n +1) ./ (2 * n .* (4 * n + 1 + lambda));
    B = (3 * (2 * n + 3) .* (4 * n + 3) .* (4 * n + 1 + lambda) ./ ((2 * n + 1) .* (4 * n + 1) .* (4 * n + 3 + 3 * lambda))) .* (xe / a); 

    g = A .* (1 - (x).^(2*n)) .* (1 + lambda .* (x).^(2 * n)) .* (1 - B .* (x));
end

function w_x = w(x, a, b)
   w_x =  b;
end

function Pxy = f_z_y(Fz, x, y, a, b, n, lambda, xe)
    Px = f_z(Fz, x, a, n(1), lambda(1), xe(1));
    w_x = w([], [], b); % Calculate w_x as a function of x and y
    Pxy = Px ./ (2 .* w_x) .* g(y./w_x, n(2), lambda(2), xe(2), w_x);
end
