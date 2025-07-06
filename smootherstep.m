function y = smootherstep(edge0, edge1, x)
   y = clamp((x - edge0) ./ (edge1 - edge0));
   y = y.^3 .* (y .* (6.0 * y - 15.0) + 10.0); 
end


function y = clamp(x, lowerlimit, upperlimit)
    if (nargin < 2)
        lowerlimit = 0;
    end
    if (nargin < 3)
        upperlimit = 1;
    end
    
    y = x;
    lower_ind = x < lowerlimit;
    upper_ind = x > upperlimit;

    y(lower_ind) = lowerlimit;
    y(upper_ind) = upperlimit;
    
end