function y = smoothstep(edge0, edge1, x)
   y = clamp((x - edge0) ./ (edge1 - edge0));
   y = y.^2 .* (3.0 - 2.0 * y); 
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