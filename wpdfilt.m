function [r] = wpdfilt(min, max, x)
    r = x(x >= min & x <= max);
end