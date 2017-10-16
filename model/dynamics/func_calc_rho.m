function [rho, c1, c2] = func_calc_rho(y)

if y < 9144
    c1 = 1.227;
    c2 = 1.093 * 1e-4;
else
    c1 = 1.754;
    c2 = 1.49 * 1e-4;
end

rho = c1 * exp(-c2 * y);