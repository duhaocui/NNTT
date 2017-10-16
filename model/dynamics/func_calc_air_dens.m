function f = func_calc_air_dens(X, g, beta)

x = X(1);
dx = X(2);
y = X(3);
dy = X(4);

rho = func_calc_rho(y);

f = - g * rho * sqrt(dx^2 + dy^2) * [dx; dy] / (2 * beta);

