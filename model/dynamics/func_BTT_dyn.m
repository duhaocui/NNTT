function x_curr = func_BTT_dyn(x_prev, u, w, k)

T = 2;
g = 9.81;
beta_tgt = 4 * 1e4;
F = [1, T, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, T;
    0, 0, 0, 1];
G = [T^2 / 2, 0
    T,       0;
    0,       T^2 / 2;
    0,       T];

f = func_calc_air_dens(x_prev, g, beta_tgt);

x_curr = F * x_prev + G * f + G * [0; -g] + w;