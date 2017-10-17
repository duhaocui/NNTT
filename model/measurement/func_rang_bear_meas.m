function  z = func_rang_bear_meas(x_curr, u, v, k)

x_R = 0;
y_R = 0;

x = x_curr(1);
y = x_curr(3);

z_r = sqrt( (x - x_R)^2 + (y - y_R)^2);
z_theta = atan2( (y - y_R), (x - x_R) );

z = zeros(2, 1);
z(1) = z_r;
z(2) = z_theta;

z = z + v;