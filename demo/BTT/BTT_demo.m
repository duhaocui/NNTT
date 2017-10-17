clear
clc
close all

% time steps
T = 2;
nSteps = 120 / T;
disp('Creating system and calculating its trajectory ...')
% function f in state equation
nStates = 4;
f = nefHandleFunction(@func_BTT_dyn, [nStates, 0, 4, 0], 'isAdditive', 1);
% function h in measurement equation
nMeas = 2;
h = nefHandleFunction(@func_rang_bear_meas, [nStates, 0, 2, 0], 'isAdditive', 1);
% state noise
q = 1;
Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
Q = q * blkdiag(Theta, Theta);
w = nefGaussianRV(zeros(nStates, 1), Q);
% measurement noise
sigma_r = 100;
sigma_theta = 0.05;
R = diag([sigma_r^2, sigma_theta^2]);
v = nefGaussianRV(zeros(nMeas, 1), R);
% initial condition
x0 = nefGaussianRV([232 * 1e3; 2290 * cosd(190); 88 * 1e3; 2290 * sind(190)], ...
    diag([1000^2, 20^2, 1000^2, 20^2]) );
% creating system and simulating its trajectory
system = nefEqSystem(f, h, w, v, x0);
[z, x] = simulate(system, nSteps, []);
% testing
figure
subplot(3, 1, 1)
plot(x(1, :), x(3, :) )
subplot(3, 1, 2)
plot(1:nSteps, z(1, :) )
subplot(3, 1, 3)
plot(1:nSteps, z(2, :) )
% setting up filters
disp('Setting up filters: UKF')
UKF = nefUKF(system);
% filterign
disp('Running UKF filter ...')
t = cputime;
[val_UKF] = estimate(UKF, z, []);
UKFtime = cputime - t;
% evaluating means and MSEM
disp('Computing statistics of the obtained results ...')
t = T * (1:nSteps);
xest_UKF = zeros(nStates, nSteps);
msem_UKF = zeros(nStates, nSteps);
for i = 1:nSteps
    xest_UKF(:, i) = evalMean(val_UKF{i});
    msem_UKF(:, i) = (xest_UKF(:, i) - x(:, i) ).^2;
end
% show statistics
fprintf('Stats : MSEM\t\t time\n')
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF) ), UKFtime);
% show plots
figure
subplot(2, 2, 1)
plot(t, x(1, :), 'b', t, xest_UKF(1, :), 'r--')
ylabel('state')
xlabel('time')
legend('true', 'UKF')
subplot(2, 2, 2)
plot(t, msem_UKF(1, :) )
subplot(2, 2, 3)
plot(t, x(3, :), 'b', t, xest_UKF(3, :), 'r--')
ylabel('state')
xlabel('time')
legend('true', 'UKF')
subplot(2, 2, 4)
plot(t, msem_UKF(3, :) )