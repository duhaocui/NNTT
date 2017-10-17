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
mw = nefGaussianSumRV(1, w, 'parameters', 'wnefgaussianrv');
% measurement noise (glint)
epsilon = 0.1;
sigma_r = [0.2, 2];
sigma_theta = [0.015, 0.15];
R1 = diag([sigma_r(1)^2, sigma_theta(1)^2]);
R2 = diag([sigma_r(2)^2, sigma_theta(2)^2]);
vModel1 = nefGaussianRV(zeros(nMeas, 1), R1);
vModel2 = nefGaussianRV(zeros(nMeas, 1), R2);
mv = nefGaussianSumRV( (1 - epsilon), vModel1, epsilon, vModel2, 'parameters', 'wnefgaussianrv');
% initial condition
x0 = nefGaussianRV([232 * 1e3; 2290 * cosd(190); 88 * 1e3; 2290 * sind(190)], ...
    diag([1000^2, 20^2, 1000^2, 20^2]) );
mx0 = nefGaussianSumRV(1, x0, 'parameters', 'wnefgaussianrv');
% creating system and simulating its trajectory
system = nefEqSystem(f, h, mw, mv, mx0);
[z, x] = simulate(system, nSteps, []);
% setting up filters
disp('Setting up filters: GSM')
GSM = nefGSM(system, 'localFilter', 'ukf', ...
                     'pruningMethod', 'Threshold', ...
                     'percentOfMaxWeight', 0.1);
% filtering
disp('Runnign GSM filter ...')
t = cputime;
[val_GSM] = estimate(GSM, z, []);
GSMtime = cputime - t;
% evaluating means and MSEM
disp('Computing statistics of the obtained results ...')
t = T * (1:nSteps);
xest_GSM = zeros(nStates, nSteps);
msem_GSM = zeros(nStates, nSteps);
for i = 1:nSteps
    xest_GSM(:, i) = evalMean(val_GSM{i});
    msem_GSM(:, i) = (xest_GSM(:, i) - x(:, i) ).^2;
end
% show statistics
fprintf('Stats : MSEM\t\t time\n')
fprintf('GSM   : %f\t%f\n',mean(mean(msem_GSM) ), GSMtime);
% show plots
% testing
figure
subplot(3, 1, 1)
plot(x(1, :), x(3, :) )
subplot(3, 1, 2)
plot(t, z(1, :) )
subplot(3, 1, 3)
plot(t, z(2, :) )
figure
subplot(2, 2, 1)
plot(t, x(1, :), 'b', t, xest_GSM(1, :), 'r--')
ylabel('state')
xlabel('time')
legend('true', 'GSM')
subplot(2, 2, 2)
plot(t, msem_GSM(1, :) )
subplot(2, 2, 3)
plot(t, x(3, :), 'b', t, xest_GSM(3, :), 'r--')
ylabel('state')
xlabel('time')
legend('true', 'GSM')
subplot(2, 2, 4)
plot(t, msem_GSM(3, :) )