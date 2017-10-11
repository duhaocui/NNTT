%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking a ship with unknown control
% State: position and velocity in the x and y directions
% Measuring bearing and range
% x(k+1) = F*x(k) + G*u(k) + w(k)
% z(k) = [atan(x2(k)/x1(k); sqrt(x1(k)^2 + x2(k)^2)] + v(k)
% nx = 4, nz = 2
% p(w(k)) = N{w(k):0,Q}
% p(v(k)) = N{w(k):0,R}
% p(x(0)) = N{x(0):m0,P0}
% for the estimation purpose the control is added as a new state component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% EKF - Extended Kalman filter
% UKF - Unscented Kalman filter
% DD1 - Divide difference filter 1st-order
% DD2 - Divide difference filter 2nd-order

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 7 (Tracking a ship with unknown control)\n')
fprintf('##################################################\n')

% #time steps
K = 455;

disp('Creating system and calculating its trajectory ...')
% STRUCTURAL DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.02;
dT = 72; % 0.02*3600
dVx0 = [4.8 -10.3 22.17 -17.2];
dVy0 = [8 1.5 -10.83 3.33];
th = [0 1 3.5 7.5 8.1 9.1];
tsec = th*3600;
tk = tsec/dT;
F = [1 0 dt 0; 
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];
G = [0 0; 
     0 0;
     1 0;
     0 1];
f = nefLinFunction(F,G,eye(4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = nefHandleFunction(@(x,u,v,t) [atan(x(2)/x(1)); sqrt(x(1)^2+x(2)^2)] + v,[4 0 2 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 1e-6*eye(4);
w = nefGaussianRV(zeros(4,1),Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 1e0*[4e-4*(pi/180)^2 0; 0 1e-4];
v = nefGaussianRV(zeros(2,1),R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = [20 50 0 -12]';
P0 = diag([1e1 1e1 1e1 1e1]);
x0 = nefGaussianRV(m0,P0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = nefEqSystem(f,h,w,v,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(2,K);
u(:,tk(2:5)) = [dVx0;dVy0];
[z,x] = simulate(system,K,u);

% MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm = [F G;zeros(2,4) eye(2)];
fm = nefLinFunction(Fm,[],eye(6));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dNoise_handle = @(x,u,v,k) eye(2);
dState_handle = @(x,u,v,k) [1/(1+(x(2)/x(1))^2)*(-x(2)/x(1)^2), 1/(1+(x(2)/x(1))^2)*(1/x(1)), 0, 0, 0, 0;...
                            x(1)/sqrt(x(1)^2+x(2)^2),           x(2)/sqrt(x(1)^2+x(2)^2),     0, 0, 0, 0];
hm = nefHandleFunction(@(x,u,v,t) [atan(x(2)/x(1)); sqrt(x(1)^2+x(2)^2)] + v,[6 0 2 0],'Diff1Noise',dNoise_handle,'Diff1State',dState_handle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qm = [Q zeros(4,2); zeros(2,4) 1e-3*eye(2)];
wm = nefGaussianRV(zeros(6,1),Qm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0m = [m0; 0; 0]; 
P0m = [P0 zeros(4,2); zeros(2,4) 1e-1*eye(2)];
x0m = nefGaussianRV(m0m,P0m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = nefEqSystem(fm,hm,wm,v,x0m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: KALMAN, UKF, DD1, DD2 ...')
EKF = nefKalman(model);
UKF = nefUKF(model);
DD1 = nefDD1(model);
DD2 = nefDD2(model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running EKF filter ...')
t = cputime;
[val_EKF] = estimate(EKF,z,u);
EKFtime = cputime-t;
disp('Running UKF filter ...')
t = cputime;
[val_UKF] = estimate(UKF,z,u);
UKFtime = cputime-t;
disp('Running DD1 filter ...')
t = cputime;
[val_DD1] = estimate(DD1,z,u);
DD1time = cputime-t;
disp('Running DD2 filter ...')
t = cputime;
[val_DD2] = estimate(DD2,z,u);
DD2time = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = t
  xest_EKF(:,i) = evalMean(val_EKF{i});
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_DD1(:,i) = evalMean(val_DD1{i});
  xest_DD2(:,i) = evalMean(val_DD2{i});

  msem_EKF(:,i) = (xest_EKF(1:4,i)-x(:,i)).^2; 
  msem_UKF(:,i) = (xest_UKF(1:4,i)-x(:,i)).^2; 
  msem_DD1(:,i) = (xest_DD1(1:4,i)-x(:,i)).^2; 
  msem_DD2(:,i) = (xest_DD2(1:4,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('EKF   : %f\t%f\n',mean(mean(msem_EKF)),EKFtime);
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('DD1   : %f\t%f\n',mean(mean(msem_DD1)),DD1time);
fprintf('DD2   : %f\t%f\n',mean(mean(msem_DD2)),DD2time);


subplot(4,1,1)
plot(t,x(1,:),'b',t,xest_DD1(1,:),'r*',t,xest_DD2(1,:),'go',t,xest_UKF(1,:),'k+',t,xest_EKF(1,:),'r+')
ylabel('x - position')
xlabel('time')
legend('true','DD1','DD2','UKF','EKF')

subplot(4,1,2)
plot(t,x(2,:),'b',t,xest_DD1(2,:),'r*',t,xest_DD2(2,:),'go',t,xest_UKF(2,:),'k+',t,xest_EKF(2,:),'r+')
ylabel('y - position')
xlabel('time')
legend('true','DD1','DD2','UKF','EKF')

subplot(4,1,3)
plot(t,x(3,:),'b',t,xest_DD1(3,:),'r*',t,xest_DD2(3,:),'go',t,xest_UKF(3,:),'k+',t,xest_EKF(3,:),'r+')
ylabel('x - velocity')
xlabel('time')
legend('true','DD1','DD2','UKF','EKF')

subplot(4,1,4)
plot(t,x(4,:),'b',t,xest_DD1(4,:),'r*',t,xest_DD2(4,:),'go',t,xest_UKF(4,:),'k+',t,xest_EKF(4,:),'r+')
ylabel('y - velocity')
xlabel('time')
legend('true','DD1','DD2','UKF','EKF')
