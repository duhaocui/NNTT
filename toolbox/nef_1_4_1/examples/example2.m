%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predator-prey model
%
%      | T*(1-a)*x1 + T*b*x1*x2 + w1 |
%  f = |                             |
%      | T*(1+c)*x2 + T*d*x1*x2 + w2 |
%
%           | x1 |                  | w1 |
% where x = |    | is state and w = |    | noise
%           | x2 |                  | w2 |
%
% Constants are given as follows:
%  T = 1
%  a = 0.04
%  b = c = d = 0.08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% UKF - Unscented Kalman filter
% KALMAN - Kalman filter
% PF - Particle filter
% RTSS - Rauch Tung Striebel smoother
% URTSS - Unscented Rauch Tung Striebel smoother

% NEF version 1.2.0
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all

fprintf('##################################################\n')
fprintf('# EXAMPLE 2 (Predator - prey model)\n')
fprintf('##################################################\n')

% #time steps
K = 20;

disp('Creating system and calculating its trajectory ...')
% STRUCTURAL DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = nefHandleFunction(@(x,u,w,k) [0.96*x(1)+0.08*x(1)*x(2)+w(1);1.08*x(2)-0.08*x(1)*x(2)+w(2)],[2 0 2 0],...
  'diff1Noise',@(x,u,v,k) eye(2),...
  'diff1State',@(x,u,w,k) [0.96+0.08*x(2) 0.08*x(1);-0.08*x(2) 1.08-0.08*x(1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [1 0];
h = nefLinFunction(H,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = eye(2)*0.005;
w = nefGaussianRV([0 0]',Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 0.01;
v = nefGaussianRV(0,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([0.9;-0.85],1e-3*eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system and simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);
[z,x] = simulate(system,K,[]);
% PROBABILISTIC DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmean =  nefHandleFunction(@(x,u,w,k) [0.96*x(1)+0.08*x(1)*x(2);1.08*x(2)-0.08*x(1)*x(2)],[2 0 0 0]);
xRV = nefGaussianRV(xmean,Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmean = nefLinFunction(H,[],[]);
zRV = nefGaussianRV(zmean,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model (for PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system_pdf = nefPDFSystem(xRV,zRV,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: UKF, KALMAN, PF ...')
UKF = nefUKF(system);
UKF_no_tests = nefUKF(system,'disableChecks',1);

KALMAN = nefKalman(system);
PF = nefPF(system_pdf,'sampleSize',1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up smoothers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up smoothers: KALMAN, UKF ...')
RTSS = nefKalman(system,'taskType','fixedLagSmoothing','taskPar',4);
URTSS = nefUKF(system,'scalingParameterType','constant','parameterValue',2,'taskType','fixedLagSmoothing','taskPar',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running UKF filter ...')
t = cputime;
for i = 1:100
[val_UKF] = estimate(UKF,z,[]);
end
UKFtime = cputime-t
t = cputime;
for i = 1:100
[val_UKF] = estimate(UKF_no_tests,z,[]);
end
UKFtime_nochecks = cputime-t
disp('Running KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,[]);
KALMANtime = cputime-t;
disp('Running PF filter ...')
t = cputime;
[val_PF] = estimate(PF,z,[]);
PFtime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running KALMAN smoother ...')
t = cputime;
[val_RTSS] = estimate(RTSS,z,[]);
RTSStime = cputime-t;
disp('Running UKF smoother ...')
t = cputime;
[val_URTSS] = estimate(URTSS,z,[]);
URTSStime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = 1:K
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_PF(:,i) = evalMean(val_PF{i});

  msem_UKF(:,i) = (xest_UKF(:,i)-x(:,i)).^2; 
  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_PF(:,i) = (xest_PF(:,i)-x(:,i)).^2; 
end
t4 = [1:K-4];
for i = 1:K-4
  xest_RTSS(:,i) = evalMean(val_RTSS{i});
  xest_URTSS(:,i) = evalMean(val_URTSS{i});
  
  msem_RTSS(:,i) = (xest_RTSS(:,i)-x(:,i)).^2; 
  msem_URTSS(:,i) = (xest_URTSS(:,i)-x(:,i)).^2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Stats : MSEM\t\t time\n')
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('PF    : %f\t%f\n',mean(mean(msem_PF)),PFtime);
fprintf('RTSS  : %f\t%f\n',mean(mean(msem_RTSS)),RTSStime);
fprintf('URTSS : %f\t%f\n',mean(mean(msem_URTSS)),URTSStime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,1)
plot(t,x(1,:),'b',t,xest_UKF(1,:),'r--',t,xest_KALMAN(1,:),'k--',t4,xest_RTSS(1,:),'k*',t4,xest_URTSS(1,:),'ro')
ylabel('Number of preys')
xlabel('time')
legend('true','UKF','KALMAN','RTSS','URTSS')
subplot(2,1,2)
plot(t,x(2,:),'b',t,xest_UKF(2,:),'r--',t,xest_KALMAN(2,:),'k--',t4,xest_RTSS(2,:),'k*',t4,xest_URTSS(2,:),'ro')
ylabel('Number of predators')
xlabel('time')
legend('true','UKF','KALMAN','RTSS','URTSS')
