
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear system with unknown parameter (x2) in state equation
% x(k+1) = [x1(k)^2*x2(k);x2(k)] + w(k)
% z(k) = x1(k) + v(k)
% p(w(k)) = N{w(k):0,Q}
% p(v(k)) = N{w(k):0,R}
% p(x(0)) = N{x(0):0,Q}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% UKF - Unscented Kalman filter
% KALMAN - Kalman filter
% DD1 - Divide difference filter 1st-order
% PF - Particle filter

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all

PWD = pwd;
addpath([PWD,'/support/']);

fprintf('##################################################\n')
fprintf('# EXAMPLE 1 (Nonlinear system with unknown parameter in state equation)\n')
fprintf('##################################################\n')
% #time steps
K = 20;

disp('Creating system and calculating its trajectory ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = nefHandleFunction(@(x,u,w,k) [x(1)^2*x(2)+w(1);x(2)+w(2)],[2 0 2 0],'diff1Noise',@(x,u,v,k) eye(2),'diff1State',@(x,u,w,k) [2*x(1)*x(2) x(1)^2;0 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [1 0];
h = nefLinFunction(H,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = eye(2)*0.05;
w = nefGaussianRV([0 0]',Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 0.01;
v = nefGaussianRV(0,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([0.9;-0.85],diag([1e-2 1e-5]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system and simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);
[z,x] = simulate(system,K,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% probabilistic description for PF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmean =  nefHandleFunction(@(x,u,w,k) [x(1)^2*x(2);x(2)],[2 0 0 0]);
xRV = nefGaussianRV(xmean,Q);
zmean = nefLinFunction(H,[],[]);
zRV = nefGaussianRV(zmean,R);
system_pdf = nefPDFSystem(xRV,zRV,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: UKF, KALMAN, DD1, PF ...')
UKF = nefUKF(system);
KALMAN = nefKalman(system);
DD1 = nefDD1(system);
PF = nefPF(system_pdf,'sampleSize',100,'resamplingSched','dynamic','samplingDensity','optimal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running UKF filter ...')
t = cputime;
[val_UKF] = estimate(UKF,z,[]);
UKFtime = cputime-t;
disp('Running KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,[]);
KALMANtime = cputime-t;
disp('Running DD1 filter ...')
t = cputime;
[val_DD1] = estimate(DD1,z,[]);
DD1time = cputime-t;
disp('Running PF filter ...')
t = cputime;
[val_PF] = estimate(PF,z,[]);
PFtime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = 1:K
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_DD1(:,i) = evalMean(val_DD1{i});
  xest_PF(:,i) = evalMean(val_PF{i});

  msem_UKF(:,i) = (xest_UKF(:,i)-x(:,i)).^2; 
  msem_DD1(:,i) = (xest_DD1(:,i)-x(:,i)).^2; 
  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_PF(:,i) = (xest_PF(:,i)-x(:,i)).^2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Stats : MSEM\t\t time\n')
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('DD1   : %f\t%f\n',mean(mean(msem_DD1)),DD1time);
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('PF    : %f\t%f\n',mean(mean(msem_PF)),PFtime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)
plot(t,x(1,:),'b',t,xest_UKF(1,:),'r--',t,xest_KALMAN(1,:),'k--',t,xest_DD1(1,:),'g--',t,xest_PF(1,:),'g:')
ylabel('state')
xlabel('time')
legend('true','UKF','KALMAN','DD1','PF')
subplot(2,1,2)
plot(t,x(2,:),'b',t,xest_UKF(2,:),'r--',t,xest_KALMAN(2,:),'k--',t,xest_DD1(2,:),'g--',t,xest_PF(2,:),'g:')
ylabel('parameter')
xlabel('time')
legend('true','UKF','KALMAN','DD1','PF')


figure
subplot(2,1,1)
pdfTimePlot(val_UKF,'trajectory',x(1,:),'index',1,'fill',1);
subplot(2,1,2)
pdfTimePlot(val_UKF,'trajectory',x(2,:),'index',2,'fill',1);

rmpath([PWD,'/support/']);
