
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear system with unknown parameter (x2) in state equation
% x(k+1) = x(k)-0.2*x(k)^2 + w(k)
% z(k) = x(k)^2 + v(k)
% p(w(k)) = N{w(k):0,Q}
% p(v(k)) = N{w(k):0,R}
% p(x(0)) = N{x(0):0,Q}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% KALMAN - Kalman filter
% GSM - Gaussian sum filter

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all

fprintf('##################################################\n')
fprintf('# EXAMPLE 11 (Scalar nonlinear system)\n')
fprintf('##################################################\n')

addpath([pwd,'/support/']);

% #time steps
K = 5;

disp('Creating system and calculating its trajectory ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = nefHandleFunction(@(x,u,w,k) x-0.2*x^2+w,[1 0 1 0],'diff1Noise',@(x,u,v,k) 1,'diff1State',@(x,u,w,k) 1-0.4*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = nefHandleFunction(@(x,u,v,k) x^2+v,[1 0 1 0],'diff1Noise',@(x,u,v,k) 1,'diff1State',@(x,u,w,k) 2*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = nefGaussianRV(0,0.25);
wm = nefGaussianSumRV(1,0,0.25,'parameters','triplets');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = nefGaussianRV(0,1);
vm = nefGaussianSumRV(1,0,1,'parameters','triplets');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV(0,0.5);
x0m = nefGaussianSumRV(0.5,-1,0.25,0.5,1,0.25,'parameters','triplets');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system and simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = nefEqSystem(f,h,wm,vm,x0m);
[z,x] = simulate(system,K,[]);
model = nefEqSystem(f,h,w,v,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: KALMAN, GSM ...')
KALMAN = nefKalman(model);
GSM = nefGSM(system);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,[]);
KALMANtime = cputime-t;
disp('Running GSM filter ...')
t = cputime;
[val_GSM] = estimate(GSM,z,[]);
GSMtime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = 1:K
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_GSM(:,i) = evalMean(val_GSM{i});

  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_GSM(:,i) = (xest_GSM(:,i)-x(:,i)).^2; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Stats : MSEM\t\t time\n')
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('GSM   : %f\t%f\n',mean(mean(msem_GSM)),GSMtime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
pdfTimePlot(val_GSM,'trajectory',x(1,:),'index',1,'fill',1);
subplot(2,1,2)
pdfTimePlot(val_KALMAN,'trajectory',x(1,:),'index',1,'fill',1);


rmpath([pwd,'/support/']);