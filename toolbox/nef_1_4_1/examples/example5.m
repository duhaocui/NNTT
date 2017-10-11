%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear Gaussian system with unknown parameter in state noise, i.e. nonlinear Gaussian system
% x(k+1) = F*x(k) + a*w(k)
% z(k) = H*x(k) + v(k)
% p(w(k)) = N(w(k):0,Q)
% p(v(k)) = N(v(k):0,R)
% p(x(0)) = N(x(0):m0,P0)
%
% the unknown parameter is used as a new state component
% i.e. x(2) = a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% KALMAN - Kalman filter
% PF - Particle filter with prior sampling density

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 5 (Linear system with unknown parameter in state noise)\n')
fprintf('##################################################\n')

% #time steps
K = 10;

disp('Creating system and calculating its trajectory ...')
% STRUCTURAL DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = 0.9;
d1state = @(x,u,w,t) [F 0;0 1];
d1noise = @(x,u,w,t) [x(2) 0;0 1];
f = nefHandleFunction(@(x,u,w,t) [F*x(1)+x(2)*w(1);x(2)+w(2)],[2,0,2,0],'Diff1State',d1state,'Diff1Noise',d1noise);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [1 0];
h = nefLinFunction(H,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = [0.1 0;0 1e-10];
w = nefGaussianRV([0;0],Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 0.01;
v = nefGaussianRV([0],R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([1;9],1e1*eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);

% PROBABILISTIC DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmean = nefHandleFunction(@(x,u,w,t) [F*x(1);x(2)],[2 0 0 0]);
xvar = nefHandleFunction(@(x,u,w,t) [x(2)^2*Q(1,1) Q(1,2);Q(2,1) Q(2,2)],[2 0 0 0]);
xRV = nefGaussianRV(xmean,xvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmean = nefLinFunction(H,[],[]);
zvar = R;
zRV = nefGaussianRV(zmean,zvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model (for PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_pdf = nefPDFSystem(xRV,zRV,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = [];
[z,x] = simulate(system,K,u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters:KALMAN, PF ...')
KALMAN = nefSDD1(system);
PF = nefPF(model_pdf,'sampleSize',1000,'resamplingSched','dynamic','samplingDensity','optimal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running (extended) KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,u);
KALMANtime = cputime-t;
disp('Running PF filter ...')
t = cputime;
[val_PF] = estimate(PF,z,u);
PFtime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = t
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_PF(:,i) = evalMean(val_PF{i});

  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_PF(:,i) = (xest_PF(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('PF    : %f\t%f\n',mean(mean(msem_PF)),PFtime);

subplot(2,1,1)
plot(t,x(1,:),'b',t,xest_KALMAN(1,:),'r*',t,xest_PF(1,:),'k+')
ylabel('state')
xlabel('time')
legend('true','KALMAN','PF')

subplot(2,1,2)
plot(t,x(2,:),'b',t,xest_KALMAN(2,:),'r*',t,xest_PF(2,:),'k+')
ylabel('parameter')
xlabel('time')
legend('true','KALMAN','PF')
