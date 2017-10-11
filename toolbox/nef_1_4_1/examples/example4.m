%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear Gaussian system with unknown parameter, i.e. nonlinear Gaussian system
% x(k+1) = a*x(k) + w(k)
% z(k) = H*x(k) + v(k)
% p(w(k)) = N(w(k):0,Q)
% p(v(k)) = N(v(k):0,R)
% p(x(0)) = N(x(0):m0,P0)
%
% the unknown parameter is used as a new state component
% i.e. x(2) = a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% UKF - Unscented Kalman filter
% SUKF - square-root version of Unscented Kalman filter
% KALMAN - Kalman filter
% DD1 - Divide difference filter 1st-order
% sDD1 - square-root version of Divide difference filter 1st-order
% ITEK - iteration Kalman filter

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 4 (Linear system with unknown parameter in state equation)\n')
fprintf('##################################################\n')

% #time steps
K = 100;

disp('Creating system and calculating its trajectory ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1state = @(x,u,w,t) [x(2) x(1);0 1];
d1noise = @(x,u,w,t) [1 0;0 1];
f = nefHandleFunction(@(x,u,w,t) [x(1)*x(2)+w(1);x(2)+w(2)],[2,0,2,0],'Diff1State',d1state,'Diff1Noise',d1noise);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [1 0];
h = nefLinFunction(H,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = [0.1 0;0 0];
w = nefGaussianRV([0;0],Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 0.01;
v = nefGaussianRV([0],R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = [1;0.8];
P0 = [0.1 0;0 1e-10];
x0 = nefGaussianRV(m0,P0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system and simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);
u = [];
[z,x] = simulate(system,K,u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: UKF, SUKF, DD1, SDD1, KALMAN, ITEK ...')
UKF = nefUKF(system,'sqrt','svd');   % SVD decomposition necessary because Q is possitive semidefinite
SUKF = nefSUKF(system);
DD1 = nefDD1(system);
SDD1 = nefSDD1(system);
KALMAN = nefKalman(system);
ITEK = nefItFilter(system,'localFilter','kalman','epsilon',0.1,'maxIter',200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running UKF filter ...')
t = cputime;
[val_UKF] = estimate(UKF,z,u);
UKFtime = cputime-t;
disp('Running SUKF filter ...')
t = cputime;
[val_SUKF] = estimate(SUKF,z,u);
SUKFtime = cputime-t;
disp('Running DD1 filter ...')
t = cputime;
[val_DD1] = estimate(DD1,z,u);
DD1time = cputime-t;
disp('Running SDD1 filter ...')
t = cputime;
[val_SDD1] = estimate(SDD1,z,u);
SDD1time = cputime-t;
disp('Running KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,u);
KALMANtime = cputime-t;
disp('Running ITEK filter ...')
t = cputime;
[val_ITEK] = estimate(ITEK,z,u);
ITEKtime = cputime-t;
disp('Summarizing results ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = t
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_SUKF(:,i) = evalMean(val_SUKF{i});
  xest_DD1(:,i) = evalMean(val_DD1{i});
  xest_SDD1(:,i) = evalMean(val_SDD1{i});
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_ITEK(:,i) = evalMean(val_ITEK{i});

  msem_UKF(:,i) = (xest_UKF(:,i)-x(:,i)).^2; 
  msem_SUKF(:,i) = (xest_SUKF(:,i)-x(:,i)).^2; 
  msem_DD1(:,i) = (xest_DD1(:,i)-x(:,i)).^2; 
  msem_SDD1(:,i) = (xest_SDD1(:,i)-x(:,i)).^2; 
  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_ITEK(:,i) = (xest_ITEK(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('SUKF  : %f\t%f\n',mean(mean(msem_SUKF)),SUKFtime);
fprintf('DD1   : %f\t%f\n',mean(mean(msem_DD1)),DD1time);
fprintf('SDD1  : %f\t%f\n',mean(mean(msem_SDD1)),SDD1time);
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('ITEK  : %f\t%f\n',mean(mean(msem_ITEK)),ITEKtime);

