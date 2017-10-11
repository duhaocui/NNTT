%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear Gaussian system
% x(k+1) = F*x(k) + G*u(k) + w(k)
% z(k) = H*x(k) + v(k)
% nx = 2, nz = 2
% p(w(k)) = N{w(k):0,Q}
% p(v(k)) = N{w(k):0,R}
% p(x(0)) = N{x(0):m0,P0}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% UD - UD version of Kalman filter
% UKF - Unscented Kalman filter
% KALMAN - Kalman filter
% DD1 - Divide difference filter 1st-order
% DD2 - Divide difference filter 2nd-order
% GSM - Gaussian Sum Method
% PF - Particle filter with prior sampling density
% PF_EKF - Particle filter with EKF sampling density

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 3 (linear Gaussian system)\n')
fprintf('##################################################\n')

% #time steps
K = 20;

disp('Creating system and calculating its trajectory ...')
% STRUCTURAL DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [0.9 0.01;0.01 0.9];
G = eye(2);
f = nefLinFunction(F,G,eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [5 0;0 1];
h = nefLinFunction(H,eye(2),eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 0.05*[1 2;2 5]; %eye(2);
w = nefGaussianRV([0;0],Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 0.01*[3 1;1 4];
v = nefGaussianRV([0;0],R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([1;10],1e-1*eye(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system and simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);
u = [ones(1,10) -1*ones(1,10);-1*ones(1,10) ones(1,10)];
% example of using pregenerated values of state noise in system simulation
wValues = drawSample(w,K);
[z,x] = simulate(system,K,u,'stateNoise',wValues);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative model - for PF_EKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ltpdf = @(nx,x,u,t) - log(2*pi) - 0.5*log(det(Q)) - 0.5*sum(((nx-F*x-G*u)'*inv(Q))'.*(nx-F*x-G*u));
llpdf = @(z,x,u,t)  - log(2*pi) - 0.5*log(det(R)) - 0.5*sum(((z-H*x)'*inv(R))'.*(z-H*x));
model2=nefEqSystem(f,h,w,v,x0,'logLikelihood',llpdf,'logTransitionPDF',ltpdf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model for GSM (all RV's must me in the form of Gaussian sums)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mx0 = nefGaussianSumRV(1,x0,'parameters','wnefgaussianrv');
mw = nefGaussianSumRV(1,w,'parameters','wnefgaussianrv');
mv = nefGaussianSumRV(1,v,'parameters','wnefgaussianrv');
model_gs=nefEqSystem(f,h,mw,mv,mx0);
% PROBABILISTIC DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmean =  nefLinFunction(F,G,[]);
xRV = nefGaussianRV(xmean,Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmean = nefLinFunction(H,[],[]);
zRV = nefGaussianRV(zmean,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model (for PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_pdf = nefPDFSystem(xRV,zRV,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: UD, UKF, DD1, DD2, KALMAN, GSM, PF_EKF, PF ...')
UD = nefDD1(system);
UKF = nefUKF(system);
DD1 = nefDD1(system);
DD2 = nefDD2(system);
KALMAN = nefKalman(system);
GSM = nefGSM(model_gs);
PF_EKF = nefPF(model2,'sampleSize',100,'resamplingSched','dynamic','samplingDensity','ekf');
PF = nefPF(model_pdf,'sampleSize',1000,'resamplingSched','dynamic','samplingDensity','prior');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running UD filter ...')
t = cputime;
[val_UD] = estimate(UD,z,u);
UDtime = cputime-t;
disp('Running UKF filter ...')
t = cputime;
[val_UKF] = estimate(UKF,z,u);
UKFtime = cputime-t;
t = cputime;
disp('Running DD1 filter ...')
[val_DD1] = estimate(DD1,z,u);
DD1time = cputime-t;
t = cputime;
disp('Running DD2 filter ...')
[val_DD2] = estimate(DD2,z,u);
DD2time = cputime-t;
disp('Running KALMAN filter ...')
t = cputime;
[val_KALMAN] = estimate(KALMAN,z,u);
KALMANtime = cputime-t;
disp('Running PF_EKF filter ...')
t = cputime;
[val_PF_EKF] = estimate(PF_EKF,z,u);
PF_EKFtime = cputime-t;
disp('Running PF filter ...')
t = cputime;
[val_PF] = estimate(PF,z,u);
PFtime = cputime-t;
disp('Running GSM filter ...')
t = cputime;
[val_GSM] = estimate(GSM,z,u);
GSMtime = cputime-t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means and MSEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing statistics of the obtained results ...')
t = [1:K];
for i = 1:K
  xest_UD(:,i) = evalMean(val_UD{i});
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_DD1(:,i) = evalMean(val_DD1{i});
  xest_DD2(:,i) = evalMean(val_DD2{i});
  xest_KALMAN(:,i) = evalMean(val_KALMAN{i});
  xest_GSM(:,i) = evalMean(val_GSM{i});
  xest_PF_EKF(:,i) = evalMean(val_PF_EKF{i});
  xest_PF(:,i) = evalMean(val_PF{i});

  msem_UD(:,i) = (xest_UD(:,i)-x(:,i)).^2; 
  msem_UKF(:,i) = (xest_UKF(:,i)-x(:,i)).^2; 
  msem_DD1(:,i) = (xest_DD1(:,i)-x(:,i)).^2; 
  msem_DD2(:,i) = (xest_DD2(:,i)-x(:,i)).^2; 
  msem_KALMAN(:,i) = (xest_KALMAN(:,i)-x(:,i)).^2; 
  msem_GSM(:,i) = (xest_GSM(:,i)-x(:,i)).^2; 
  msem_PF_EKF(:,i) = (xest_PF_EKF(:,i)-x(:,i)).^2; 
  msem_PF(:,i) = (xest_PF(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('UD    : %f\t%f\n',mean(mean(msem_UD)),UDtime);
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('DD1   : %f\t%f\n',mean(mean(msem_DD1)),DD1time);
fprintf('DD2   : %f\t%f\n',mean(mean(msem_DD2)),DD2time);
fprintf('KALMAN: %f\t%f\n',mean(mean(msem_KALMAN)),KALMANtime);
fprintf('PF_EKF: %f\t%f\n',mean(mean(msem_PF_EKF)),PF_EKFtime);
fprintf('PF    : %f\t%f\n',mean(mean(msem_PF)),PFtime);
fprintf('GSM   : %f\t%f\n',mean(mean(msem_GSM)),GSMtime);

