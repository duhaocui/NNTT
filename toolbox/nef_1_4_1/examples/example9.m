%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple models with linear state and measurement equations
% x(k+1) = F_i x(k) + w(k)
% z(k) = H x(k) + v(k)
% with corresponding probabilities alpha_i
%
% described by GS models
% p(x(k+1)|x(k)) = sum_{i=1}^{N_i} alpha_i N{x(k+1):F_i x(k),Q}
% p(z(k)|x(k)) = N{z(k):H x(k),R}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% GSM - Gaussian sum method

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 9 (Linear multiple models)\n')
fprintf('##################################################\n')
K = 10;

disp('Creating system and calculating its trajectory ...')
% PROBABILISTIC DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 2; % #models
% 1st model
F1 = eye(nx);
xm1 = nefLinFunction(F1,[],[]); % mean
Q1 = eye(nx)*1e-2; % variance
xmodel1 = nefGaussianRV(xm1,Q1); % model
% 2nd model
F2 = eye(nx) * 0.9;
xm2 = nefLinFunction(F2,[],[]); % mean
Q2 = eye(nx)*1e-2; % variance
xmodel2 = nefGaussianRV(xm2,Q2); % model
% weights of the model
w1 = 0.2;
w2 = 0.8;
xRV = nefGaussianSumRV(w1,xmodel1,w2,xmodel2,'parameters','wnefgaussianrv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz = 2;
H = eye(nz);
R = 1e-2*eye(nz);
h = nefLinFunction(H,[],[]);
zmodel = nefGaussianRV(h,R);
zRV = nefGaussianSumRV(1,zmodel,'parameters','wnefgaussianrv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = [0.6;-1];
P0 = 1*eye(nx);
x0RV = nefGaussianSumRV(1,m0,P0,'parameters','triplets');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = nefPDFSystem(xRV,zRV,x0RV);
[z,x] = simulate(system,K,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: GSM ...')
GSM = nefGSM(system,'localFilter','ukf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  xest_GSM(:,i) = evalMean(val_GSM{i});
  msem_GSM(:,i) = (xest_GSM(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('GSM   : %f\t%f\n',mean(mean(msem_GSM)),GSMtime);
