%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural networks MLP - 20 hidden neurons in one layer
% 4 input neurons
% data obtained from 
% [1] I-Cheng Yeh, "Modeling of strength of high performance concrete using artificial neural networks,
% " Cement and Concrete Research, Vol. 28, No. 12, pp. 1797-1808 (1998). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% EKF - extended Kalman filter

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 8 (Neural network training)\n')
fprintf('##################################################\n')

addpath([pwd,'/support/']);
load('data/concrete_normed')
Ndata = size(data,1);
Ntest = 200;
Ktrain = Ndata - Ntest;
idxTest = Ktrain+1:Ndata;
disp('Creating system and calculating its trajectory ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nh = 20; % #hidden neurons
ni = 9; % #input neurons
np = ni*nh+nh+1; % #unknown parameters
F = eye(np);
f = nefLinFunction(F,[],eye(np));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = nefHandleFunction(@(x,u,v,t) x(end) + x(np-nh:np)'*[tanh(reshape(x(1:nh*ni),ni,nh)'*u);1] + v,[np ni 1 0],...
    'diff1Noise',@(x,u,v,k) 1,'diff1State',@(x,u,w,k) derivace_nn(x,u,w,k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = eye(np)*1e-5;
w = nefGaussianRV(zeros(np,1),Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 1e-2;
v = nefGaussianRV(0,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m0 = rand(np,1)*2-1;
P0 = zeros(np,np);
x0 = nefGaussianRV(m0,P0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = nefEqSystem(f,h,w,v,x0);
u = [data(1:Ktrain,1:8)';ones(1,Ktrain)];
z = data(1:Ktrain,9)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: EKF ...')
EKF = nefKalman(system);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running EKF filter ...')
t = cputime;
[val_EKF] = estimate(EKF,z,u);
EKFtime = cputime-t;

for i = 1:Ktrain
  xest_EKF(:,i) = evalMean(val_EKF{i});
end

utest = [data(idxTest,1:8)';ones(1,Ntest)];
ztest = data(idxTest,9)';
for i = 1:Ntest
  zNN(i) =  evaluate(h,xest_EKF(:,end),utest(:,i),0,[]);
end

t = [1:Ntest];
plot(t,ztest,'b',t,zNN,'r',t,(ztest-zNN).^2,'k')
xlabel('time')
ylabel('output')
legend('true output','NN model output','error')

rmpath([pwd,'/support/']);
