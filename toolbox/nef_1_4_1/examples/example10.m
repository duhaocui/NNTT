%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fifth-order two-phase nonlinear model of induction motor
% state: x(k) = [i_sa i_sb fi_ra fi_rb om]
% i_sa - stator current a
% i_sb - stator current b
% fi_ra - rotor flux a
% fi_rb - rotor flux b
% om - angular speed
% measurement z = [i_sa i_sb]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% UKF - Unscented Kalman filter
% DD1 - Divide difference filter 1st-order
% DD2 - Divide difference filter 2nd-order

% NEF version 1.2.0
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 10 (fifth-order two-phase nonlinear model of induction motor)\n')
fprintf('##################################################\n')

% #time steps
N = 1000;

disp('Creating system and calculating its trajectory ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rs = 0.18; % [Ohm]; stator per-phase resistance
Rr = 0.15; % [Ohm]; rotor per-phase resistance
M = 0.068; % [H]
Ls = 0.0699; % [H]; stator per-phase inductance
Lr = 0.0699; % [H]; rotor per-phase inductance
J = 0.0586; % [kg*m^2]; rotor moment of inertia
Tl = 10; % [Nm]
p = 1;
h = 1e-4; % [s]
Tr = Lr/Rr; % rotor time constant
sigma = 1-M^2/(Ls*Lr); 
K = M/(sigma*Ls*Lr);
gamma = Rs/(sigma*Ls)+Rr*M^2/(sigma*Ls*Lr^2);

f = nefHandleFunction(@(x,u,w,t) [x(1)+h*(-gamma*x(1)+K/Tr*x(3)+K*p*x(5)*x(4)+1/(sigma*Ls)*u(1)); ...
                                  x(2)+h*(-gamma*x(2)+K/Tr*x(4)-K*p*x(5)*x(3)+1/(sigma*Ls)*u(2)); ...
                                  x(3)+h*(M/Tr*x(1)-1/Tr*x(3)-p*x(5)*x(4)); ...
                                  x(4)+h*(M/Tr*x(2)-1/Tr*x(4)+p*x(5)*x(3)); ...
                                  x(5)+h*(p*M/(J*Lr)*(x(3)*x(2)-x(4)*x(1))-Tl/J)] + w,[5 2 5 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = nefHandleFunction(@(x,u,v,t) [x(1); x(2)] + v,[5 0 2 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 1e-1*eye(5);
w = nefGaussianRV(zeros(5,1),Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurment noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 5e-1*eye(2);
v = nefGaussianRV(zeros(2,1),R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%m0 = [200 200 50 50 300]';
m0 = [0 0 0 0 0]';
P0 = 1e-1*eye(5);
x0 = nefGaussianRV(m0,P0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = nefEqSystem(f,h,w,v,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = [350*cos(0.3*(0:N-1)); 300*sin(0.3*(0:N-1))];
[z,x] = simulate(system,N,u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters: UKF, DD1, DD2 ...')
UKF = nefUKF(system,'scalingParameterType','constant','parameterValue',1);
DD1 = nefDD1(system);
DD2 = nefDD2(system,'scalingParameterType','constant','parameterValue',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
t = [1:N];
for i = 1:N
  xest_UKF(:,i) = evalMean(val_UKF{i});
  xest_DD1(:,i) = evalMean(val_DD1{i});  
  xest_DD2(:,i) = evalMean(val_DD2{i});  
end

for i = 1:N
  msem_UKF(:,i) = (xest_UKF(:,i)-x(:,i)).^2; 
  msem_DD1(:,i) = (xest_DD1(:,i)-x(:,i)).^2; 
  msem_DD2(:,i) = (xest_DD2(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('UKF   : %f\t%f\n',mean(mean(msem_UKF)),UKFtime);
fprintf('DD1   : %f\t%f\n',mean(mean(msem_DD1)),DD1time);
fprintf('DD2   : %f\t%f\n',mean(mean(msem_DD2)),DD2time);

subplot(3,2,1)
plot(t,x(1,:),'b',t,xest_DD1(1,:),'r*',t,xest_UKF(1,:),'g+',t,xest_DD2(1,:),'mo')
xlabel('time')
ylabel('i_sa')
legend('true','DD1','UKF','DD2')
subplot(3,2,2)
plot(t,x(2,:),'b',t,xest_DD1(2,:),'r*',t,xest_UKF(2,:),'g+',t,xest_DD2(2,:),'mo')
xlabel('time')
ylabel('i_sb')
legend('true','DD1','UKF','DD2')
subplot(3,2,3)
plot(t,x(3,:),'b',t,xest_DD1(3,:),'r*',t,xest_UKF(3,:),'g+',t,xest_DD2(3,:),'mo')
xlabel('time')
ylabel('fi_ra')
legend('true','DD1','UKF','DD2')
subplot(3,2,4)
plot(t,x(4,:),'b',t,xest_DD1(4,:),'r*',t,xest_UKF(4,:),'g+',t,xest_DD2(4,:),'mo')
xlabel('time')
ylabel('fi_rb')
legend('true','DD1','UKF','DD2')
subplot(3,2,5)
plot(t,x(5,:),'b',t,xest_DD1(5,:),'r*',t,xest_UKF(5,:),'g+',t,xest_DD2(5,:),'mo')
xlabel('time')
ylabel('om')
legend('true','DD1','UKF','DD2')
