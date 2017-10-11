%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% track-before-detect application with target appearing at each time instant
% x(k) = [xpos xvel ypos yvel I]
% xpos - position in x-direction
% xvel - velocity in x-direction
% ypos - position in y-direction
% yvel - velocity in y-direction
% I - intensity of the target
% x(k+1) = F*x(k) + w(k)
% measuring intensity of signal within each of nx*ny cells with dimensions dx and dy
% target contribution to intensity is modelled by a two-dimensional Gaussian distribution
% for NEF purposes, the measurement is stacked up as a column
% h(ix,iy) = dx*dy*I/2/pi/Sigma^2 * exp(-((ix*dx-xpos)^2(iy*dy-ypos)^2)/2/Sigma^2)
% z(k) = h(:) + v;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimators:
% PF - Particle filter with prior sampling density

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

clear all
fprintf('##################################################\n')
fprintf('# EXAMPLE 6 (track-before-detect)\n')
fprintf('##################################################\n')

% #time steps
K = 20;

disp('Creating system and calculating its trajectory ...')
% STRUCTURAL DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 1;% revisit time
F = [1 T 0 0 0;
     0 1 0 0 0;
     0 0 1 T 0;
     0 0 0 1 0;
     0 0 0 0 1];

f = nefLinFunction(F,[],eye(5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 1;
dy = 1;
nx = 20;
ny = 20;
Sigma = 0.7;
xaxis = [1:nx]'*dx;
yaxis = [1:ny]*dy;
h_handle = @(x,u,v,t) reshape(dx*dy*x(5)/2/pi/Sigma^2 * exp(-(repmat((xaxis-x(1)).^2,1,ny)+repmat((yaxis-x(3)).^2,nx,1))/2/Sigma^2),nx*ny,1)+v;
h = nefHandleFunction(h_handle,[5 0 nx*ny 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1 = 0.001;
q2 = 0.01;
Q = [q1/3*T^3 q1/2*T^2 0        0        0;
     q1/2*T^2 q1*T     0        0        0;
     0        0        q1/3*T^3 q1/2*T^2 0;
     0        0        q1/2*T^2 q1*T     0;
     0        0        0        0        q2*T*100];
w = nefGaussianRV(zeros(5,1),Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 3;
v = nefGaussianRV(zeros(nx*ny,1),eye(nx*ny)*sigma^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([4.2;0.45;7.2;0.25;20],zeros(5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system=nefEqSystem(f,h,w,v,x0);

% PROBABILISTIC DESCRIPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmean = nefLinFunction(F,[],[]);
xvar = Q;
xRV = nefGaussianRV(xmean,xvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurement pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmean_handle = @(x,u,v,t) reshape(dx*dy*x(5)/2/pi/Sigma^2 * exp(-(repmat((xaxis-x(1)).^2,1,ny)+repmat((yaxis-x(3)).^2,nx,1))/2/Sigma^2),nx*ny,1);
zmean = nefHandleFunction(zmean_handle,[5 0 0 0]);
zvar = eye(nx*ny)*sigma^2;
zRV = nefGaussianRV(zmean,zvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model (for PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_pdf = nefPDFSystem(xRV,zRV,x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating its trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = [];
[z,x] = simulate(system,K,u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting up filters:PF ...')
PF = nefPF(model_pdf,'sampleSize',100,'resamplingSched','dynamic','samplingDensity','prior');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  xest_PF(:,i) = evalMean(val_PF{i});

  msem_PF(:,i) = (xest_PF(:,i)-x(:,i)).^2; 
end

fprintf('Stats : MSEM\t\t time\n')
fprintf('PF    : %f\t%f\n',mean(mean(msem_PF)),PFtime);


subplot(3,1,1)
plot(t,x(1,:),'b',t,xest_PF(1,:),'r*')
ylabel('x - position')
xlabel('time')
legend('true','PF')

subplot(3,1,2)
plot(t,x(3,:),'b',t,xest_PF(3,:),'r*')
ylabel('y - position')
xlabel('time')
legend('true','PF')

subplot(3,1,3)
plot(t,x(5,:),'b',t,xest_PF(5,:),'r*')
ylabel('intensity')
xlabel('time')
legend('true','PF')
