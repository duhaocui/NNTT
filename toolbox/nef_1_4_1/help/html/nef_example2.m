%% Estimation experiment for probabilistically described system
%
% In this example a ship is moving in x-y plane with the stationary
% observer at the origin of the plane, The ship is assumed to accelerate
% and decelerate randomly over time. The state $x_k$ is given as ${{\bf
% x}}_k=\left(x_k,\dot{x}_k,y_k,\dot{y}_k\right)^T$, where $x_k,y_k$
% represent position of the ship and $\dot{x}_k,\dot{y}_k$ are
% corresponding velocities. The state dynamics is described by the
% following transition pdf
%
% <<./graphics/ex2-transition-pdf.png>>
%
% where $I_4$ is the 4-by-4 identity matrix.  The measurement of the model
% is the direction given by $\tan^{-1}(\frac{y_k}{x_k})$ and the
% corresponding measurement pdf is given as
%
% <<./graphics/ex2-measurement-pdf.png>>
%
% and the pdf of the initial state ${\bf x}_{0}$ is
%
% <<./graphics/ex2-x0-pdf.png>>
%
% The classes necessary for description of this problem is depicted below.
%
% <<./graphics/system-example2.png>>
%
%%
% For specification of the system, the |nefPDFSystem| class will be
% used. First, the transition pdf is set using its mean and covariance matrix as
F = [1 1 0 0;0 1 0 0 ;0 0 1 1 ;0 0 0 1];
xMean = nefLinFunction(F,[],[]);
xVariance = 0.0001*eye(4);
xPdf = nefGaussianRV(xMean,xVariance);
%%
% Then, the measurement pdf is also given by its mean and variance as
mFun = @(x,u,v,t) atan(x(3)/x(1));
zMean = nefHandleFunction(mFun,[4 0 0 0]);
zVariance = 0.0001;
zPdf = nefGaussianRV(zMean,zVariance);
%%
% Now, the initial condition is specified as 
x0Pdf = nefGaussianRV([-0.05 0.001 2 -0.055]',0.01*eye(4));
%%
% Consequently, the system is created using probabilistic description as
model = nefPDFSystem(xPdf,zPdf,x0Pdf);
%%
% and simulated using commands
nSteps=20;
[z,x] = simulate(model,nSteps,[]);
%%
% where the measurement values are stored in |z| and values of the
% state $x_k$ in |x|. Finally, to estimate the state $x_k$
% using the measurement $z^k$, the auxiliary particle filter will be
% used. It will be built by a single command specifying system description
% (|system|), lag |0| and type of the sampling density:
pfEstimator = nefPF(model,'samplingDensity','pointAuxiliary');
%%
% Note that no other design parameters concerning the particle filter were
% specified and thus the default values were used, such as default sample
% size is $100$ samples, systematic resampling executed at each time
% instant.
%% 
% The actual estimation is issued by
[estimates] = estimate(pfEstimator,z,[]);
%%
% As a result, estimates of the filtering pdf $p(x_k|z_k)$ given by
% the empirical pdf |nefEmpiricalRV| are stored in the |estimates|.

