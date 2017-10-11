%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate of Noise Covariance Matrices - ACD_2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear system with unknown noise covariance matricis (Q, R)
% x(k+1) = F*x(k) + N*u(k) + M*w(k)
% z(k) = H*x(k) + O*v(k)
% p(w(k)) = N{w(k):0,Q}
% p(v(k)) = N{v(k):0,R}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used method:
% ALS - Autocovariance least-square method
% KF - Kalman filter
%
% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia
clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
F = 0.9;
N = 0.5;
M = 2;
f = nefLinFunction(F,N,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2;
H = [1;0.5];
O = [2 0.1;0.1 1];
h = nefLinFunction(H,[],O);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV(zeros(n,1),eye(n));

L=3;
% varargin: 'pointsCrit','numberCovEq','beginPointCrit','omittedInitPeriod'
myALS=nefALS(f,h,x0,'numberCovEq',L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 1.5;
w = nefGaussianRV(zeros(size(Q,1),1),Q);
%w = nefUniformRV(-sqrt(Q*12)/2,sqrt(Q*12)/2);
%w =nefGaussianSumRV(0.5,-0.7071,1,0.5,0.7071,1);
%Q = evalVariance(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = [2 .5;.5 3];
v = nefGaussianRV(zeros(size(R,1),1),R);
%v = nefUniformRV(-[sqrt(R(1,1)*12)/2;sqrt(R(end,end)*12)/2],[sqrt(R(1,1)*12)/2;sqrt(R(end,end)*12)/2]);
%R = evalVariance(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel=nefEqSystem(f,h,w,v,x0);


T=2e3;% Simulation horizont
MC=1e4;% Number of Monte-Carlo simulations


u=[sin(0:pi/(T-1):pi)];% Input of system

% Four Examples from article: O. Kost, O. Straka, J. Duník, "Identification of state and measurement noise covariance matrices using Nonlinear Estimation Framework", The 12th European Workshop on Advanced Control and Diagnosis, 2015
for example=1:4
    startOld=now;
    for iMC=1:MC
        start=now;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulating trajectory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z = myModel.simulate(T,u);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identify of matrices Q and R using method ALS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % varargin of identifyQR: 'knownQ','knownR','ALSAprioriSetting'
        % example 1
        if(example==1)[wEst,vEst] = myALS.identifyQR(z,u,'ALSAprioriSetting','fast');

        % example 2
        elseif(example==2)[wEst,vEst] = myALS.identifyQR(z,u,'ALSAprioriSetting','automatic');

        % example 3 
        elseif(example==3)
            aprioriMatrix.Q=eye(n)*10;
        	aprioriMatrix.R=eye(p)*6;
            [wEst,vEst] = myALS.identifyQR(z,u,'ALSAprioriSetting',aprioriMatrix);

        % example 4 
        elseif(example==4)
            elemensR = [2 0.5;0.5 NaN];
            [wEst,vEst] = myALS.identifyQR(z,u,'knownR',elemensR,'ALSAprioriSetting','fast');
        end
        MCestW(:,:,iMC)=wEst;
        MCestV(:,:,iMC)=vEst;
        
        if(abs(mod(((iMC)/(MC))*10000,1))<(10^-10) | abs(mod(((iMC)/(MC))*10000,1)-1)<(10^-10))
            disp([datestr((now-startOld)*((MC-iMC)/iMC),'dd:HH:MM:SS'),' (',datestr(now+(now-startOld)*((MC-iMC)/iMC),'dd.mm.yy HH:MM'),') | ',datestr((now-start+1)*60,'HHMM:SS'),' | ',num2str((iMC)/(MC)*100),'%'])
        end
    end
    disp(['Total time of Simulation+Identification of Example ',num2str(example),': ', num2str(datestr((now-startOld),'dd:HH:MM:SS'))]);
    disp(' ')
    
    % saving a calculation data
    if(example==1)save('MC_ACD_I.mat');
    elseif(example==2)save('MC_ACD_II.mat');
    elseif(example==3)save('MC_ACD_III.mat'); 
    elseif(example==4)save('MC_ACD_IV.mat');
        
    end  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluating means, MSEMs and making histograms of MC simulations for all examplex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for example=1:4
    disp(' ')
    if(example==1)
        load('MC_ACD_I.mat');
        disp('Example I');
    elseif(example==2)
        disp('Example II');        
        load('MC_ACD_II.mat');
    elseif(example==3)
        disp('Example III');
        load('MC_ACD_III.mat'); 
    elseif(example==4)
        disp('Example IV');
        load('MC_ACD_IV.mat');
    end    

    wEst=mean(MCestW,3);
    vEst=mean(MCestV,3);
    MSE_Q_Est=sqrt(var(MCestW,0,3))
    MSE_R_Est=sqrt(var(MCestV,0,3))

    disp('true | estimated - Q');
    xx=[' | '];
    for i=2:size(F,1)
      xx=[xx;' | '];
    end
    disp([num2str(Q),xx,num2str(wEst)]) 
    disp(' ')
    disp('true | estimated - R');
    xx=[' | '];
    for i=2:size(H,1)
      xx=[xx;' | '];
    end
    disp([num2str(R),xx,num2str(vEst)])  

    for i=1:iMC-1
        histmc(1,i)=(MCestW(1,1,i));
        histmc(2,i)=(MCestV(1,1,i));
        histmc(3,i)=(MCestV(1,2,i));
        histmc(4,i)=(MCestV(2,2,i));
    end

    figure  
    n=50;
    subplot(1,4,1);
    hist([histmc(1,:)],n)
    hold on
    plot(Q(1,1),0,'*r')
    xlabel('Q')
    ylabel('number of occurrences')
    title(['Histogram of Estimates from Example ', num2str(example)]);
    grid

    subplot(1,4,2)
    hist([histmc(2,:)],n)
    hold on
    plot(R(1,1),0,'*r')
    xlabel('R(1,1)')
    ylabel('number of occurrences')
    grid

    subplot(1,4,3)
    hist([histmc(3,:)],n)
    hold on
    plot(R(1,2),0,'*r')
    xlabel('R(1,2)')
    ylabel('number of occurrences')
    grid

    subplot(1,4,4)
    hist([histmc(4,:)],n)
    hold on
    plot(R(2,2),0,'*r')
    xlabel('R(2,2)')
    ylabel('number of occurrences')
    grid
    legend('estimates','true')
end

