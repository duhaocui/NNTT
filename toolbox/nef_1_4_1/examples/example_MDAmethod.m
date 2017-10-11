%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate of Noise Covariance Matrices by MDA Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear system with unknown noise covariance matricis (Q, R)
% x(k+1) = F(k)*x(k) + N(k)*u(k) + M(k)*w(k)
% z(k) = H*x(k) + O*v(k) / = h(x(k),k) + O(k)*v(k)
% p(w(k)), mean = 0 , covariance = Q
% p(v(k)), mean = 0 , covariance = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used method:
% MDA - Measurement difference autocovariance method
%
% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia
clear all
clc
close all

numberMeasurement=1e4;% min 2*N+2*L-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f in state equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%f = nefHandleFunction(@(x,u,w,t) [x(1);x(2)]+u+[1 0;0 1]*w,[2 2 2 0]);%LTI

f = nefLinFunction([1 0;0 1],[1 1;0 1],[2 1;1 -1]);%LTI
nx=f.dimState;
nq=f.dimNoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h in measurement equation (with 1st derivatives)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zeme=6371*1e3;
druzice=zeme+20350*1e3;

% s1=[-sqrt(0.4);sqrt(0.6)].*druzice;
% s2=[sqrt(0.6);sqrt(0.5)].*druzice;
% s3=[-sqrt(0.7);sqrt(0.4)].*druzice;
% h = nefHandleFunction(@(x,u,v,t)[sqrt((x(1)-s1(1))^2+(x(2)-s1(2))^2);sqrt((x(1)-s2(1))^2+(x(2)-s2(2))^2);sqrt((x(1)-s3(1))^2+(x(2)-s3(2))^2)]+v,[2 0 3 0],...
%   'diff1State',@(x,u,v,k) [(x(1)-s1(1))/sqrt((x(1)-s1(1))^2+(x(2)-s1(2))^2) (x(2)-s1(2))/sqrt((x(1)-s1(1))^2+(x(2)-s1(2))^2);...
%                            (x(1)-s2(1))/sqrt((x(1)-s2(1))^2+(x(2)-s2(2))^2) (x(2)-s2(2))/sqrt((x(1)-s2(1))^2+(x(2)-s2(2))^2);...
%                            (x(1)-s3(1))/sqrt((x(1)-s3(1))^2+(x(2)-s3(2))^2) (x(2)-s3(2))/sqrt((x(1)-s3(1))^2+(x(2)-s3(2))^2)],...
%   'diff1Noise',@(x,u,v,k) eye(3));%NTV


% h = nefHandleFunction(@(x,u,v,t)[2 0;0 1;1 1]*x+v,[2 0 3 0],...
%   'diff1State',@(x,u,v,k) [2 0;0 1;1 1],...
%   'diff1Noise',@(x,u,v,k) eye(3));%LTI

h = nefLinFunction([0.5 0;0 1;1 1],[],[1 1 0;0 2 0;0 0 1]);%LTI the same LTI

nz=size(h,1);
nr=h.dimNoise;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condiditon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = nefGaussianRV([0;zeme],zeros(nx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown state noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = [1 0;0 2];
w = nefGaussianRV(zeros(size(Q,1),1),Q);
%Q = evalVariance(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unknown measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = diag(1:nr);
v = nefGaussianRV(zeros(size(R,1),1),R);
%R = evalVariance(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myModel=nefEqSystem(f,h,w,v,x0);

ux=10*[cos(2*pi*(1:numberMeasurement)/numberMeasurement);sin(2*pi*(1:numberMeasurement)/numberMeasurement)]*pi/numberMeasurement;
%Q = eye(nx);
%w = nefGaussianRV(zeros(size(Q,1),1),Q);
%w=drawSample(w,numberMeasurement);

[z,xt] = myModel.simulate(numberMeasurement,ux);

ux=[ux(:,2:end),zeros(size(ux,1),1)];

L=1;
N=1;
myMDA=nefMDA(f,h,'numberObs',L,'numberPred',N);

xLin=x0.Mean;
xLin=xLin(:,ones(1,numberMeasurement));%'xLin',xLin(:,:)
%xLin=xt;

% a priory Q/R = [1st 2nd 4th ; 2nd 3rd 5th ; 4th 5th 6th];
% aprQR.Q/.R = [1st;2nd;3rd;4th;5th;6th];
% aprQR.Q/.R is upper triangular matrix of apriory matrices Q/R
% aprQR.P is block diag matrix where first matrix is a priory covariance of Q
%                               and second matrix is a priory covariance of R

aprQR.Q=zeros(nr*(nr+1)/2,1);
aprQR.R=zeros(nr*(nr+1)/2,1);
aprQR.P=blkdiag(1e5*eye(nq*(nq+1)/2),1e5*eye(nr*(nr+1)/2));%,'aprioryQR',aprQR


% Q = [1st 2nd 4th ; 2nd 3rd 5th ; 4th 5th 6th]
% knownElementsQ = [1st;2nd;3rd;4th;5th;6th];
% knownElementsQ is upper triangular matrix Q
% Q Elements which will estimated: 1st, 2nd, and 5th
% Known Q elements: 3rd=1 , 4th=0 , and 6th=2
% 'knownQ' - knownElementsQ=[NaN,NaN,1,0,NaN,2]
% matrix R is defined equally

knownElementsQ = [1;NaN;NaN]; %,'knownQ',elemensQ
knownElementsR = [NaN;0;2;NaN;NaN;NaN]; %,'knownR',elemensR

[Qest,Rest] = myMDA.identifyQR(z,ux,'knownQ',knownElementsQ,'knownR',knownElementsR);

for i=1:nr
    for j=i:nr
        plot(squeeze(Rest(i,j,:)),'r')
        hold on
    end
end
figure
for i=1:nq
    for j=i:nq
        plot(squeeze(Qest(i,j,:)),'b')
        hold on
    end
end

meanQEst=mean(Qest(:,:,2*N+L:end),3)
meanREst=mean(Rest(:,:,2*N+L:end),3)

lastQEst=Qest(:,:,end)
lastREst=Rest(:,:,end)

