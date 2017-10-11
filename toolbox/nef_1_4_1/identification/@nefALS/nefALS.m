classdef nefALS < nefIdentification
   %file @nefALS/nefALS.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  % Identification and Decision Making Research Group, Department of Cybernetics,
  % University of West Bohemia
    
  properties (SetAccess = 'protected')
    f=[];% Equation of system, object of nefLinFunction
    h=[];% Equation of measurement, object of nefLinFunction
    Q=[];% Matrix of matrices which are used for calculating the optimum(criterium) gain L in ALS method
    R=[];% Matrix of matrices which are used for calculating the optimum(criterium) gain L in ALS method
    x0;% Initial condiditon of system, object of nefGaussianRV
    
    J=[];% Matrix of values of criteria
    K=[];% Matrix of kalman gains
    index;% Index of optimal(criterium) matrix used in ALS method
  
    aT;% Number of points of criteria
    min;% Begin of points criteria (End of points criteria is 1/min)
    L;% Number of covariance equation used in criterion
    Kappa;% Number of covariance equation used in criterion
  end 
  methods
    function [obj] = nefALS(f,h,x0,varargin)
      %   nefALS Creates nefALS object.
      %
      %   OBJ = nefALS(f,h,x0,varargin) creates a nefALS object OBJ 
      %   representing model of system for identification method Autocovariance Least Squares (ALS)
      %
      %   f must be object of nefLinFunction
      %   h must be object of nefLinFunction
      %   x0 must be object of nefGaussianRV or empty value
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %  PARAMETER           DEFAULT VAL DESCRIPTION                          VALID VALUES
      % ======================================================================================
      %  pointsCrit          100     sample size(integers)                    x > 2 
      %  beginPointCrit      0.001   begin of points criteria                 0 < x < 1
      %  numberCovEq         4       number of covariance equation(integers)  x > 1
      %  omittedInitPeriod   50      share of omitted innovation sequence     0 <= x < 100
    
      p = inputParser;
      p.KeepUnmatched = false;
      p.FunctionName = 'NEFALS';
      p.addRequired('f',@(x) isa(x, 'nefLinFunction'));
      p.addRequired('h',@(x) isa(x, 'nefLinFunction'));
      p.addRequired('x0',@(x) isa(x, 'nefGaussianRV') || isempty(x));
      p.addParameter('pointsCrit'    ,100  ,@(x) ~mod(x,1) && x>2 && isnumeric(x));
      p.addParameter('beginPointCrit',1e-3 ,@(x) validateattributes(x,{'numeric'},{'>',0,'<',1}));      
      p.addParameter('numberCovEq'   ,-1    ,@(x) ~mod(x,1) && x>1 && isnumeric(x));
      p.addParameter('omittedInitPeriod'   ,50    ,@(x)  validateattributes(x,{'numeric'},{'>=',0,'<',100}));
      p.parse( f, h, x0 , varargin{:});

      obj.aT = p.Results.pointsCrit;
      obj.min = p.Results.beginPointCrit;
      obj.Kappa = p.Results.omittedInitPeriod;
      
      obj.f=f;
      if(size(f,1)~=size(h.F,2))
          error(['NEF:nefALS:nefALS - ',...
          'f and h must be matrix relevant dimension'])   
      end
      obj.h=h;
      if(isempty(x0))
          obj.x0 = nefGaussianRV(zeros(size(obj.f,1),1),eye(size(obj.f,1)));
      else
          if(~prod(size(x0)==size(f)))
              error(['NEF:nefALS:nefALS - ',...
              'x0 must be matrix relevant dimension'])   
          end
          obj.x0=x0; 
      end  
      obj.L = p.Results.numberCovEq;
      if(obj.L==-1)
          obj.L=ceil((2*size(f,1)+size(h,1)+1)/(2*size(h,1)));
      end
    end
    function [w,v,sys] = identifyQR(sys,z,u,varargin)
        % Identification of matrices of noises from system and measurement
        %
        % [w,v] = identifyQR(model,measurement,input)
        %
        % model - must be object of nefALS
        % measurement - measurement of system, must be matrix relevant dimension
        % input - control of system, must be matrix relevant dimension
        % w - matrix of identified state noises
        % v - matrix of identified measurement noises              
        %
        % the following user-defined parameters can be changed
        % via the standard Parameter-Value MATLAB mechanism
        %  PARAMETER           DEFAULT VAL DESCRIPTION                                                 VALID VALUES
        % ======================================================================================================================================
        %  ALSAprioriSetting   optimal(criterium) gain    string or vector of objects nefGaussianRV    'fast' or 'automatic' or AprioriQR
        %  knownQ              matrix full of NaNs        matrix size n x n                            matrix with numbers and NaNs
        %  knownR              matrix full of NaNs        matrix size p x p                            matrix with numbers and NaNs        

        q = inputParser;
        q.KeepUnmatched = false;
        q.FunctionName = 'NEFALSidentify';
        q.addRequired('sys',@(x) isa(x, 'nefALS'));
        q.addRequired('z',@(x) isnumeric(x));
        q.addRequired('u',@(x) isnumeric(x) || isempty(x));
        q.addParameter('ALSAprioriSetting',[],@(x) isstruct(x) || isstr(x)); 
        q.addParameter('knownQ',[],@(x) isnumeric(x)); 
        q.addParameter('knownR',[],@(x) isnumeric(x));
        q.parse(sys, z, u, varargin{:});
        
        F=sys.f.F;
        H=sys.h.F;    
        n=size(F,1);
        p=size(H,1);
        nu=(n*n-n)/2+n;
        pu=(p*p-p)/2+p;
        Kappa=sys.Kappa;
        L=sys.L; 

        if(rank(obsv(F,H))~=n)
            error(['NEF:nefALS:identify - ',...
           'system must be observable'])  
        end             
        if(~(size(z,1)==p))
           error(['NEF:nefALS:identify - ',...
           'measurement must be matrix relevant dimension'])   
        end
        if(isempty(u))
           u=zeros(n,size(z,2)); 
        elseif(~(size(u,1)==n))
           error(['NEF:nefALS:identify - ',...
           'input must be matrix relevant dimension'])   
        end        
        if(~(size(u,2)==size(z,2)))
            error(['NEF:nefALS:identify - ',...
           'measurement and inputs must be same length']) 
        end
        
        if((size(z,2)*(1-Kappa/100))<L)
            error(['NEF:nefALS:identify - ',...
           'Kappa must be bigger (shortened inovation sequence is overly short)'])  
        end
        
        if([1,4]==size(q.Results.ALSAprioriSetting)) 
            if('fast'~=(q.Results.ALSAprioriSetting)) 
                error(['NEF:nefALS:identify - ',...
               'ALSAprioriSetting have incorrect parameters']) 
            end
            sc=1;
        elseif(isstruct(q.Results.ALSAprioriSetting))
            if(~issymmetric(q.Results.ALSAprioriSetting.Q) && ~issymmetric(q.Results.ALSAprioriSetting.R)) 
                error(['NEF:nefALS:identify - ',...
               'Apriori Setting matrices Q and R must be symmetric'])
            elseif(~prod((abs(eig(q.Results.ALSAprioriSetting.Q))>=0)) && ~prod((abs(eig(q.Results.ALSAprioriSetting.R))>=0))) 
               error(['NEF:nefALS:identify - ',...
               'Apriori Setting matrices Q and R must be positive semidefinite']) 
            end
            wv = q.Results.ALSAprioriSetting;
            wa = nefGaussianRV(zeros(n,1),wv.Q);
            va = nefGaussianRV(zeros(p,1),wv.R);
            sc=2;           
        elseif(prod([1,9]==size(q.Results.ALSAprioriSetting)) || isempty(q.Results.ALSAprioriSetting))
            if(~isempty(q.Results.ALSAprioriSetting))
                if(prod('automatic'~=(q.Results.ALSAprioriSetting)) && isempty(q.Results.ALSAprioriSetting))
                    error(['NEF:nefALS:identify - ',...
                   'ALSAprioriSetting have incorrect parameters']) 
                end
            end
            sc=0;
        else
            error(['NEF:nefALS:identify - ',...
            'ALSAprioriSetting have incorrect parameters'])             
        end
               
        Qa=q.Results.knownQ;
        if(isempty(Qa))
            Qa=NaN*ones(n,n);
        elseif(~prod(size(Qa)==[n,n]))
           error(['NEF:nefALS:identify - ',...
           'Apriori Q must be matrix relevant dimension'])  
        end
        for i=1:n 
            for j=1:n
               if(isnan(Qa(i,j))) 
                   Qax(i,j)=0;
               else
                   Qax(i,j)=Qa(i,j);
               end
            end
        end
        if(sum(sum(Qax~=Qax')) || sum(sum(isnan(Qa)==isnan(Qa)'))~=n*n)
            error(['NEF:nefALS:identify - ',...
           'Apriori Q must be symetric matrix']) 
        end
        
        Ra=q.Results.knownR;
        if(isempty(Ra))
            Ra=NaN*ones(p,p);
        elseif(prod(size(Ra)~=[p,p]))
           error(['NEF:nefALS:identify - ',...
           'Apriori R must be matrix relevant dimension'])  
        end
        for i=1:p 
            for j=1:p
               if(isnan(Ra(i,j))) 
                   Rax(i,j)=0;
               else
                   Rax(i,j)=Ra(i,j);
               end
            end
        end
        if(sum(sum(Rax~=Rax')) || sum(sum(isnan(Ra)==isnan(Ra)'))~=p*p)
            error(['NEF:nefALS:identify - ',...
           'Apriori R must be symetric matrix']) 
        end

        Qa=nefALS.reductionMatrix(Qa);
        Ra=nefALS.reductionMatrix(Ra);
        noiseAp=[Qa;Ra];         

        % 1 - Calculate state prediction
        if(sc==0)% ALSAprioriSetting = 'automatic'  
            K=nefIdentification.gainCalculation(sys.f,sys.h,nefGaussianRV(zeros(n,1),eye(n)),nefGaussianRV(zeros(p,1),eye(p)));
            An = rank(nefALS.reductionLSDesignMatrix(sys,nefALS.calculateLSDesignMatrix(sys,K)))-pu;
            if(An<sum(sum(isnan(Qa))))        
                warning(['NEF:nefALS:identify - '...
                ,'can be identified only ',num2str(An)...
                ,'  unique elements of matrix state noises and'...
                ,', you want ',num2str(sum(sum(isnan(Qa))))...
                ,' (if you want identified all matrix must'...
                ,' set more elements Q [knownQ])']) 
                Qa=diag(NaN*ones(n,1));
                Qa=nefALS.reductionMatrix(Qa);
                noiseAp=[Qa;Ra];
            end 
            % I - Criterium calculation
            [sys.Q,sys.R,sys.J,sys.K]=nefALS.calculateCrit(sys,noiseAp);
            
            % II - Minimizing criterium 
            [w,v,K,sys.index]=nefALS.minimalizationCrit(sys);
            
            % III - Simulating xpred for optimal Q and R
            system=nefEqSystem(sys.f,sys.h,w,v,sys.x0);
            KALMAN = nefKalman(system);
            [val_KALMANfil,val_KALMANpred] = KALMAN.estimate(z,u);
            for i = 1:size(z,2)
                xpred(:,i) = evalMean(val_KALMANpred{i});
            end 
        elseif(sc==1)% ALSAprioriSetting = 'fast'  
            if(~prod(abs(eig(F))<1))
                error(['NEF:nefALS:identify - ',...
                'Matrix F must be stable for choice ''fast'' '])             
            end
            sys.Q=zeros(size(F,1));   
            sys.R=eye(size(H,1));           
            K=zeros(size(F,1),size(H,1)); 
            sys.K=K; 
            
%             w = nefGaussianRV(zeros(n,1),sys.Q);
%             v = nefGaussianRV(zeros(p,1),sys.R);           
%             system=nefEqSystem(sys.f,sys.h,w,v,sys.x0);
%             KALMAN = nefKalman(system);
%             predEst.RV=sys.x0;           
%             pred=predEst;    
%             xpred=sys.x0.Mean;
%             for i=2:size(z,2) 
%               pred=KALMAN.timeUpdate(pred,u(i-1),i);
%               xpred(:,i) = evalMean(pred.RV);                                            
%             end
            
            pred=sys.x0.Mean;
            for i=2:size(z,2)          
              pred=nefKalman.KalmanTimeUpdate(pred,sys.x0.Var,sys.f,u(:,i-1),zeros(n,1),sys.Q,i,F,sys.f.H);                         
              xpred(:,i)=pred;                                  
            end

        elseif(sc==2)% ALSAprioriSetting = aprioriQR   
            system=nefEqSystem(sys.f,sys.h,wa,va,sys.x0);
            KALMAN = nefKalman(system);
            [val_KALMANfil,val_KALMANpred] = KALMAN.estimate(z,u);
            for i = 1:size(z,2)
                xpred(:,i) = evalMean(val_KALMANpred{i});
            end 
            K = nefIdentification.gainCalculation(sys.f,sys.h,wa,va);
            sys.K=K;
            sys.R=va.Var;
            sys.Q=wa.Var;
        end    
        
        % 2 - Calculate matrix LSDM (LSDM*x(Q,R)=bnc)
        LSDM = nefALS.calculateLSDesignMatrix(sys,K); 
        reducLSDM = nefALS.reductionLSDesignMatrix(sys,LSDM);      
        
        reducLSDMx = rank(reducLSDM)-pu; % Control, if all elements matrix x are estimable
        if(reducLSDMx<sum(sum(isnan(Qa))))        
            warning(['NEF:nefALS:identify - '...
            ,'can be identified only ',num2str(An)...
            ,'  unique elements of matrix state noises and'...
            ,', you want ',num2str(sum(sum(isnan(Qa))))...
            ,' (if you want identified all matrix must'...
            ,' set more elements Q [knownQ])'])           
            Qa=diag(NaN*ones(n,1));
            Qa=nefALS.reductionMatrix(Qa);
            noiseAp=[Qa;Ra];
        end 
        
        % 3 - Calculate matrix bnc (Anc*x(Q,R)=bnc)
        [RN] = nefALS.calcInnovationCov(sys,xpred,z);            
        bnc = reshape(RN,(sys.L*p^2),1);
        
%         if(false)% true matrix bnc (Anc*x(Q,R)=bnc)
%             Qt=[1];% true Q
%             Rt=[1];% true R
%             Ap = F - F*K*H;%(1)
%             Oo=H;%(7)
%             for i=1:sys.aN-1
%                 Oo = [Oo;H*Ap^(i)];%(7)
%             end
%             gamma=eye(p);%(7)
%             for i=1:sys.aN-1
%                 gamma=[gamma;-H*Ap^(i-1)*F*K];
%             end
%             %nef xk+1 = F*xk + G*uk + H*wk;
%             %
%             %ACD xk+1 = F*xk + N*uk + M*wk;
%             %    zk   = H*xk        + O*vk;
%             M=sys.f.H;
%             O=sys.h.H;
%             RNteor=(kron(H,Oo)*pinv(eye(n^2)-kron(Ap,Ap)))*reshape(M*Qt*M',n^2,1)+(kron(H,Oo)*pinv(eye(n^2)-kron(Ap,Ap))*kron(F*K,F*K)+kron(eye(p),gamma))*reshape(O*Rt*O',p^2,1);
%             bncTeor=reshape(RNteor,p^2*sys.aN,1);
%             %bnc=bncTeor;
%         end
     
        % 4 - Estimate matrices Q and R (LSDM*x(Q,R)=bnc)
        if(isnan(noiseAp)==ones(size(reducLSDM,2),1))
            QRest=pinv(reducLSDM)*bnc;                
        else
            xindex=[];
            Azvol=[];
            for j=1:length(noiseAp)  
               if(isnan(noiseAp(j))) 
                   xindex=[xindex;0];
                   noiseAp(j)=0;
                   Azvol=[Azvol,reducLSDM(:,j)];
               else
                   xindex=[xindex;1];
               end
            end
            bzvol=reducLSDM*noiseAp;                  
            xest=pinv(Azvol)*(bnc-bzvol);   
            QRest=[];
            k=1;
            for j=1:length(noiseAp)
               if(xindex(j)) 
                   QRest=[QRest;noiseAp(j)];
               else  
                   QRest=[QRest;xest(k)];
                   k=k+1;
               end
            end
        end    
        [Qe,Re]=nefALS.sortMatrix(sys,QRest);
        
        % Control, if matrices Q and R are positive semidefinite
        if(~prod(eig(Qe) >= 0))
            disp(' ');
            disp('Identified Q:');
            disp(num2str(Qe))
            error(['NEF:nefALS:identify - ',...
            'Identified matrix Q is not positive semidefinite'])             
        else
            %w = nefGaussianRV(zeros(n,1),Qe); 
            w = Qe;
        end 
        if(~prod(eig(Re) >= 0))
             disp(' ');
             disp('Identified R:');
             disp(num2str(Re))
             error(['NEF:nefALS:identify - ',...
            'Identified matrix R is not positive semidefinite'])
        else
            %v = nefGaussianRV(zeros(p,1),Re); 
            v = Re;
        end
       
    end      
  end % methods  
  methods(Static)
    function [Qx,Rx,Jx,Kx] = calculateCrit(sys,noiseAp)
        % Calculating of criterim for ALS method
        %
        %   [w,v,L,index] = calculateCrit(model,noiseAp) 
        %
        %   model - must be object of nefALS
        %   noiseAp - apriory set noises
        %
        %   Qx - covariance of state noises use in criterium
        %   Rx - covariance of measurement noises use in criterium
        %   Jx - value of criteria
        %   Kx - Kalman gain

        F=sys.f.F;
        M=sys.f.H;
        H=sys.h.F;
        O=sys.h.H;
        n=size(F,1);
        p=size(H,1);
        nu=(n*n-n)/2+n;
        pu=(p*p-p)/2+p;
        
        if(isnan(noiseAp)==ones(pu+nu,1))
            setJ=eye(n*n+p*p);
        else          
            [Qe,Re] = nefALS.sortFnc(sys,noiseAp);
            setJ=diag([reshape(isnan(Qe),n*n,1);reshape(isnan(Re),p*p,1)]);
        end
        
        Jx=zeros(n*n+p*p,n*n+p*p,sys.aT);
        Qx=zeros(n,n,sys.aT);
        Rx=zeros(p,p,sys.aT);
        Kx=zeros(n,p,sys.aT);
        L=sys.L;
        rozsahQ=(0:(1/sys.min-sys.min)/(sys.aT-1):1/sys.min);  
        rozsahR=(sys.min:(1/sys.min-sys.min)/(sys.aT-1):1/sys.min);  
        for i=1:sys.aT     
            Qx(:,:,i)=eye(n)*rozsahQ(i);
            Q=Qx(:,:,i);
            Rx(:,:,i)=eye(p)*rozsahR(end-i+1);
            R=Rx(:,:,i);          
            K = nefIdentification.gainCalculation(sys.f,sys.h,nefGaussianRV(zeros(n,1),Q),nefGaussianRV(zeros(p,1),R));
            Kx(:,:,i)=K;
            
            Fp = F - F*K*H;
            G = [M, -F*K*O];
            Anc = nefALS.calculateLSDesignMatrix(sys,K);

            Jt = zeros(p^2,n^2+n*p+n);
            k=0;
            for l=1:2*p-1
                if(mod(l,2)>0.5)
                    Jtp=zeros(p^2,p);        
                    for j=1:p
                        Jtp(j+k*p,j)=1;
                    end
                    Jt=[Jt,Jtp];
                    k=k+1;
                else
                    Jt=[Jt,zeros(p^2,n)];
                end
            end
            Mkrit=kron(H,H)*inv(eye(n^2)-kron(Fp,Fp))*kron(G,G)+kron(O,O)*Jt;

            %JN=ones(aN);
            %Jx(:,:,i)=pinv(Anc)*kron((M*M'),JN)*pinv(Anc)';          

            Pc=[];
            k=0;
            for l=1:p*p
                if(mod(l,2)==1)
                    k=k+1;
                    x=[zeros(L,p*p*L-(p*p*L-L*(k-1))),eye(L),zeros(L,p*p*L-L*k)];
                else
                    x=[zeros(L,round(p*p/2)*L+L*(k-1)),eye(L),zeros(L,fix(p*p/2)*L-L*k)];
                end
                Pc=[Pc;x];
            end
            Pr=Pc';

            JN=ones(L);
            Jx(:,:,i)=setJ*pinv(Anc)*Pr*kron(JN,(Mkrit*Mkrit'))*Pc*pinv(Anc)'*setJ;
        end    
    end   
    function [w,v,K,index] = minimalizationCrit(sys)
        %   Minimalization of criterim for ALS method
        %
        %   [w,v,K,index] = minimalizationCrit(model) 
        %
        %   model - must be object of nefALS
        %
        %   w - objekt nefGaussianRV, state noises
        %   v - objekt nefGaussianRV, measurement noises
        %   K - optimal(criterum) gain
        %   index - index of optimal(criterium) values
          
        min=trace(sys.J(:,:,1)); 
        index=1;
        for i=2:sys.aT;
            if(min>trace(sys.J(:,:,i)))
                min=trace(sys.J(:,:,i)); 
                index=i;
            end
        end
     
        w = nefGaussianRV(zeros(size(sys.f.F,1),1),sys.Q(:,:,index));
        v = nefGaussianRV(zeros(size(sys.h.F,1),1),sys.R(:,:,index));
        K = sys.K(:,:,index);  
        
    end   
    function [LSDM] = calculateLSDesignMatrix(sys,K)
        %   Calculating of matrices which are used in the method of ALS
        %
        %   [LSDM] = calculateLSDesignMatrix(model,Gain) 
        %
        %   model - must be object of nefALS
        %   Gain - Kalman gain      
        
        F=sys.f.F; % xk+1 = F*xk + G*uk + H*wk;
        M=sys.f.H;
        H=sys.h.F;
        O=sys.h.H;
        n=size(F,1);
        p=size(H,1);
        L=sys.L;
        Fp = F - F*K*H;
        Oo=H;
        for i=1:L-1
            Oo = [Oo;H*Fp^(i)];
        end
        D = kron(H,Oo)*pinv(eye(n^2)-kron(Fp,Fp));
        gamma=O;%(7)
        for i=1:L-1
            gamma=[gamma;-H*Fp^(i-1)*F*K*O];
        end
        % calculating Fnc (Fnc*X=bnc, where X=[Q;R])
        Aa1=D*kron(M,M);
        Aa2=D*kron(F*K*O,F*K*O)+kron(O,gamma);
        LSDM=[Aa1 Aa2];      
    end 
    function [RN] = calcInnovationCov(sys,x,z)
        %   Calculating covariance matrix of innovation sequence, which is used in the method of ALS
        %
        %   [RN] = calcInnovationCov(model,prediction,measurement) creates a nefALS object OBJ 
        %
        %   model - must be object of nefALS
        %   prediction - prediction of state of system       
        %   measurement - measurement of system
        %   RN - covariance matrix of innovation sequence     
        n=size(sys.f.F,1);
        p=size(sys.h.F,1);
        E = z-sys.h.F*x;
        Kappa = sys.Kappa;
        E = E(:,(fix(end*(Kappa/100))+1:end));
        L=sys.L;
        N=size(E,2);
        RN=[];
        for i=1:L   
            RNi=(E(:,i:N)*E(:,1:N-i+1)')/(N-i+1);
            RN=[RN;RNi]; 
        end;
    end 
    function [reducLSDM] = reductionLSDesignMatrix(sys,LSDM)
        %   Reduction of matrices which are used in the method of ALS for faster calculation
        %
        %   [reducLSDM] = reductionLSDesignMatrix(model,LSDM)
        %
        %   model - must be object of nefALS
        %   LSDM - matrix which is reduction
        %   reducLSDM - reduced matrix
        n=size(sys.f.F,1);
        p=size(sys.h.F,1);
        
        for i=1:n%reduction number of duplicate columns in Fnc because of inversion
            for j=i:n%Q part
                if(i==1 & j==1)
                   reducLSDM=LSDM(:,1);
                elseif(i==j)
                   reducLSDM=[reducLSDM,(LSDM(:,(j+(i-1)*n)))]; 
                else
                   reducLSDM=[reducLSDM,(LSDM(:,(j+(i-1)*n))+LSDM(:,(i+(j-1)*n)))];
                end
            end
        end
        idxReducLSDM=j+(i-1)*n;
        for i=1:p%R part
            for j=i:p
                if(i==j)
                   reducLSDM=[reducLSDM,(LSDM(:,(j+(i-1)*p+idxReducLSDM)))];      
                else
                   reducLSDM=[reducLSDM,(LSDM(:,(j+(i-1)*p+idxReducLSDM))+LSDM(:,(i+(j-1)*p+idxReducLSDM)))];
                end
            end
        end    
    end
    function [Xx] = reductionMatrix(X)
        %   Reduction of matrix on the unique elements 
        %
        %   [Xx] = reductionMatrix(X)
        %
        %   X - must be square matrix
        %   Xx - vector of unique elements
        n=size(X,1);
        Xx=X(find(tril(ones(n))));
    end 
    function [Qe,Re] = sortMatrix(sys,QRest)
        %   Sorting of unique elements matrix back in all matrix 
        %
        %   [Qe,Re] = sortMatrix(model,Vector)
        %
        %   model - must be object of nefALS
        %   Vector - vector of unique elements two matrices
        %   Qe - first matrices of Vector 
        %   Re - second matrices of Vector
        n=size(sys.f.F,1);
        p=size(sys.h.F,1);
        nu=(n*n-n)/2+n;
        pu=(p*p-p)/2+p;    
        tmp = triu(ones(n));
        tmp(tmp>0) = QRest(1:nu);
        Qe=(tmp+tmp')./(eye(n)+1);
        tmp = triu(ones(p));
        tmp(tmp>0) = QRest(nu+1:end);
        Re=(tmp+tmp')./(eye(p)+1);
    end   
  end % static methods
end % classdef
