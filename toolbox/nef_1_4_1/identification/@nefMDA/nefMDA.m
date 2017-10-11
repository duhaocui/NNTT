classdef nefMDA < nefIdentification
   %file @nefMDA/nefMDA.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  % Identification and Decision Making Research Group, Department of Cybernetics,
  % University of West Bohemia
    
  properties (SetAccess = 'protected')
    f=[];% Equation of system, object of nefLinFunction or nefHandleFunction
    h=[];% Equation of measurement, object of nefLinFunction or nefHandleFunction
  
    L;% Number of elements in observability matrix O
    N;% Number of prediction step
  end 
  methods
    function [obj] = nefMDA(f,h,varargin)
      %   nefMDA Creates nefMDA object.
      %
      %   OBJ = nefMDA(f,h,varargin) creates a nefMDA object OBJ 
      %   representing model of system for identification method Measurement difference autocovariance (MDA)
      %
      %   f must be object of nefLinFunction or nefHandleFunction
      %   h must be object of nefLinFunction or nefHandleFunction
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %  PARAMETER           DEFAULT VAL DESCRIPTION                          VALID VALUES
      % ======================================================================================
      %  numberObs           1                                                x > 0 
      %  numberPred          1                                                x > 0
    
      p = inputParser;
      p.KeepUnmatched = false;
      p.FunctionName = 'NEFMDA';
      p.addRequired('f',@(x) isa(x, 'nefHandleFunction') || isa(x, 'nefLinFunction'));
      p.addRequired('h',@(x) isa(x, 'nefHandleFunction') || isa(x, 'nefLinFunction'));
      p.addParameter('numberObs',1,@(x) ~mod(x,1) && x>0 && isnumeric(x));
      p.addParameter('numberPred',1,@(x) ~mod(x,1) && x>0 && isnumeric(x));      
      p.parse( f, h , varargin{:});
      

      obj.L = p.Results.numberObs;
      obj.N = p.Results.numberPred;
      
      obj.f=f;
      obj.h=h;
      
    end
    function [Qest,Rest] = identifyQR(sys,z,u,varargin)
        % Identification of matrices of noises from system and measurement
        %
        % [Q,R] = identifyQR(model,measurement,input)
        %
        % model - must be object of nefMDA
        % measurement - measurement of system, must be matrix relevant dimension
        % input - control of system, must be matrix relevant dimension
        % Q - matrix of identified state noises
        % R - matrix of identified measurement noises              
        %
        % the following user-defined parameters can be changed
        % via the standard Parameter-Value MATLAB mechanism
        %  PARAMETER           DEFAULT VAL DESCRIPTION                                                                   VALID VALUES
        % ======================================================================================================================================
        %  xLin                linearisation points                    matrix size nx x numberOfMeasurements    vector
        %  knownQ              vector full of NaNs                     vector size nx*(nx+1)/2                  vector (with unique elements corresponding the unique elements of covariance matrices) with numbers(known covariance matrix elements) and NaNs(estimated covariance matrix elements)
        %  knownR              vector full of NaNs                     vector size nz*(nz+1)/2                  vector (with unique elements corresponding the unique elements of covariance matrices) with numbers(known covariance matrix elements) and NaNs(estimated covariance matrix elements)
        %  aprioryQR           structure.Q(.R)=0, structure.P=1e5      structure

        q = inputParser;
        q.KeepUnmatched = false;
        q.FunctionName = 'NEFMDAidentify';
        q.addRequired('sys',@(x) isa(x, 'nefMDA'));
        q.addRequired('z',@(x) isnumeric(x));
        q.addRequired('u',@(x) isnumeric(x) || isempty(x));
        q.addParameter('xLin',[],@(x) isnumeric(x));
        q.addParameter('knownQ',[],@(x) isnumeric(x)); 
        q.addParameter('knownR',[],@(x) isnumeric(x));
        q.addParameter('aprioryQR',[],@(x) isstruct(x));
        q.parse(sys, z, u, varargin{:});
        
        nx=sys.f.dimState;
        nu=sys.f.dimInput;
        nq=sys.f.dimNoise;
        nz=size(sys.h,1);
        nr=sys.h.dimNoise;
        L=sys.L;
        N=sys.N;
        eq=L+N;
        numberMeasurement=size(z,2);

        if(isempty(u))
           u=zeros(nu,numberMeasurement); 
        elseif(~(size(u,1)==nu))
           error(['NEF:nefMDA:identify - ',...
           'input must be matrix relevant dimension'])   
        end 
        
        if(~(size(u,2)==size(z,2)))
            error(['NEF:nefMDA:identify - ',...
           'measurement and inputs must be same length']) 
        end
        
        if isempty(q.Results.aprioryQR)
            Papr=1e5*eye(nq*(nq+1)/2+nr*(nr+1)/2);
            Qapr=zeros(nq*(nq+1)/2,1);
            Rapr=zeros(nr*(nr+1)/2,1);
        else             
            Papr=q.Results.aprioryQR.P;
            if(~prod(size(Papr)==[nq*(nq+1)/2+nr*(nr+1)/2,nq*(nq+1)/2+nr*(nr+1)/2]))
               error(['NEF:nefMDA:identify - ',...
               'Apriori P must be matrix relevant dimension'])  
            end
            Qapr=q.Results.aprioryQR.Q;
            if(~prod(size(Qapr)==[nq*(nq+1)/2,1]))
               error(['NEF:nefMDA:identify - ',...
               'Apriori Q must be matrix relevant dimension'])  
            end
            Rapr=q.Results.aprioryQR.R;
            if(~prod(size(Rapr)==[nr*(nr+1)/2,1]))
               error(['NEF:nefMDA:identify - ',...
               'Apriori Q must be matrix relevant dimension'])  
            end
        end
        
        Qk=q.Results.knownQ;
        if(isempty(Qk))
            Qk=NaN*ones(nq*(nq+1)/2,1);
        elseif(~prod(size(Qk)==[nq*(nq+1)/2,1]))
           error(['NEF:nefMDA:identify - ',...
           'Known Q must be matrix relevant dimension'])  
        end
        
        Rk=q.Results.knownR;
        if(isempty(Rk))
            Rk=NaN*ones(nr*(nr+1)/2,1);
        elseif(~prod(size(Rk)==[nr*(nr+1)/2,1]))
           error(['NEF:nefMDA:identify - ',...
           'Known R must be matrix relevant dimension'])  
        end

        if isa(sys.f, 'nefHandleFunction')% LTV or NTV model
            [F,Fu,Fw]=sys.matrixF(sys,numberMeasurement);
        elseif isa(sys.f, 'nefLinFunction')% LTI model
            F=sys.f.F(:,:,ones(1,numberMeasurement));
            Fu=sys.f.G(:,:,ones(1,numberMeasurement));
            Fw=sys.f.H(:,:,ones(1,numberMeasurement));
        end
        
        if isa(sys.h, 'nefHandleFunction')% LTV or NTV model
            
            xLin=q.Results.xLin;
            if(isempty(xLin))
                xLin=zeros(nx,numberMeasurement);
            elseif(~prod(size(xLin)==[nx,numberMeasurement]))
               error(['NEF:nefMDA:identify - ',...
               'Linearisation points must be relevant dimension'])  
            end            
            
            [zL,H,Hv]=nefMDA.linearization(sys,z,xLin);
        elseif isa(sys.h, 'nefLinFunction')% LTI model      
            H=sys.h.F(:,:,ones(1,numberMeasurement));
            Hv=sys.h.H(:,:,ones(1,numberMeasurement));
            zL=z;
        end
        
        isLTI=(isa(sys.f, 'nefLinFunction') && isa(sys.h, 'nefLinFunction'));% model is LTI
 
        O = nefMDA.observabilityMatrix(sys,F,H);% eq 5, Computing "observability" matrix O
        
        Gamma = nefMDA.computationGamma(sys,F,Fw,H,Hv);% eq 7, Computing matrix Gamma 
        GammaU = nefMDA.computationGammaU(sys,F,H,Hv);
        extendedZ = nefMDA.extendZwithoutControl(sys,zL,u,Gamma,GammaU,Fu);% eq 5/7/8, Computing of extended measurement without control
        [MPE,B] = nefMDA.computeMPE(sys,extendedZ,u,F,Fu,O);% eq 9, Computing of MPE and matrix B
        covMPE = nefMDA.computeAutocovarianceMPE(sys,MPE);% eq 11, Computing of MPE covariance
        for i=1+N:numberMeasurement-L+1, G(:,:,i)=B(:,:,i)*Gamma(:,:,i-N);, end% eq 10, Only substitutions        
        [Tqc,Trc] = nefMDA.matrixConnectingQandR(sys);% eq 64, Connecting the same elements of Q and R       
        [Q_Symetric,R_Symetric] = nefMDA.matrixConnectingUniqueQandR(sys);% eq 63, Connecting the unique elements of Q and R
              
        
        allKnownEnlements=[Qk;Rk];                                      % known elements (estimating elemests are NaN)
        estimate=isnan(allKnownEnlements);                              % index of known elements 
        estimateNeg=~isnan(allKnownEnlements);                          % index of UNknown elements 
        QRknown=allKnownEnlements;QRknown(isnan(allKnownEnlements))=0;  % known elements (estimating elemests are 0)
       
        
        %n=(nz*(L+1)+nx*(L+1))^2;       % number of rows matrix A(eq 10)
        p=(nz*L)^2*eq;                  % number of columns matrix A(eq 10)
        ne=sum(estimate);               % number of estimated elements
        

        startTime=2*N+L;              % Start time of method MDA
  
        F12Sym=blkdiag(Q_Symetric,R_Symetric);  % Matrix which connect symetric elements of Q ad R
        ey=eye(nr*(nr+1)/2+nq*(nq+1)/2);EstimteMatrix=ey(:,estimate);   % Matrix which select only estimated elements of Q ad R

        
        % initializing zeros matrices
        P=zeros(ne,ne,numberMeasurement-L+1);     % Covariance of estimates
        QRest=zeros(ne,numberMeasurement-L+1);    % Initial knowledge of estimates     
        Qest=zeros(nq,nq,numberMeasurement-L+1);  % estimates Q
        Rest=zeros(nr,nr,numberMeasurement-L+1);  % estimates R
        
        % set a priori knowledge
        P(:,:,startTime-1)=EstimteMatrix'*Papr*EstimteMatrix;             
        QRapr=[Qapr;Rapr];QRest(:,startTime-1)=QRapr(estimate);
        
        for k=startTime:numberMeasurement-L+1

            
            if (isLTI & k==startTime) | ~isLTI % zrychlení, ještì by šli omezit výpoèty všech matic Gamma,GammaU,B,... 
                A_MDA = nefMDA.computeA_MDA(sys,F,Fw,H,O,B,G,Gamma,k);
                
                % Connecting the same elements of Q and R
                Ar=zeros(p,nr^2);
                Aq=zeros(p,nq^2);
                for i=0:L-1
                    Aq = Aq + A_MDA(:,Tqc(:,i+1));
                    Ar = Ar + A_MDA(:,Trc(:,i+1));
                end          

                FTestSym=[Aq,Ar]*F12Sym; % Connect symetric elements of Q ad R

                FTestSymEstimate=FTestSym*EstimteMatrix; % Sellecting only estimating elements of Q and R

                FTestSymKnown=FTestSym*diag(estimateNeg);% Matrix A of known elements Q and R
                
                if rank(FTestSymEstimate)<ne       
                    error(['NEF:nefMDA:identify - '...
                    ,'can be identified only ',num2str(rank(FTestSymEstimate)-sum(isnan(Rk)))...
                    ,' unique elements of state noise matrix and'...
                    ,', you want estimate ',num2str(sum(isnan(Qk)))...
                    ,' (if you want identified all state noise matrix must'...
                    ,' set more (known) elements of state noise matrix']) 
                end
            end
          
            covMPEEstimate=covMPE(:,k)-FTestSymKnown*QRknown;% Compensation of known elements Q and R           

            % Recursive mean-squares Method
            e=covMPEEstimate-FTestSymEstimate*QRest(:,k-1);
            invL=inv(eye(p)+FTestSymEstimate*P(:,:,k-1)*FTestSymEstimate');
            K=P(:,:,k-1)*FTestSymEstimate'*invL;
            QRest(:,k)=QRest(:,k-1)+K*e;
            P(:,:,k)=P(:,:,k-1)-K*FTestSymEstimate*P(:,:,k-1);
            
            QRestimate=EstimteMatrix*QRest(:,k);% estimated elements Q and R
            QRestimate=QRestimate+QRknown;% estimated elements Q and R + know elements Q and R

            % estimated elemnts
            Qest(:,:,k)=reshape(QRestimate(1:size(Q_Symetric,2))'*Q_Symetric',nq,nq);
            Rest(:,:,k)=reshape(QRestimate(size(Q_Symetric,2)+1:end)'*R_Symetric',nr,nr);
        end  
  end  
     
  end % methods  
  methods(Static)
    
    function [zLinearized,H,Hv] = linearization(sys,z,linearizationPoint)  
        [nz,numberMeasurement]=size(z);
        nr=sys.h.dimNoise;
        nx=sys.f.dimState;
        T2=zeros(nz,nx,numberMeasurement);
        for j=1:numberMeasurement
            T1=sys.h.evaluate(linearizationPoint(:,j),[],zeros(nr,1),j);
            T2(:,:,j)=sys.h.Diff1State(linearizationPoint(:,j));
            
            zLinearized(:,j)=z(:,j)-T1+T2(:,:,j)*linearizationPoint(:,j);

            Hv(:,:,j)=sys.h.Diff1Noise(zeros(nr));
        end
        H=T2;
    end     
    function [F,Fu,Fw] = matrixF(sys,numberMeasurement)  
        nx=sys.f.dimState;
        nu=sys.f.dimInput;
        nq=sys.f.dimNoise;
        for j=1:numberMeasurement
            for i=1:nx
                vectorState=zeros(nx,1);vectorState(i)=1;
                F(:,i,j)=sys.f.evaluate(vectorState,zeros(nu,1),zeros(nq,1),j);
            end
            for i=1:nu
                vectorState=zeros(nu,1);vectorState(i)=1;
                Fu(:,i,j)=sys.f.evaluate(zeros(nx,1),vectorState,zeros(nq,1),j);
            end
            for i=1:nq
                vectorState=zeros(nq,1);vectorState(i)=1;
                Fw(:,i,j)=sys.f.evaluate(zeros(nx,1),zeros(nu,1),vectorState,j);
            end
        end
    end    
    function [S] = MatrixS(nx,nz,L,t)
        abst=abs(t);
        if t==0
            S=eye((L+1)*nz+(L+1)*nx);
        elseif(abst<L+1)
            X=[zeros((L+1-abst)*nz,abst*nz),eye((L+1-abst)*nz);zeros(abst*nz,(L+1)*nz)];
            Y=[zeros((L+1-abst)*nx,abst*nx),eye((L+1-abst)*nx);zeros(abst*nx,(L+1)*nx)];
            S=blkdiag(X,Y);
        else
            S=zeros((L+1)*nz+(L+1)*nx,(L+1)*nz+(L+1)*nx);
        end

        if t<0
            S=S';
        end
    end    
    function O = observabilityMatrix(sys,F,H)
        [nz,nx,numberMeasurement]=size(H);
        nq=sys.f.dimState;
        O=zeros(nz*(sys.L),nx,numberMeasurement-sys.L+1);
        for i=1:numberMeasurement-sys.L+1
            for j=0:sys.L-1
                O(1+j*nz:(j+1)*nz,:,i)=nefMDA.mprod(H,F,i+j,j);
            end  
            if rank(O(:,:,i))<nq  
                error(['NEF:nefMDA:identify - '...
                ,'observability matrix have inadequate rank ('...
                ,num2str(rank(O(:,:,i))),')']) 
            end
        end
    end    
    function Gamma = computationGamma(sys,F,Fw,H,Hv);%
        [nz,nx,numberMeasurement]=size(H);
        nr=size(Hv,2);  
        nq=size(Fw,2);
        Gamma=zeros((sys.L)*nz,(sys.L)*nr+(sys.L)*nq,numberMeasurement-sys.L+1);
        for i=1:numberMeasurement-sys.L+1
            for j=0:sys.L-1
                Gamma(1+j*nz:(j+1)*nz,1+j*nr:(j+1)*nr,i)=Hv(:,:,i+j);
                for k=1:j        
                   Gamma(1+j*nz:(j+1)*nz,(sys.L)*nr+nq+(1+(k-1)*nq:k*nq),i) = nefMDA.mprod(H,F,i+j,j-k)*Fw(:,:,i+k-1);
                end
            end
        end
    end  
    function GammaU = computationGammaU(sys,F,H,Hv);%
        [nz,nx,numberMeasurement]=size(H);
        nr=size(Hv,2);       
        GammaU=zeros((sys.L)*nz,(sys.L)*nr+(sys.L)*nx,numberMeasurement-sys.L);
        for i=1:numberMeasurement-sys.L+1
            for j=0:sys.L-1
                GammaU(1+j*nz:(j+1)*nz,1+j*nr:(j+1)*nr,i)=Hv(:,:,i+j);
                for k=1:j            
                   GammaU(1+j*nz:(j+1)*nz,(sys.L)*nr+nx+(1+(k-1)*nx:k*nx),i) = nefMDA.mprod(H,F,i+j,j-k);
                end
            end
        end
    end  
    function [XYprod] = mprod(X,Y,k,i)   
        XYprod=X(:,:,k);
        for i=1:i
            XYprod=XYprod*Y(:,:,k-i);
        end
    end
    function extendedZ = extendZwithoutControl(sys,z,u,Gamma,GammaU,Fu)
        [nz,numberMeasurement]=size(z);
        nx=sys.f.dimState;
        nq=sys.f.dimNoise;
        nu=sys.f.dimInput;
        nr=sys.h.dimNoise;
        for i=0:sys.L-1
            extendedZ(1+i*nz:(i+1)*nz,:)=z(:,(1+i:numberMeasurement-sys.L+1+i));
        end
        for i=1:numberMeasurement-sys.L+1  
            ux=[];
            for j=i:i+sys.L-2
                ux=[ux;Fu(:,:,j)*u(:,j)];
            end
            extendedZ(:,i)=extendedZ(:,i)-GammaU(:,:,i)*[zeros(nr*(sys.L),1);zeros(nx,1);ux(:)];
        end
    end   
    function [MPE,B] = computeMPE(sys,extendedZ,u,F,Fu,O);  
        [nzE,numberMeasurementE]=size(extendedZ);
        MPE=zeros(nzE,numberMeasurementE);
        for i=1+sys.N:numberMeasurementE
            ux=zeros(nzE,1);
            for j=0:sys.N-1
                ux=ux+nefMDA.mprod(O,F,i,j)*Fu(:,:,i-j-1)*u(:,i-j-1);
            end
            B(:,:,i)=nefMDA.mprod(O,F,i,sys.N)*pinv(O(:,:,i-sys.N));
            MPE(:,i)=extendedZ(:,i)-B(:,:,i)*extendedZ(:,i-sys.N)-ux;
        end
    end  
    function covMPE = computeAutocovarianceMPE(sys,MPE)
        [nze,numberMPE]=size(MPE);
        nx=sys.f.dimState;
        nz=size(sys.h,1);
        covMPE=zeros((nz*(sys.L))^2*(sys.L+sys.N),1);
        for k=2*sys.N+sys.L:numberMPE
            for i=0:sys.L+sys.N-1
               abc=MPE(:,k)*MPE(:,k-i)';
               covMPE(1+i*(nz*(sys.L))^2:(i+1)*(nz*(sys.L))^2,k)=abc(:);
            end
        end
        
    end 
    function [Tqc,Trc] = matrixConnectingQandR(sys)   
        nx=sys.f.dimState;
        nq=sys.f.dimNoise;
        nr=sys.h.dimNoise;
        n=nr*(sys.L)+nq*(sys.L);
        
        Index_Matrix=zeros(n);Index_Matrix(:)=1:n^2;
        for i=0:sys.L-1
            Index_Matrix_R = Index_Matrix(1:nr,1:nr)+(n*nr+nr)*(i);

            Index_Matrix_R(Index_Matrix_R==0)=n;
            Trc(:,i+1)=Index_Matrix_R(:);

            Index_Matrix_Q= Index_Matrix(1:nq,1:nq)+(n*nq+nq)*(i)+(n*nr+nr)*(sys.L);    

            Index_Matrix_Q(Index_Matrix_Q==0)=n;
            Tqc(:,i+1)=Index_Matrix_Q(:);
        end
    end  
    function [A_MDA] = computeA_MDA(sys,F,Fq,H,O,B,G,Gamma,k)
        nx=sys.f.dimState;
        nz=size(sys.h,1);
        nr=sys.h.dimNoise;
        nq=sys.f.dimNoise;
        eq=sys.L+sys.N;
        n=(nr*(sys.L)+nq*(sys.L))^2;
        p=(nz*(sys.L))*(nz*(sys.L));
        
        C=zeros(nq,nr*(sys.L)+nq*(sys.L));
        C(:,nr*(sys.L)+1:nr*(sys.L)+nq)=eye(nq);        
        A_MDA=zeros(p*eq,n+1);%poslední slupec bude nulový
        for j=0:eq-1 
            sum1=zeros(p,n);
            for i=0:sys.N-1
                sum1=sum1+kron(nefMDA.mprod(O,F,k-j,i)*Fq(:,:,k-j-i-1)*C,Gamma(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,i+j));
            end

            sum2=zeros(p,n);
            for i=0:sys.N-1
                sum2=sum2+kron(Gamma(:,:,k-j),nefMDA.mprod(O,F,k,i)*Fq(:,:,k-i-1)*C*nefMDA.MatrixS(nq,nr,sys.L-1,j-i));
            end

            sum3=zeros(p,n);
            for i=0:sys.N-1
                for l=0:sys.N-1
                    sum3=sum3+kron(nefMDA.mprod(O,F,k-j,l)*Fq(:,:,k-j-l-1)*C,nefMDA.mprod(O,F,k,i)*Fq(:,:,k-i-1)*C*nefMDA.MatrixS(nq,nr,sys.L-1,l-i+j));
                end
            end

            sum4=zeros(p,n);
            for i=0:sys.N-1
                sum4=sum4+kron(G(:,:,k-j),nefMDA.mprod(O,F,k,i)*Fq(:,:,k-i-1)*C*nefMDA.MatrixS(nq,nr,sys.L-1,sys.N-i+j));
            end

            sum5=zeros(p,n);
            for i=0:sys.N-1
                sum5=sum5+kron(nefMDA.mprod(O,F,k-j,i)*Fq(:,:,k-j-i-1)*C,G(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,i+j-sys.N));
            end             

            A_MDA(1+j*p:(j+1)*p,1:end-1)=...
             kron(Gamma(:,:,k-j),Gamma(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,j))...
            -kron(G(:,:,k-j),Gamma(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,j+sys.N))...
            -kron(Gamma(:,:,k-j),G(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,j-sys.N))...
            +kron(G(:,:,k-j),G(:,:,k)*nefMDA.MatrixS(nq,nr,sys.L-1,j))...
            +sum1...
            +sum2...
            +sum3...
            -sum4...
            -sum5;
        end      
    end   
    function [Q_Symetric,R_Symetric] = matrixConnectingUniqueQandR(sys) 
        nq=sys.f.dimNoise;
        nr=sys.h.dimNoise;
        
        Select_Index_Matrix = logical(eye(nq^2));
        Index_Matrix = zeros(nq); Index_Matrix(:) = 1:nq^2; Index_Matrix_T = Index_Matrix';% index matrices A and B
        Q_Symetric = Select_Index_Matrix(:,Index_Matrix(triu(ones(nq))~= 0))|Select_Index_Matrix(:,Index_Matrix_T(triu(ones(nq))~=0));

        Select_Index_Matrix = logical(eye(nr^2));
        Index_Matrix = zeros(nr); Index_Matrix(:) = 1:nr^2; Index_Matrix_T = Index_Matrix';% index matrices A and B
        R_Symetric = Select_Index_Matrix(:,Index_Matrix(triu(ones(nr))~=0))|Select_Index_Matrix(:,Index_Matrix_T(triu(ones(nr))~=0));
        
    end
    
  end % static methods
end % classdef
