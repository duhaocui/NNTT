classdef nefSUKF < nefEstimator
  %file @nefSUKF/nefSUKF.m (Square-Root Unscented Kalman Estimator)
  % nefSUKF Properties:
  %   nefUKF_optsDef - definition of the nefSUKF options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefSUKF Methods:
  %   nefSUKF - class constructor (For list of options type 'help nefSUKF/nefSUKF')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterCheck - check of chosen scaling parameter and warning spec.
  %   smsp - computation of set of deterministically chosen weighted points
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   SUKFTimeUpdate - time update of the SUKF
  %   SUKFMeasurementUpdate - measurement update of the SUKF
  %   SURTSSmoothUpdate - smoothing based on the SUKF and Rauch-Tung-Striebel smoother
  %
  % References:
  % S. J. Julier, J. K. Uhlmann and H. F. Durrant-White (2000):
  % A new method for the nonlinear transformation of means and covariances in filters and estimators.
  % IEEE Trans. On AC 45(3), 477-482.
  %
  % M. Simandl and J. Dunik (2009):
  % Derivative-Free Estimation Methods: New Results and Performance Analysis
  % Automatica 45(7), 1749-1757.
  %
  % M. Simandl, and J. Dunik (2005):
  % Sigma point gaussian sum filter design using square root unscented filters.
  % In: Preprints of 16th IFAC World Congress. Prague, Czech Republic.
  %
  % M. Simandl and J. Dunik (2006):
  % Design of derivative-free smoothers and predictors.
  % In: Preprints of the 14th IFAC Symposium on System Identification.
  % Newcastle, Australia.
  %
  % J. Dunik, M. Simandl and O. Straka (2010): 
  % Adaptive choice of scaling parameter in derivative-free local filters. 
  % In: Proceedings of the 13th International Conference on Information Fusion, Edinburgh, UK.

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % definition of estimator options    
    nefSUKF_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
    2,'scalingParameterType:constant','parameterValue',@(x) max(3-x.sys.dimState,0),'constant scaling parameter', @(x) isfloat(x) && length(x)==1 && x>=0;
    2,'scalingParameterType:adaptive','parameterDomain',[0 4 1],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
    2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE','MMPE'}));};
  end % constant properties

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%
    % MATRICES
    %%%%%%%%%%%%%%%%%%%%%
    Sr = [];
    Sq = [];
    wSVar = [];
    vSVar = [];
    %%%%%%%%%%%%%%%%%%%%%%%
    %PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    %kappa;
    isRConst;
    isQConst;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSUKF(system,varargin)
      % NEFSUKF Creates NEFSUKF object.
      %
      %   OBJ = NEFSUKF(system,varargin) creates a NEFSUKF object OBJ 
      %   representing square-root unscented Kalman filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed 
      %   via the standard Parameter-Value MATLAB mechanism
      %  
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  scalingParameterType        constant         scaling parameter can be             constant, adaptive
      %                                               constant or adaptively chosen   
      %  parameterValue              3-dim(x),0       constant scaling parameter           real, non-negative
      %  parameterDomain             [0 4 1]          adaptive grid search of scaling      nonnegative real, Kmin<Kmax
      %                                               parameter K within domain 
      %                                               [Kmin, Kmax, Kstep] for 
      %                                               adaptive setting
      %  parameterCriterion          MLE              selection of criterion for           MLE, MMPE
      %                                               adaptive search
      %
      %  See also ESTIMATE.      

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFSUKF';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for  i = 1:size(nefSUKF.nefSUKF_optsDef,1)
        p.addParamValue(nefSUKF.nefSUKF_optsDef{i,3},[],nefSUKF.nefSUKF_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefSUKF.nefSUKF_optsDef,p.Results);

      % scaling parameter check
      obj = nefSUKF.ScalingParameterCheck(obj);


      % change x0 to x0 with S factored covariance matrix
      obj.x0 = nefGaussianRV(obj.x0.evalMean,nefCholFactorFunction(evalVariance(obj.x0),'matrixType','cov'),'check',0);

      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end
      obj.estQueue(end+1) = data;

      obj.Sq = [];
      obj.Sr = [];
      obj.wSVar = [];
      obj.vSVar = [];

      obj.isRConst =  obj.sys.v.VarNumeric;
      obj.isQConst =  obj.sys.w.VarNumeric;
    end % NEFUKF constructor

    function disp(obj)
      % DISP Displays nefUKF object parameters and options.

      fprintf('A nefUKF object with parameters\n')
      showOptions(obj)
    end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time,varargin)
      % TIMEUPDATE Performs time update step with respect to system description.

      p = inputParser;
      p.addParamValue('smoothingPurpose',0,@(x) (x==0) || (x==1) );
      p.parse(varargin{:});

      % parameters of filtering or previous prediction estimate
      filtMean = evalMean(filtEst.RV,[],[],[],[]);
      filtVar = filtEst.RV.Var;
      %filtVar = evalVariance(filtEst.RV,[],[],[],[]);      
      % mean and variance of state noise
      wMean = evalMean(obj.sys.w,[],Input,[],Time);
      wVar = evalVariance(obj.sys.w,[],Input,[],Time);

      % obtain S factor of Q = wVar
      % for constant Q,  the factor is computed at the first step only
      % for t-variant Q,  the factor is computed at each time step
      if ((obj.isQConst == 1) && (isempty(obj.wSVar))) || (obj.isQConst == 0)
        obj.wSVar = nefCholFactorFunction.CholFactor(wVar);
      end

      % evaluating predictive mean and variance
      % - if state eq. is linear (with additive noise), then the predictive part of the square-root KF is used
      if obj.sys.f.isLinear
        F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
        Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
        [predMean,predVar] = nefSKalman.SKalmanTimeUpdate(filtMean,filtVar.S,obj.sys.f,Input,wMean,obj.wSVar,Time,F,Gamma);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen - value from the last filtering step
          aux_domain = getOption(obj,'parameterDomain');
          kappa = aux_domain(4);
        end
        [predMean, predVar] = nefSUKF.SUKFTimeUpdate(filtMean,filtVar.S,obj.sys.f,Input,wMean,obj.wSVar,Time,kappa);
      end

      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);

      % is timUpdate is used for smoothing
      if p.Results.smoothingPurpose
        predEstimate.auxData.Input = Input;
        predEstimate.auxData.Time = Time;
        predEstimate.auxData.wMean = wMean;
        predEstimate.auxData.wVar = wVar;
        if ~obj.sys.f.isLinear && strcmp(getOption(obj,'scalingParameterType'),'adaptive')
          predEstimate.auxData.kappa = kappa;
        end
        if obj.sys.f.isLinear
          predEstimate.auxData.F = F;
          predEstimate.auxData.Gamma = Gamma;
        end
      end % if

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
      % MEASUREMENTUPDATE Performs measurement update step with respect to system description.

      % parameters of prediction estimate
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = predEst.RV.Var;
      %predVar = evalVariance(predEst.RV,[],[],[],[]);      
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);

      % obtain S factor of R = vVar
      % for constant R,  the factor is computed at the first step only
      % for t-variant R,  the factor is computed at each time step
      if ((obj.isRConst == 1) && (isempty(obj.vSVar))) || (obj.isRConst == 0)
        obj.vSVar = nefCholFactorFunction.CholFactor(vVar);
      end

      % evaluating filtering mean and variance
      % - if measurement eq. is linear (with additive noise), then the predictive part of the square-root KF is used
      if obj.sys.h.isLinear
        H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
        Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
        [filtMean,filtVar] = nefSKalman.SKalmanMeasurementUpdate(predMean,predVar.S,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,H,Delta);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen
          aux_domain = getOption(obj,'parameterDomain');
          aux_criterion = getOption(obj,'parameterCriterion');
          kappa = nefSUKF.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,aux_domain(1:3),aux_criterion);
          setOption(obj,'parameterDomain',[aux_domain(1:3), kappa]);
        end
        [filtMean,filtVar] = nefSUKF.SUKFMeasurementUpdate(predMean,predVar.S,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,kappa);
      end

      % returning nefGaussianRV with symmetrized variance
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
      % return also by-products if required
    end % function measurementUpdate


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SMOOTHUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothEstimate] = smoothUpdate(obj,initEstimate,filtEstimate,predEstimate)
      % SMOOTHUPDATE Performs smoothing step with respect to system
      % description.

      initM = evalMean(initEstimate.RV);
      %initV = evalVariance(initEstimate.RV);
      initVfac = initEstimate.RV.Var.S;
      filtM = evalMean(filtEstimate.RV);
      %filtV = evalVariance(filtEstimate.RV);
      filtVfac = filtEstimate.RV.Var.S;
      predM = evalMean(predEstimate.RV);
      %predV = evalVariance(predEstimate.RV);
      predVfac = predEstimate.RV.Var.S;
      Input = predEstimate.auxData.Input;
      Time = predEstimate.auxData.Time;
      wMean = predEstimate.auxData.wMean;
      wVar = predEstimate.auxData.wVar;

      % obtain S factor of Q = wVar
      % for constant Q,  the factor is computed at the first step only
      % for t-variant Q,  the factor is computed at each time step
      %if ((obj.isQConst == 1) && (Time == 0)) || (obj.isQConst == 0)
      factor = nefCholFactorFunction(wVar,'matrixType','cov');
      obj.Sq = factor.S;
      %end

      % evaluating smoothing mean and variance
      % - if state eq. is linear (with additive noise), then the predictive part of the square-root KF is used
      if obj.sys.f.isLinear
        F = predEstimate.auxData.F;
        Gamma = predEstimate.auxData.Gamma;
        [smoothM,smoothV] = nefSKalman.SKalmanRTSSmoothUpdate(initM,initVfac,filtM,filtVfac,predM,predVfac,F,Gamma,wMean,obj.Sq,obj.sys.f,Time,Input);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen - stored from the prediction
          kappa = predEstimate.auxData.kappa;
        end
        [smoothM,smoothV] = nefSUKF.SURTSSmoothUpdate(initM,initVfac,filtM,filtVfac,predM,predVfac,wMean,obj.Sq,obj.sys.f,Time,Input,kappa);
      end

      % returning nefGaussianRV with symmetrized variance
      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate

  end % methods
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ScalingParameterCheck
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = ScalingParameterCheck(obj)
      % SCALINGPARAMETERCHECK Checks the chosen scaling parameter against
      % the system description.

      if strcmp(getOption(obj,'scalingParameterType'),'constant')
        % cannot be negative
      else
        % if adaptive choice of parameter, then create a fourth parameter
        % in the "domain" for saving the actual value                
        aux = getOption(obj,'parameterDomain');
        aux = [aux, aux(1)];        
        setOption(obj,'parameterDomain',aux)        
        % if function in measurement eq. is linear then it is independent
        % of the scaling parameter and optimisation is useless        
        if obj.sys.h.isLinear
          disp('**************')
          disp('nefSUKF WARNING: Function in measurement equation is linear and independent of the scaling parameter. Its optimisation is useless.')
          %disp('Better choice of the scaling parameter is "constant" (instead of adaptive).')
          disp('**************')
        end                
      end
    end % function ScalingParameterCheck

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SigmaPointComputation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [chi,w] = smsp(x, S, kappa)
      %   SMSP Method smsp for the calculation of sigma-points, called from methods for prediction, filtration and smoothing.
      %   For standard (non-square-root) estimation algorithms.
      %
      %   function [ch,w]=msp(x,S,kappa)
      %   input parameters                   - x - mean value,
      %                                      - S - covariance matrix factor,
      %                                      - kappa - parameter for distribution of sigma-points,
      %   output parameters                  - chi - calculated point
      %                   		             - w -  weight of the point

      nx = size(x,1);

      pom = sqrt(nx+kappa)*S;

      chi(:,1) = x;
      chi(:,2:2*nx+1) = x(:,ones(1,2*nx)) + [pom, -pom];

      w(1)=kappa/(nx+kappa);
      w(2:2*nx+1)=0.5/(nx+kappa);% Nonlinear Filtering Toolbox version 3.0


      % check
      %chi_ = chi - repmat(x,1,2*nx+1);
      %chi_w = chi_.*w(ones(1,size(x,1)),:);
      %Pap = chi_w*chi_';

    end % function smsp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scaling Parameter Adaptive Choice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function kappa_act = ScalingParameterAdaptiveChoice(predMean,predVar,h,Input,vMean,Sr,Time,Measurement,kdomain,kcriterion)
      % SCALINGPARAMETERADAPTIVECHOICE Performs adaptive choice of the
      % scaling parameter within predefined domain in filtering step.

      % number of sigma points
      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);
      if additive
        Jns = 2*(nx) + 1;
      else
        Jns = 2*(nx+nv) + 1;
      end

      % grid points
      kappa = kdomain(1):kdomain(3):kdomain(2);

      Sp = predVar.S;

      if strcmpi(kcriterion,'MLE')
        % kappa with the MAXIMUM LIKELIHOOD is selected           
        psi = zeros(1,length(kappa));
        for i=1:length(kappa)        
          % predictive sigma-points
          dzetap = zeros(nz,Jns);
          if h.isAdditive
            [chip,w] = nefSUKF.smsp(predMean,Sp,kappa(i));
            % points for computation prediction of measurement
            dzetap = evaluateOverAllStates(h,chip,Input,vMean,Time);
          else
            [chip,w] = nefSUKF.smsp([predMean; vMean],[Sp, zeros(nx,nv); zeros(nv,nx), Sr],kappa(i));
            % points for computation prediction of measurement
            dzetap = evaluate(h,chip(1:nx,:),Input(:,ones(1,Jns)),chip(nx+1:nx+nv,:),Time(:,ones(1,Jns)));
          end


          % measurement prediction mean
          zPredMean = dzetap*w'; % zp          

          % measurement prediction variance
          dzetap_ = bsxfun(@minus,dzetap,zPredMean);
          if h.isAdditive
            dzetap_sw = dzetap_.*sqrt(w(ones(1,nz),:));
            Szp = [dzetap_sw, Sr];% Szp; Pzzp = Szp*Szp'
            zPredVar = Szp*Szp'; % Pzzp
          else
            dzetap_sw = dzetap_.*sqrt(w(ones(1,nz),:)); % Szp; Pzzp = Szp*Szp'
            zPredVar = dzetap_sw*dzetap_sw'; % Pzzp
          end

          %psi(i) = mvnpdf(Measurement,zPredMean,zPredVar);                
          exponent = -0.5*((Measurement-zPredMean)'/zPredVar)*(Measurement-zPredMean);
          psi(i)  = 1/(2*pi)^(nz/2)/sqrt(det(zPredVar))*exp(exponent);
        end      
        kappa_act = kappa(logical(psi==max(psi)));
      elseif strcmpi(kcriterion,'MMPE')
        % kappa providing the MINIMUM MEASUREMENT PREDICTION ERROR is selected
        mpe = zeros(1,length(kappa));      
        for i=1:length(kappa)
          % predictive sigma-points
          dzetap = zeros(nz,Jns);
          if h.isAdditive
            [chip,w] = nefSUKF.smsp(predMean,Sp,kappa(i));
            % points for computation prediction of measurement
            dzetap = evaluateOverAllStates(h,chip,Input,vMean,Time);
          else
            [chip,w] = nefSUKF.smsp([predMean; vMean],[Sp, zeros(nx,nv); zeros(nv,nx), Sr],kappa(i));
            % points for computation prediction of measurement
            dzetap = evaluate(h,chip(1:nx,:),Input(:,ones(1,Jns)),chip(nx+1:nx+nv,:),Time(:,ones(1,Jns)));
          end

          % measurement prediction mean
          zPredMean = dzetap*w'; % zp        

          % measurement prediction error
          mpe(i)  = (Measurement-zPredMean)'*(Measurement-zPredMean);
        end % for  
        kappa_act = kappa(logical(mpe==min(mpe)));
      end % if
      % for the case of several values of kappa with same likelihood
      kappa_act = kappa_act(1);
    end %function scaling parameter adaptive choice


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UKFTimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = SUKFTimeUpdate(filtMean,filtSVar,f,Input,wMean,Sq,Time,kappa)
      % SUKFTIMEUPDATE Implements time update of the SUKF.

      %numval = size(filtMean,2);

      %  number of sigma-points
      nx = size(filtMean,1);
      nw = size(wMean,1);
      if f.isAdditive
        Jns = 2*nx+ 1;
      else
        Jns = 2*(nx+nw)+ 1;
      end


      % filtering sigma-points
      chip = zeros(nx,Jns);
      if f.isAdditive
        [chi,w] = nefSUKF.smsp(filtMean,filtSVar,kappa);
        % computing predictive sigma-points
        chip = evaluateOverAllStates(f,chi,Input,wMean,Time);
      else
        [chi,w] = nefSUKF.smsp([filtMean; wMean],[filtSVar, zeros(nx,nw); zeros(nw,nx), Sq],kappa);
        % computing predictive sigma-points
        chip = evaluate(f,chi(1:nx,:),Input(:,ones(1,Jns)),chi(nx+1:nx+nw,:),Time(:,ones(1,Jns)));
      end

      % evaluating prediction mean
      predMean = chip*w';

      % evaluating prediction variance
      chip_ = bsxfun(@minus,chip,predMean);
      chip_sw = chip_.*sqrt(w(ones(1,nx),:));

      if f.isAdditive
        Sp = nefDD1.triag([chip_sw, Sq]);
      else
        Sp = nefDD1.triag(chip_sw);
      end
      predVar = nefCholFactorFunction(Sp,'matrixType','fact');
      %end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UKFMeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = SUKFMeasurementUpdate(predMean,predSVar,h,Input,vMean,Sr,Time,Measurement,kappa)
      % SUKFMEASUREMENTUPDATE Implements measurement update of the SUKF. 

      %      numval = size(predMean,2);

      % number of sigma points
      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);
      if h.isAdditive
        Jns = 2*(nx) + 1;
      else
        Jns = 2*(nx+nv) + 1;
      end

      %       filtMean = zeros(nx,numval);
      %       zPredMean = zeros(nz,numval);
      %       invzPredVar = zeros(nz,nz,numval);
      %       detzPredVar = zeros(1,numval);
      %       filtVar = zeros(nx,nx,numval);

      %for in = 1:numval
      Sp = predSVar;
      % predictive sigma-points
      dzetap = zeros(nz,Jns);
      if h.isAdditive
        [chip,w] = nefSUKF.smsp(predMean,Sp,kappa);
        % points for computation prediction of measurement
        dzetap = evaluateOverAllStates(h,chip,Input,vMean,Time);
      else
        [chip,w] = nefSUKF.smsp([predMean; vMean],[Sp, zeros(nx,nv); zeros(nv,nx), Sr],kappa);
        % points for computation prediction of measurement
        dzetap = evaluate(h,chip(1:nx,:),Input(:,ones(1,Jns)),chip(nx+1:nx+nv,:),Time(:,ones(1,Jns)));
      end

      % measurement prediction mean and prediction error
      zPredMean = dzetap*w'; % zp
      dz = Measurement-zPredMean;

      % measurement prediction variance
      dzetap_ = bsxfun(@minus,dzetap,zPredMean);
      if h.isAdditive
        dzetap_sw = dzetap_.*sqrt(w(ones(1,nz),:));
        Szp = [dzetap_sw, Sr];% Szp; Pzzp = Szp*Szp'
      else
        dzetap_sw = dzetap_.*sqrt(w(ones(1,nz),:)); % Szp; Pzzp = Szp*Szp'
      end

      % cross-covariance predictive matrix of x and z
      chip_ = bsxfun(@minus,chip(1:nx,:),predMean);
      chip_sw = chip_.*sqrt(w(ones(1,nx),:));
      Pxzp = chip_sw*dzetap_sw';

      % evaluating filtering mean and variance
      if h.isAdditive
        K = Pxzp/Szp'/Szp;
        Sf = nefDD1.triag([chip_sw-K*dzetap_sw, K*Sr]);
      else
        K = Pxzp/dzetap_sw'/dzetap_sw;
        Sf = nefDD1.triag(chip_sw-K*dzetap_sw);
      end
      filtMean = predMean + K*dz;

      filtVar = nefCholFactorFunction(Sf,'matrixType','fact');
      %end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        if h.isAdditive
          zPredVar = Szp*Szp'; % Pzzp
        else
          zPredVar = dzetap_sw*dzetap_sw'; % Pzzp
        end
        % used for GSM filter
        % not for a particle filter - TODO
        varargout{3} = zPredMean;
        varargout{4} = inv(zPredVar);
        varargout{5} = det(zPredVar);
      end
    end %function measurement update

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % URTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = SURTSSmoothUpdate(initMean,initVarfac,filtMean,filtVarfac,predMean,predVarfac,wMean,Sq,f,Time,Input,kappa)
      % SURTSSMOTHUPDATE Implements smoothing step based on the SUKF and Rauch-Tung-Striebel smoother.

      nx = size(filtMean,1);
      nw = size(wMean,1);
      if additive
        Jns = 2*(nx) + 1;
      else
        Jns = 2*(nx+nw) + 1;
      end

      Sf = filtVarfac;

      % filtering sigma-points
      chip_v = zeros(nx, Jns);
      if f.isAdditive
        [chi_v,w] = nefSUKF.smsp(filtMean,Sf,kappa);
        % computing predictive sigma-points
        chip_v = evaluateOverAllStates(f,chi_v,Input,wMean,Time);
      else
        [chi_v,w] = nefSUKF.smsp([filtMean; wMean],[Sf, zeros(nx,nw); zeros(nw,nx), Sq],kappa);
        % computing predictive sigma-points
        chip_v = evaluate(f,chi_v(1:nx,:),Input(:,ones(1,Jns)),chi_v(nx+1:nx+nw,:),Time(:,ones(1,Jns)));
      end


      % predictive statistics for smoothing - mean, covariance matrix
      % mean
      xpvv = chip_v*w';

      % covariance matrix Pxx = E{(x(k+1)-x_pred(k+1))(x(k)-x_filt(k))}
      if f.isAdditive
        chi_v_ = bsxfun(@minus,chi_v,filtMean);
      else
        chi_v_ = bsxfun(@minus,chi_v,[filtMean; wMean]);
      end
      chip_v_ = bsxfun(@minus,chip_v,xpvv);
      if f.isAdditive
        chi_v_sw = chi_v_.*sqrt(w(ones(1,nx),:));
        chip_v_sw = chip_v_.*sqrt(w(ones(1,nx),:));
      else
        chi_v_sw = chi_v_.*sqrt(w(ones(1,nx+nw),:));
        chip_v_sw = chip_v_.*sqrt(w(ones(1,nx),:));
      end

      % Evaluation of the covariance cov(x_predictive, x_filrering),
      % i.e. E{(x(k+1)-x_pred(k+1))(x(k)-x_filt(k))}
      Pxx = chi_v_sw(1:nx,:)*chip_v_sw';


      % it should hold chip_v_*chip_v_w' =  predVar

      % smoother gain
      Kv = Pxx/predVarfac'/predVarfac;

      % smoothed mean and covariance matrix
      smoothMean = filtMean - Kv * (predMean - initMean);
      if additive
        Sv = nefDD1.triag([chi_v_sw(1:nx,:)-Kv*chip_v_sw, Kv*Sq, Kv*initVarfac]);
      else
        Sv = nefDD1.triag([chi_v_sw(1:nx,:)-Kv*chip_v_sw, Kv*initVarfac]);
      end
      smoothVar = nefCholFactorFunction(Sv,'matrixType','fact');
    end

  end % methods
end % the class
