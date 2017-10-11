classdef nefSDD1 < nefEstimator
  %file @nefSDD1/nefSDD1.m (Square-Root Divided Difference Estimator 1st Order)
  % nefSDD1 Properties:
  %   nefDD1_optsDef - definition of the nefSDD1 options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefSDD1 Methods:
  %   nefSDD1 - class constructor (For list of options type 'help nefSDD1/nefSDD1')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   SDD1TimeUpdate - time update of the SDD1
  %   SDD1MeasurementUpdate - measurement update of the SDD1
  %   SD1RTSSmoothUpdate - smoothing based on the SDD1 and Rauch-Tung-Striebel smoother
  %
  % References:
  % M. Norgaard, N. K. Poulsen and O. Ravn (2000):
  % New developments in state estimation for nonlinear systems.
  % Automatica 36(11), 1627-1638.
  %
  % M. Simandl and J. Dunik (2009):
  % Derivative-Free Estimamation Methods: New Results and Performance Analysis
  % Automatica 45(7), 1749-1757.
  %
  % M. Simandl and J. Dunik (2006):
  % Design of derivative-free smoothers and predictors.
  % In: Preprints of the 14th IFAC Symposium on System Identification.
  % Newcastle, Australia.
  %
  % J. Dunik, M. Simandl and O. Straka (2010): 
  % Adaptive choice of scaling parameter in derivative-free local filters.
  % In: Proceedings of the 13th International Conference on Information
  % Fusion, Edinburgh, UK.
  
  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia
  
  properties (Constant)
    % definition of estimator options    
    nefSDD1_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
                      2,'scalingParameterType:constant','parameterValue',sqrt(3),'constant scaling parameter', @(x) isfloat(x) && length(x)==1 && x>0;
                      2,'scalingParameterType:adaptive','parameterDomain',[0.1 2.1 0.5],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && x(1)>0 && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
                      2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE'}));};
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
    %hpar;
    isRConst;
    isQConst;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSDD1(system,varargin)
      % NEFSDD1 Creates NEFSDD1 object.
      %
      %   OBJ = NEFSDD1(system,varargin) creates a NEFSDD1 object OBJ 
      %   representing square-root divided difference filter 1st order for model SYSTEM.
      %
      %   the following user-defined parameters can be changed 
      %   via the standard Parameter-Value MATLAB mechanism
      %  
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  scalingParameterType        constant         scaling parameter can be             constant, adaptive
      %                                               constant or adaptively chosen   
      %  parameterValue              sqrt(3)          constant scaling parameter           real, positive
      %  parameterDomain             [0.1 2.1 0.5]    adaptive grid search of scaling      positive real, Kmin<Kmax
      %                                               parameter K within domain 
      %                                               [Kmin, Kmax, Kstep] for 
      %                                               adaptive setting
      %  parameterCriterion          MLE              selection of criterion for           MLE                       
      %                                               adaptive search
      %
      %  See also ESTIMATE. 

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFSDD1';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for i = 1:size(nefSDD1.nefSDD1_optsDef,1)
        p.addParamValue(nefSDD1.nefSDD1_optsDef{i,3},[],nefSDD1.nefSDD1_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefSDD1.nefSDD1_optsDef,p.Results);      
      
      % scaling parameter check
      obj = nefDD1.ScalingParameterCheck(obj);

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

      obj.Sr = [];
      obj.Sq = [];
      obj.wSVar = [];
      obj.vSVar = [];

      obj.isRConst =  obj.sys.v.VarNumeric;
      obj.isQConst =  obj.sys.v.VarNumeric;
    end % NEFSDD1 constructor
    
    function disp(obj)
      % DISP Displays nefSDD1 object parameters and options.
      
      fprintf('A nefSDD1 object with parameters\n')
      showOptions(obj)
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time,varargin)
    % TIMEUPDATE Performs time update step with respect to system
    % description.

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
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen - value from the last filtering step
          aux_domain = getOption(obj,'parameterDomain');
          hpar = aux_domain(4);
        end
        [predMean, predVar] = nefSDD1.SDD1TimeUpdate(filtMean,filtVar.S,obj.sys.f,Input,wMean,obj.wSVar,Time,hpar);
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
          predEstimate.auxData.hpar = hpar;
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
    % MEASUREMENTUPDATE Performs measurement update step with respect to
    % system description.

      % parameters of prediction estimate
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = predEst.RV.Var;
      %predVar = evalVariance(predEst.RV,[],[],[],[]);      
      % mean and variancec of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);

      % obtain S factor of R = vVar
      % for constant R,  the factor is computed at the first step only
      % for t-variant R,  the factor is computed at each time step
      if ((obj.isRConst == 1) && ((isempty(obj.vSVar)))) || (obj.isRConst == 0)
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
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen
          aux_domain = getOption(obj,'parameterDomain');
          aux_criterion = getOption(obj,'parameterCriterion');
          hpar = nefSDD1.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,aux_domain(1:3),aux_criterion);                         
          setOption(obj,'parameterDomain',[aux_domain(1:3), hpar]);
        end
        [filtMean,filtVar] = nefSDD1.SDD1MeasurementUpdate(predMean,predVar.S,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,hpar);
      end

      % returning nefGaussianRV with symmetrized variance
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
      % return also by-products if required
    end % function measurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SMOOTHUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothEstimate] = smoothUpdate(obj,initEstimate,filtEstimate,predEstimate)
    % SMOOTHUPDATE Performs smoothing step with respect to system description.  
    
      initM = evalMean(initEstimate.RV);
      initVfac = initEstimate.RV.Var.S;
      filtM = evalMean(filtEstimate.RV);
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
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen - stored from the prediction
          hpar = predEstimate.auxData.hpar;
        end
        [smoothM,smoothV] = nefSDD1.SD1RTSSmoothUpdate(initM,initVfac,filtM,filtVfac,predM,predVfac,wMean,obj.Sq,obj.sys.f,Time,Input,hpar);
      end

      % returning nefGaussianRV with symmetrized variance
      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate


  end % methods
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ScalingParameterCheck
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % included in nefDD1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MatrixTriangularization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % included in nefDD1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scaling Parameter Adaptive Choice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hpar_act = ScalingParameterAdaptiveChoice(predMean,predVar,h,Input,vMean,Sr,Time,Measurement,kdomain,kcriterion)
    % SCALINGPARAMETERADAPTIVECHOICE Performs adaptive choice of the
    % scaling parameter within predefined domain in filtering step.
      
      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);
      
      % grid points
      hpar = kdomain(1):kdomain(3):kdomain(2);
      
      Sp = predVar.S;      
      
      if strcmpi(kcriterion,'MLE')
        % hpar with the MAXIMUM LIKELIHOOD is selected     
        psi = zeros(1,length(hpar));
        for i=1:length(hpar)
          % Auxiliary matrices Szx1, Szv1 containing divided
          % differences of first order
          xd = evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar(i),predMean),Input,vMean,Time);
          xd_ = evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar(i),predMean),Input,vMean,Time);
          Szx1 = (xd - xd_)/(2*hpar(i));
          if h.isAdditive
            Szv1 = Sr;
          else
            xd = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            xd_ = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            Szv1 = (xd - xd_)/(2*hpar(i));
          end
          % measurement prediction mean and prediction error
          zPredMean = evaluate(h,predMean,Input,vMean,Time); % zp
          % measurement prediction variance        
          zPredVar = [Szx1 Szv1]*[Szx1 Szv1]'; % Pzzp

          %psi(i) = mvnpdf(Measurement,zPredMean,zPredVar);        
          exponent = -0.5*((Measurement-zPredMean)'/zPredVar)*(Measurement-zPredMean);
          psi(i)  = 1/(2*pi)^(nz/2)/sqrt(det(zPredVar))*exp(exponent);
        end % for
        hpar_act = hpar(logical(psi==max(psi)));
      %elseif strcmpi(kcriterion,.....)
      end % if
        
      % for the case of several values of hpar with same likelihood
      hpar_act = hpar_act(1);
    end %function scaling parameter adaptive choice

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SDD1TimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = SDD1TimeUpdate(filtMean,filtSVar,f,Input,wMean,Sq,Time,hpar,additive)
    % SDD1TIMEUPDATE Implements time update of the SDD1.
      
      %numval = size(filtMean,2);

      nx = size(filtMean,1);
      nw = size(wMean,1);

      %for in = 1:numval
      Sf = filtSVar;

      % Auxiliary matrices Sxx1, Sxw1  containing divided
      % differences of first order
      xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean),Input,wMean,Time);
      xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean),Input,wMean,Time);
      Sxx1 = (xd - xd_)/(2*hpar);
      if f.isAdditive
        Sxw1 = Sq;
      else
        xd = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,Sq*hpar,wMean),Time(:,ones(1,nw)));
        xd_ = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,-Sq*hpar,wMean),Time(:,ones(1,nw)));
        Sxw1 = (xd - xd_)/(2*hpar);
      end;

      % evaluating prediction mean
      predMean = evaluate(f,filtMean,Input,wMean,Time);

      % evaluating prediction variance
      Sp = nefDD1.triag([Sxx1 Sxw1]);
      predVar = nefCholFactorFunction(Sp,'matrixType','fact');
      %end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SDD1MeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = SDD1MeasurementUpdate(predMean,predSVar,h,Input,vMean,Sr,Time,Measurement,hpar)
    % SDD1MEASUREMENTUPDATE Implements measurement update of the SDD1.
      
      %numval = size(predMean,2);

      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);

      %for in = 1:numval
      Sp = predSVar;

      % Auxiliary matrices Szx1, Szv1 containing divided
      % differences of first order
      xd = evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar,predMean),Input,vMean,Time);
      xd_ = evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar,predMean),Input,vMean,Time);
      Szx1 = (xd - xd_)/(2*hpar);
      if h.isAdditive
        Szv1 = Sr;
      else
        xd = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar,vMean),Time(:,ones(1,nv)));
        xd_ = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar,vMean),Time(:,ones(1,nv)));
        Szv1 = (xd - xd_)/(2*hpar);
      end

      % measurement prediction mean and prediction error
      zPredMean = evaluate(h,predMean,Input,vMean,Time); % zp
      dz = Measurement-zPredMean;

      % measurement prediction variance
      Szp = [Szx1 Szv1];

      % cross-covariance predictive matrix of x and z
      Pxzp = Sp*Szx1';

      % evaluating filtering mean and variance
      K = Pxzp/Szp'/Szp;
      %K = (Pxzp/nefDD1.triag(Szp)')/nefDD1.triag(Szp)';
      filtMean = predMean + K*dz;
      Sf = nefDD1.triag([Sp-K*Szx1, K*Szv1]);
      filtVar = nefCholFactorFunction(Sf,'matrixType','fact');
      %end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        % used for GSM filter, not for a particle filter TODO
        zPredVar = Szp*Szp'; % Pzzp
        varargout{3} = zPredMean;
        varargout{4} = inv(zPredVar);
        varargout{5} = det(zPredVar);
        %         invzPredVar = inv(zPredVar);
        %         detzPredVar = det(zPredVar);
        %         varargout{3} = zPredMean;
        %         varargout{4} = invzPredVar;
        %         varargout{5} = detzPredVar;
      end
    end %function measurement update

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SD1RTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = SD1RTSSmoothUpdate(initMean,initVarfac,filtMean,filtVarfac,predMean,predVarfac,wMean,Sq,f,Time,Input,hpar,additive)
    % SD1RTSSMOOTHUPDATE Implements smoothing step based on the SDD1 and
    % Rauch-Tung-Striebel smoother.
      
      
      nx = size(filtMean,1);
      nw = size(wMean,1);

      Sf = filtVarfac;
      initSv = initVarfac;

      %  Auxiliary matrix Sxx1 containing divided differences
      xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean),Input,wMean,Time);
      xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean),Input,wMean,Time);
      Sxx1 = (xd - xd_)/(2*hpar);
      if f.isAdditive
        Sxw1 = Sq;
      else
        sj = Sq(:,idx);
        xd = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,Sq*hpar,wMean),Time(:,ones(1,nw)));
        xd_ = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,-Sq*hpar,wMean),Time(:,ones(1,nw)));
        Sxw1 = (xd - xd_)/(2*hpar);
      end;

      % covariance matrix Pxx = E{(x(k+1)-x_pred(k+1))(x(k)-x_filt(k))}
      Pxx = Sf*Sxx1';

      % it should hold Sxx1*Sxx1' + wVar =   predVar

      % smoother gain
      %Kv = Pxx * inv(predVarfac*predVarfac');
      Kv = (Pxx/predVarfac')/predVarfac;

      % smoothed mean and covariance matrix
      smoothMean = filtMean - Kv * (predMean - initMean);
      %Sv = nefDD1.triag([Sf-Kv*Sxx1, Kv*Sq, Kv*initSv]);
      Sv = nefDD1.triag([Sf-Kv*Sxx1, Kv*Sxw1, Kv*initSv]);
      smoothVar = nefCholFactorFunction(Sv,'matrixType','fact');
    end %SmoothUpdate

  end % methods
end % the class
