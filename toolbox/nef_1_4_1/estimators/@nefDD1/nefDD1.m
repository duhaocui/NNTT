classdef nefDD1 < nefEstimator
  %file @nefDD1/nefDD1.m (Divided Difference Estimator 1st Order)
  % nefDD1 Properties:
  %   nefDD1_optsDef - definition of the nefDD1 options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefDD1 Methods:
  %   nefDD1 - class constructor (For list of options type 'help nefDD1/nefDD1')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterCheck - check of chosen scaling parameter and warning spec.
  %   triag - matrix trinagularization via QR
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   DD1TimeUpdate - time update of the DD1
  %   DD1MeasurementUpdate - measurement update of the DD1
  %   D1RTSSmoothUpdate - smoothing based on the DD1 and Rauch-Tung-Striebel smoother
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
    nefDD1_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
    2,'scalingParameterType:constant','parameterValue',sqrt(3),'constant scaling parameter', @(x) isfloat(x) && length(x)==1 && x>0;
    2,'scalingParameterType:adaptive','parameterDomain',[0.1 2.1 0.5],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && x(1)>0 && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
    2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE'}));};
  end % constant properties


  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefDD1(system,varargin)      
      % NEFDD1 Creates NEFDD1 object.
      %
      %   OBJ = NEFDD1(system,varargin) creates a NEFDD1 object OBJ 
      %   representing divided difference filter 1st order for model SYSTEM.
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
      %
      %  See also ESTIMATE.      

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFDD1';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for i = 1:size(nefDD1.nefDD1_optsDef,1)
        p.addParamValue(nefDD1.nefDD1_optsDef{i,3},[],nefDD1.nefDD1_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefDD1.nefDD1_optsDef,p.Results);      

      % scaling parameter check
      obj = nefDD1.ScalingParameterCheck(obj);

      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end
      obj.estQueue(end+1) = data;

    end % NEFDD1 constructor

    function disp(obj)
      % DISP Displays nefDD1 object parameters and options.

      fprintf('A nefDD1 object with parameters\n')
      showOptions(obj)
    end % function disp

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
      filtVar = evalVariance(filtEst.RV,[],[],[],[]);
      % mean and variance of state noise
      wMean = evalMean(obj.sys.w,[],Input,[],Time);
      wVar = evalVariance(obj.sys.w,[],Input,[],Time);

      % evaluating predictive mean and variance
      % - if state eq. is linear (with additive noise), then the predictive part of the KF is used
      if obj.sys.f.isLinear
        % this part of the KF is not in factorized (numerically more stable) form
        % but it is slightly faster
        F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
        Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
        [predMean,predVar] = nefKalman.KalmanTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,F,Gamma);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen - value from the last filtering step
          aux_domain = getOption(obj,'parameterDomain');
          hpar = aux_domain(4);
        end
        [predMean, predVar] = nefDD1.DD1TimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,hpar);
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
      predVar = evalVariance(predEst.RV,[],[],[],[]);
      % mean and variancec of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);

      % evaluating filtering mean and variance
      % - if measurement eq. is linear (with additive noise), then the predictive part of the KF is used
      if obj.sys.h.isLinear
        % this part of the KF is not in factorized (numerically more stable) form
        % but it is slightly faster
        H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
        Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
        [filtMean,filtVar] = nefKalman.KalmanMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,H,Delta);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen
          aux_domain = getOption(obj,'parameterDomain');
          aux_criterion = getOption(obj,'parameterCriterion');
          hpar = nefDD1.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,aux_domain(1:3),aux_criterion);                         
          setOption(obj,'parameterDomain',[aux_domain(1:3), hpar]);
        end
        [filtMean,filtVar] = nefDD1.DD1MeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,hpar);
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
      initV = evalVariance(initEstimate.RV);
      filtM = evalMean(filtEstimate.RV);
      filtV = evalVariance(filtEstimate.RV);
      predM = evalMean(predEstimate.RV);
      predV = evalVariance(predEstimate.RV);
      Input = predEstimate.auxData.Input;
      Time = predEstimate.auxData.Time;
      wMean = predEstimate.auxData.wMean;
      wVar = predEstimate.auxData.wVar;

      % evaluating smoothing mean and variance
      % - if state eq. is linear (with additive noise), then the predictive part of the KF is used
      if obj.sys.f.isLinear
        % this part of the KF is not in factorized (numerically more stable) form
        % but it is slightly faster        
        F = predEstimate.auxData.F;
        [smoothM,smoothV] = nefKalman.RTSKalmanSmoothUpdate(initM,initV,filtM,filtV,predM,predV,F);
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if hpar is constant during the experiment
          hpar = getOption(obj,'parameterValue');
        else
          % hpar is adaptively chosen - stored from the prediction
          hpar = predEstimate.auxData.hpar;
        end
        [smoothM,smoothV] = nefDD1.D1RTSSmoothUpdate(initM,initV,filtM,filtV,predM,predV,wMean,wVar,obj.sys.f,Time,Input,hpar);
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
        % hpar cannot be zero or negative      
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
          disp('nef(S)DD1 WARNING: Function in measurement equation is linear and independent of the scaling parameter. Its optimisation is useless.')
          %disp('Better choice of the scaling parameter is "constant" (instead of adaptive).')
          disp('**************')
        end                
      end
    end % function ScalingParameterCheck

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MatrixTriangularization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function M = triag(N)
      % TRIAG Performs matrix triangularization via QR decomposition.

      [~,M]=qr(N',0);
      M = M';
    end% function triag

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scaling Parameter Adaptive Choice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hpar_act = ScalingParameterAdaptiveChoice(predMean,predVar,h,Input,vMean,R,Time,Measurement,kdomain,kcriterion)
      % SCALINGPARAMETERADAPTIVECHOICE Performs adaptive choice of the
      % scaling parameter within predefined domain in filtering step.

      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);

      % grid points
      hpar = kdomain(1):kdomain(3):kdomain(2);

      Sp = nefCholFactorFunction.CholFactor(predVar);
      Sr = nefCholFactorFunction.CholFactor(R);

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
        end   % for   
        hpar_act = hpar(logical(psi==max(psi)));
        %elseif strcmpi(kcriterion,.....)
      end % if


      % for the case of several values of hpar with same likelihood
      hpar_act = hpar_act(1);
    end %function scaling parameter adaptive choice

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DD1TimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = DD1TimeUpdate(filtMean,filtVar,f,Input,wMean,Q,Time,hpar)
      % DD1TIMEUPDATE Implements time update of the DD1.

      numval = size(filtMean,2);

      nx = size(filtMean,1);
      nw = size(wMean,1);

      predMean = zeros(nx,numval);
      predVar = zeros(nx,nx,numval);

      for in = 1:numval
        Sf = nefCholFactorFunction.CholFactor(filtVar(:,:,in));
        Sq = nefCholFactorFunction.CholFactor(Q(:,:,in));

        % Auxiliary matrices Sxx1, Sxw1  containing divided
        % differences of first order
        xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean(:,in)),Input,wMean(:,in),Time);
        xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean(:,in)),Input,wMean(:,in),Time);
        Sxx1 = (xd - xd_)/(2*hpar);
        if f.isAdditive
          Sxw1 = Sq;
        else
          xd = evaluate(f,filtMean(:,in*ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,Sq*hpar,wMean(:,in)),Time(:,ones(1,nw)));
          xd_ = evaluate(f,filtMean(:,in*ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,-Sq*hpar,wMean(:,in)),Time(:,ones(1,nw)));
          Sxw1 = (xd - xd_)/(2*hpar);
        end;

        % evaluating prediction mean
        predMean(:,in) = evaluate(f,filtMean(:,in),Input,wMean(:,in),Time);

        % evaluating prediction variance
        %Sp = nefDD1.triag([Sxx1 Sxw1]);
        Sp = [Sxx1 Sxw1];
        predVar(:,:,in) = Sp*Sp';
      end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DD1MeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = DD1MeasurementUpdate(predMean,predVar,h,Input,vMean,R,Time,Measurement,hpar)
      % DD1MEASUREMENTUPDATE Implements measurement update of the DD1.

      numval = size(predMean,2);

      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);

      filtMean = zeros(nx,numval);
      zPredMean = zeros(nz,numval);
      zPredVar = zeros(nz,nz,numval);
      filtVar = zeros(nx,nx,numval);

      for in = 1:numval
        Sp = nefCholFactorFunction.CholFactor(predVar(:,:,in));
        Sr = nefCholFactorFunction.CholFactor(R(:,:,in));

        % Auxiliary matrices Szx1, Szv1 containing divided
        % differences of first order
        xd = evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar,predMean(:,in)),Input,vMean(:,in),Time);
        xd_ = evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar,predMean(:,in)),Input,vMean(:,in),Time);
        Szx1 = (xd - xd_)/(2*hpar);
        if h.isAdditive
          Szv1 = Sr;
        else
          xd = evaluate(h,predMean(:,in*ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar,vMean(:,in)),Time(:,ones(1,nv)));
          xd_ = evaluate(h,predMean(:,in*ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar,vMean(:,in)),Time(:,ones(1,nv)));
          Szv1 = (xd - xd_)/(2*hpar);
        end

        % measurement prediction mean and prediction error
        zPredMean(:,in) = evaluate(h,predMean(:,in),Input,vMean(:,in),Time); % zp
        dz = Measurement(:,in)-zPredMean(:,in);

        % measurement prediction variance
        Szp = [Szx1 Szv1];
        zPredVar(:,:,in) = Szp*Szp'; % Pzzp

        % cross-covariance predictive matrix of x and z
        Pxzp = Sp*Szx1';

        % evaluating filtering mean and variance
        K = Pxzp/zPredVar(:,:,in);
        filtMean(:,in) = predMean(:,in) + K*dz;
        %Sf = nefDD1.triag([Sp-K*Szx1, K*Szv1]);
        Sf = [Sp-K*Szx1, K*Szv1];
        filtVar(:,:,in) = Sf*Sf';
      end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        % used for a particle and GSM filter
        invzPredVar = zeros(nz,nz,numval);
        detzPredVar = zeros(1,numval);
        for in = 1:numval
          invzPredVar(:,:,in) = inv(zPredVar);
          detzPredVar(:,in) = det(zPredVar);
        end
        varargout{3} = zPredMean;
        varargout{4} = invzPredVar;
        varargout{5} = detzPredVar;
      end
    end %function measurement update

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % D1RTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = D1RTSSmoothUpdate(initMean,initVar,filtMean,filtVar,predMean,predVar,wMean,wVar,f,Time,Input,hpar) %#ok<INUSD>
      % D1RTSSMOOTHUPDATE Implements smoothing step based on the DD1 and Rauch-Tung-Striebel smoother.

      nx = size(filtMean,1);
      Sf = nefCholFactorFunction.CholFactor(filtVar);

      %  Auxiliary matrix Sxx1 containing divided differences
      xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean),Input,wMean,Time);
      xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean),Input,wMean,Time);
      Sxx1 = (xd - xd_)/(2*hpar);

      % covariance matrix Pxx = E{(x(k+1)-x_pred(k+1))(x(k)-x_filt(k))}
      Pxx = Sf*Sxx1';

      % it should hold Sxx1*Sxx1' +wVar =   predVar

      % smoother gain
      Kv = Pxx/predVar;

      % smoothed mean and covariance matrix
      smoothMean = filtMean - Kv * (predMean - initMean);
      smoothVar = filtVar - Kv * (predVar - initVar) * Kv';
      smoothVar = (smoothVar+smoothVar')/2;
    end %SmoothUpdate

  end % methods
end % the class
