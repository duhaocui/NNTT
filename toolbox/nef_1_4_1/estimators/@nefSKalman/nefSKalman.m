classdef nefSKalman < nefEstimator
  %file @nefSKalman/nefSKalman.m (Square-Root (Extended) Kalman Filter and Rauch-Tung-Striebel smoother - based on triangularization)
  % nefSKalman Methods:
  %   nefSKalman - class constructor
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   SKalmanTimeUpdate - time update of the square-root (E)KF
  %   SKalmanMeasurementUpdate - measurement update of the square-root (E)KF
  %   SKalmanRTSSmoothUpdate - smoothing based on the square-root (extended) Rauch-Tung-Striebel smoother
  %
  % Similar algorithms can be found e.g. in the following references
  %
  % D. Simon (2006):
  % Optimal State Estimation: Kalman, H infinity and Nonlinear Approaches.
  % John Wiley & Sons, Inc.
  %
  % M. S. Grewal, A. P. Andrews (2001):
  % Kalman Filtering: Theory and Practise Using MATLAB.
  % John Wiley & Sons, Inc.

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%
    % MATRICES
    %%%%%%%%%%%%%%%%%%%%%
    Sr = [];
    Sq = [];
    vSVar = [];
    wSVar = [];
    %%%%%%%%%%%%%%%%%%%%%%%
    %PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    isRConst;
    isQConst;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSKalman(system,varargin)
      % NEFSKALMAN Creates NEFSKALMAN object.
      %
      %   OBJ = NEFSKALMAN(system,varargin) creates a NEFSKALMAN object OBJ 
      %   representing square-root Kalman filter or Rauch-Tung-Striebel smoother for model
      %   SYSTEM.
      
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFSKALMAN';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

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
      obj.vSVar = [];
      obj.wSVar = [];

      obj.isRConst =  isVarNumeric(obj.sys.v);
      obj.isQConst =  isVarNumeric(obj.sys.w);
    end % NEFSKALMAN constructor
    
    function disp(obj)
      % DISP Displays nefSKalman object parameters and options.
      
      fprintf('A nefSKalman object with parameters\n')
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
      % mean and variance of state noise
      wMean = evalMean(obj.sys.w,[],Input,[],Time);
      wVar = evalVariance(obj.sys.w,[],Input,[],Time);
      % matrices F for state and Gamma for noise in state equation
      F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
      Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);

      % obtain S factor of Q = wVar
      % for constant Q,  the factor is computed at the first step only
      % for variable Q,  the factor is computed at each time step
      if ((obj.isQConst == 1) && (isempty(obj.wSVar))) || (obj.isQConst == 0)
        obj.wSVar = nefCholFactorFunction.CholFactor(wVar);
      end

      % evaluating predictive mean and variance
      [predMean, predVar] = nefSKalman.SKalmanTimeUpdate(filtMean,filtVar.S,obj.sys.f,Input,wMean,obj.wSVar,Time,F,Gamma);

      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);

      % is timUpdate is used for smoothing
      if p.Results.smoothingPurpose
        predEstimate.auxData.Input = Input;
        predEstimate.auxData.Time = Time;
        predEstimate.auxData.wMean = wMean;
        predEstimate.auxData.wVar = wVar;
        predEstimate.auxData.F = F;
        predEstimate.auxData.Gamma = Gamma;
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
      % mean and variancec of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);
      % matrices H for state and Delta for noise in state equation
      H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
      Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);

      % obtain S factor of R = vVar
      % for constant R,  the factor is computed at the first step only
      % for variable R,  the factor is computed at each time step
      if ((obj.isRConst == 1) && (isempty(obj.vSVar))) || (obj.isRConst == 0)
        obj.vSVar = nefCholFactorFunction.CholFactor(vVar);
      end

      % evaluating filtering mean and variance
      [filtMean,filtVar] = nefSKalman.SKalmanMeasurementUpdate(predMean,predVar.S,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,H,Delta);

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
      F = predEstimate.auxData.F;
      Gamma = predEstimate.auxData.Gamma;

      % obtain S factor of Q = wVar
      % for constant Q,  the factor is computed at the first step only
      % for variable Q,  the factor is computed at each time step
      %if ((obj.isQConst == 1) && (Time == 0)) || (obj.isQConst == 0)
      factor = nefCholFactorFunction(wVar,'matrixType','cov');
      obj.Sq = factor.S;
      %end

      % evaluating smoothed mean and variance
      [smoothM,smoothV] = nefSKalman.SKalmanRTSSmoothUpdate(initM,initVfac,filtM,filtVfac,predM,predVfac,F,Gamma,wMean,obj.Sq,obj.sys.f,Time,Input);

      % returning nefGaussianRV with symmetrized variance
      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate


  end % methods
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SKalmanTimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = SKalmanTimeUpdate(filtMean,filtSVar,f,Input,wMean,wSVar,Time,F,Gamma)
    % SKALMANTIMEUPDATE Implements time update of the square-root (E)KF.
      
      % evaluating prediction mean
      predMean = evaluate(f,filtMean,Input,wMean,Time);

      % evaluating prediction variance
      Sp = nefDD1.triag([F*filtSVar Gamma*wSVar]);
      predVar = nefCholFactorFunction(Sp,'matrixType','fact');
      %end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SKalmanMeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = SKalmanMeasurementUpdate(predMean,predSVar,h,Input,vMean,vSvar,Time,Measurement,H,Delta)
    % SKALMANMEASUREMENTUPDATE Implements measurement update of the square-root (E)KF.
      

      % measurement prediction mean and prediction error
      zPredMean = evaluate(h,predMean,Input,vMean,Time); % zp
      dz = Measurement-zPredMean;

      % measurement prediction variance
      Szp = [H*predSVar Delta*vSvar];

      % cross-covariance predictive matrix of x and z
      Pxzp = predSVar*predSVar'*H';

      % evaluating filtering mean and variance
      K = Pxzp/Szp'/Szp;
      filtMean = predMean + K*dz;
      Sf = nefDD1.triag([predSVar-K*H*predSVar, K*Delta*vSvar]);
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
    % SKalmanRTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = SKalmanRTSSmoothUpdate(initMean,initVarfac,filtMean,filtVarfac,predMean,predVarfac,F,Gamma,wMean,Sq,f,Time,Input) %#ok<INUSD,INUSL>
    % SKALMANRTSSMOOTHUPDATE Implements smoothing step by means of 
    % the square-root Rauch-Tung-Striebel smoother.
      
      Sf = filtVarfac;

      % covariance matrix Pxx = E{(x(k+1)-x_pred(k+1))(x(k)-x_filt(k))}
      Pxx = Sf*Sf'*F';

      % smoother gain
      Kv = (Pxx/predVarfac')/predVarfac;

      % evaluating smoothed mean and covariance matrix
      smoothMean = filtMean - Kv * (predMean - initMean);
      Sv = nefDD1.triag([Sf-Kv*F*Sf, Kv*Gamma*Sq, Kv*initVarfac]);
      smoothVar = nefCholFactorFunction(Sv,'matrixType','fact');
    end %SmoothUpdate

  end % methods
end % the class
