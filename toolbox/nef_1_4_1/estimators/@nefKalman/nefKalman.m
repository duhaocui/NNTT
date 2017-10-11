classdef nefKalman < nefEstimator
  %file @nefKalman/nefKalman.m ((extended) Kalman Filter and Rauch-Tung-Striebel Smoother)
  % nefKalman Methods:
  %   nefKalman - class constructor
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   KalmanTimeUpdate - time update of the (E)KF
  %   KalmanMeasurementUpdate - measurement update of the (E)KF
  %   RTSSmoothUpdate - smoothing based on the (extended) Rauch-Tung-Striebel smoother
  %
  % References:
  % A.H. Jazwinski (1970):
  % Stochastic Processes and filtering Theory, Academic Press

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia
  
  properties (SetAccess = 'protected') % protected properties
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefKalman(system,varargin) 
      % NEFKALMAN Creates NEFKALMAN object.
      %
      %   OBJ = NEFKALMAN(system,varargin) creates a NEFKALMAN object OBJ 
      %   representing Kalman filter or Rauch-Tung-Striebel smoother for model
      %   SYSTEM.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFKALMAN';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end
      obj.estQueue(end+1) = data;
    end % kalman constructor
    
    function disp(obj)
      % DISP Displays nefKalman object parameters and options.
      
      fprintf('A nefKalman object with parameters\n')
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
      % matrices F for state and Gamma for noise in state equation
      F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
      Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);

      % evaluating prediction mean and variance
      [predMean,predVar] = nefKalman.KalmanTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,F,Gamma);

      % returning nefGaussianRV
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);
      % return also by-products if required

      % is timUpdate is used for smoothing
      if p.Results.smoothingPurpose
        predEstimate.auxData.F = F; % return matrix F also as an additional field in predEstimate
      end % if
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
    % MEASUREMENTUPDATE Performs measurement update step with respect to system description.
      
      % parameters of prediction estimate
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = evalVariance(predEst.RV,[],[],[],[]);
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);
      % matrices H for state and Delta for noise in state equation
      H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
      Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);

      % evaluating filtering mean and variance
      [filtMean,filtVar] = nefKalman.KalmanMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,H,Delta);

      % returning nefGaussianRV
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
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
      F = predEstimate.auxData.F;

      [smoothM,smoothV] = nefKalman.RTSKalmanSmoothUpdate(initM,initV,filtM,filtV,predM,predV,F);

      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate
  end %methods
  methods(Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RTSKalmanSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = RTSKalmanSmoothUpdate(initMean,initVar,filtMean,filtVar,predMean,predVar,F)
      % RTSSMOOTHUPDATE Implements smoothing step by means of 
      % Rauch-Tung-Striebel smoother.
      
      Kv = filtVar * F' /predVar;
      smoothMean = filtMean - Kv * (predMean - initMean);
      smoothVar = filtVar - Kv * (predVar - initVar) * Kv';
      smoothVar = (smoothVar+smoothVar')/2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % KalmanMeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = KalmanMeasurementUpdate(predMean,predVar,h,Input,vMean,R,Time,Measurement,H,Delta)
      % KALMANMEASUREMENTUPDATE Implements measurement update of the (E)KF.
      
      numval = size(predMean,2);
      %preallocation
      zPredVar = zeros(size(Measurement,1),size(Measurement,1),numval);
      invzPredVar = zeros(size(Measurement,1),size(Measurement,1),numval);
      filtMean = zeros(size(predMean,1),numval);
      filtVar = zeros(size(predMean,1),size(predMean,1),numval);

      I = eye(size(predMean,1));
      zPredMean = evaluate(h,predMean,Input,vMean,Time);
      for in = 1:numval
        Rtilde = Delta(:,:,in)*R(:,:,in)*Delta(:,:,in)';
        dz = Measurement(:,in)-zPredMean(:,in);
        zPredVar(:,:,in) = Rtilde + H(:,:,in)*predVar(:,:,in)*H(:,:,in)';
        invzPredVar(:,:,in) = inv(zPredVar(:,:,in));
        K = predVar(:,:,in)*H(:,:,in)'*invzPredVar(:,:,in);
        filtMean(:,in) = predMean(:,in) + K*dz;
        filtVar(:,:,in) = (I - K * H(:,:,in)) * predVar(:,:,in) * (I - K * H(:,:,in))' + K * Rtilde * K';
        filtVar(:,:,in) = (filtVar(:,:,in)+filtVar(:,:,in)')/2;
      end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        % used for GSM filter
        varargout{3} = zPredMean;
        varargout{4} = invzPredVar;
        detzPredVar = zeros(1,numval);
        for in = 1:numval
          detzPredVar(in) = det(zPredVar(:,:,in));
        end
        varargout{5} = detzPredVar;
      end
    end %function KalmanFiltering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % KalmanTimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean,predVar] = KalmanTimeUpdate(filtMean,filtVar,f,Input,wMean,Q,Time,F,Gamma)
    % KALMANTIMEUPDATE Implements time update of the (E)KF.
      
      numval = size(filtMean,2);
      predMean = evaluate(f,filtMean,Input,wMean,Time);
      predVar =  zeros(size(filtMean,1),size(filtMean,1),numval);
      for in = 1:numval
        Qtilde = Gamma(:,:,in)*Q(:,:,in)*Gamma(:,:,in)';
        predVar(:,:,in) = F(:,:,in)*filtVar(:,:,in)*F(:,:,in)' + Qtilde;
        predVar(:,:,in) = (predVar(:,:,in) + predVar(:,:,in)')/2;
      end
    end % function KalmanTimeUpdate
  end % static methods
end % the class
