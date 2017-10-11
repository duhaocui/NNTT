classdef nefItFilter < nefEstimator
  %file @nefItFilter/nefItFilter.m (Iterated Kalman Filter)
  % nefItFilter Properties:
  %   nefItFilter_optsDef - definition of the nefItFilter options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefItFilter Methods:
  %   nefItFilter - class constructor (For list of options type 'help nefItFilter/nefItFilter')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   ScalingParameterCheck - check of chosen scaling parameter and warning spec.
  %
  % References:
  % A.H. Jazwinski (1970):
  % Stochastic Processes and filtering Theory, Academic Press

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % definition of estimator options
    nefItFilter_optsDef = {2,'','localFilter','kalman','type of local filter',@(x) any(strcmpi(x,{'kalman','udkalman','skalman','ukf','dd1','dd2','sukf','sdd1','sdd2'}));
                           0,'','scalingParameter',NaN,'constant scaling parameter', @(x) isfloat(x) && length(x)==1;
                           2,'','epsilon',0.1,'maximum RMSE for stopping iteration step',@(x) isfloat(x) && x>0;
                           2,'','maxIter',100,'maximum number of iterations',@(x) isfloat(x) && x>0;
                           2,'','sqrt','chol','matrix square-root used',@(x) any(strcmpi(x,{'chol','svd'}));
                           };
  end % constant properties

  properties (SetAccess = 'protected') % protected properties
    %localFilter@char; % type of local filter
    %scalingParameter
    %maxIter % maximum number of iterations
    %epsilon

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time-invariant state and mesurement noise
    isWConst
    isVConst
    % U and D factors of measurement noise variance for UD type filters
    Dr
    invUr
    % S factors of state and measurement noise variance for S type fitlers
    Sq
    Sr
    % used algoritms for filtering and prediction
    localFilterFilt
    localFilterPred
    % decomposition function for covariance matrix in UKF
    matrixDecompositionFunction;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefItFilter(system,varargin)
      % NEFITFILTER Creates NEFITFILTER object.
      %
      %   OBJ = NEFITFILTER(system,varargin) creates a NEFITFILTER object OBJ
      %   representing iteration filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  localFilter                 kalman           chosen local filter                  kalman,udkalman,skalman,ukf,dd1,dd2,sukf,sdd1,sdd2
      %                                               constant or adaptively chosen
      %  scalingParameter            depends on       constant scaling parameter           if any then real, depends on local filter
      %                              localFilter
      %  epsilon                     0.1              maximum RMSE for stopping            real, positive
      %                                               iteration step
      %  maxIter                     100              maximum number of                    real, positive
      %                                               allowed iterations
      %   See also ESTIMATE.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFITFILTER';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for i = 1:size(nefItFilter.nefItFilter_optsDef,1)
        p.addParamValue(nefItFilter.nefItFilter_optsDef{i,3},[],nefItFilter.nefItFilter_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefItFilter.nefItFilter_optsDef,p.Results);

      % scaling parameter check
      %sParameter = nefItFilter.ScalingParameterCheck(obj);
      %setOption(obj,'scalingParameter',sParameter);
      setOption(obj,'scalingParameter',nefItFilter.ScalingParameterCheck(obj))


      % computing square roots or UD of initial condition
      if any(strcmp(getOption(obj,'localFilter'), {'sukf', 'sdd1', 'sdd2', 'skalman'}))
          obj.x0 = nefGaussianRV(obj.x0.evalMean,nefCholFactorFunction(evalVariance(obj.x0),'matrixType','cov'),'check',0);
      elseif strcmp(getOption(obj,'localFilter'), 'udkalman')
          obj.x0 = nefGaussianRV(obj.x0.evalMean,nefUDFactorFunction(evalVariance(obj.x0)),'check',0);
      end

      % check of linearity (if linear measurement or state function then appropriate part of Kalman filter used instead of chosen one)
      % prediciton part
      if  obj.sys.f.isLinear && ~any(strcmp(getOption(obj,'localFilter'),{'kalman','skalman', 'udkalman'}))
        switch getOption(obj,'localFilter')
          case {'ukf', 'dd1', 'dd2'}
            obj.localFilterPred = 'kalman';
          case {'sukf', 'sdd1', 'sdd2'}
            obj.localFilterPred = 'skalman';
        end
      else
        obj.localFilterPred = getOption(obj,'localFilter');
      end
      % filtering part
      if  obj.sys.h.isLinear && ~any(strcmpi(getOption(obj,'localFilter'),{'kalman', 'skalman', 'udkalman'}))
        switch getOption(obj,'localFilter')
          case {'ukf', 'dd1', 'dd2'}
            obj.localFilterFilt = 'kalman';
          case {'sukf', 'sdd1', 'sdd2'}
            obj.localFilterFilt = 'skalman';
        end
      else
        obj.localFilterFilt = getOption(obj,'localFilter');
      end

      % setting matrix Decomposition function according to the options
      switch getOption(obj,'sqrt')
        case 'chol',
          obj.matrixDecompositionFunction = @(x) chol(x,'lower');
        case 'svd',
          obj.matrixDecompositionFunction = @(x) nefUKF.SVDDecomp(x);
      end;

      % if diff1Noise is Constant and measurement noise covariances are constant then Vtilde is also constant
      obj.isVConst =  obj.sys.h.isDiff1NoiseConst & obj.sys.v.VarNumeric;
      % if diff1Noise is Constant and state noise covariances are constant then Vtilde is also constant
      obj.isWConst =  obj.sys.f.isDiff1NoiseConst & obj.sys.w.VarNumeric;

      % predictors and, filters require prediction only in the queue
      data.RV = obj.x0;
      obj.estQueue(end+1) = data;
    end % kalman constructor

    function disp(obj)
      % DISP Displays nefItFilter object parameters and options.

      fprintf('A nefItFilter object with parameters\n')
      showOptions(obj)
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time)
    % TIMEUPDATE Performs time update step with respect to system
    % description. Time update can be performed by means of any local
    % filter.

      % parameters of filtering or previous prediction estimate
      filtMean = evalMean(filtEst.RV,[],[],[],[]);
      filtVar = filtEst.RV.Var;
      %filtVar = evalVariance(filtEst.RV,[],[],[],[]);
      % mean and variance of state noise
      wMean = evalMean(obj.sys.w,[],Input,[],Time);
      wVar = evalVariance(obj.sys.w,[],Input,[],Time);

      % parameters of state noise
      if any(strcmp(getOption(obj,'localFilter'),{'sukf','sdd1','sdd2','skalman'}))
          if ((obj.isWConst == 1) && (Time == 0)) || (obj.isWConst == 0)
              % computing square roots
              factor = nefCholFactorFunction(wVar,'matrixType','cov');
              obj.Sq = factor.S;
          end
      elseif strcmp(getOption(obj,'localFilter'), 'udkalman')
          wVar = evalVariance(obj.sys.w,[],Input,[],Time);
      end

      % scaling parameter determination
      if any(strcmp(obj.localFilterPred,{'ukf','dd1','dd2','sukf','sdd1','sdd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      switch obj.localFilterPred
        case 'kalman'
          % matrices F for state and Gamma for noise in state equation
          F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
          Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
          [predMean,predVar] = nefKalman.KalmanTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,F,Gamma);
        case 'ukf'
          [predMean, predVar] = nefUKF.UKFTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,sParameter,obj.matrixDecompositionFunction); % scalingParameter = kappa
        case 'dd1'
          [predMean,predVar] = nefDD1.DD1TimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,sParameter); % scalingParameter = hpar
        case 'dd2'
          [predMean,predVar] = nefDD2.DD2TimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,sParameter); % scalingParameter = hpar
        case 'sukf'
          [predMean, predVar] = nefSUKF.SUKFTimeUpdate(filtMean,evaluateS(filtVar,[],[],[],[]),obj.sys.f,Input,wMean,obj.Sq,Time,sParameter); % scalingParameter = kappa
        case 'sdd1'
          [predMean, predVar] = nefSDD1.SDD1TimeUpdate(filtMean,evaluateS(filtVar,[],[],[],[]),obj.sys.f,Input,wMean,obj.Sq,Time,sParameter); % scalingParameter = hpar
        case 'sdd2'
          [predMean, predVar] = nefSDD2.SDD2TimeUpdate(filtMean,evaluateS(filtVar,[],[],[],[]),obj.sys.f,Input,wMean,obj.Sq,Time,sParameter); % scalingParameter = hpar
        case 'skalman'
          % matrices F for state and Gamma for noise in state equation
          F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
          Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
          [predMean, predVar] = nefSKalman.SKalmanTimeUpdate(filtMean,evaluateS(filtVar,[],[],[],[]),obj.sys.f,Input,wMean,obj.Sq,Time,F,Gamma);
        case 'udkalman'
          F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
          Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
          [predMean,predVar] = nefUDKalman.UDKalmanTimeUpdate(filtMean,filtVar.U,filtVar.D,obj.sys.f,Input,wMean,Gamma*wVar*Gamma',Time,F);
      end
      % returning nefGaussianRV
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
    % MEASUREMENTUPDATE Performs measurement update step with respect to
    % system description. Measurement update can be performed by means of any local
    % filter.

      % parameters of prediction estimate
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = predEst.RV.Var;
      %predVar = evalVariance(predEst.RV,[],[],[],[]);
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);

      % parameters of measurement noise
      if any(strcmp(getOption(obj,'localFilter'),{'sukf','sdd1','sdd2','skalman'}))
          if ((obj.isVConst == 1) && (Time == 0)) || (obj.isVConst == 0)
              % computing square roots
              factor = nefCholFactorFunction(vVar,'matrixType','cov');
              obj.Sr = factor.S;
          end
      elseif strcmp(getOption(obj,'localFilter'), 'udkalman')
          if ((obj.isVConst == 1) && (Time == 0)) || (obj.isVConst == 0)
              % UD factors
              Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
              %Rtilde = Delta*vVar*Delta';
              factor = nefUDFactorFunction(Delta*vVar*Delta');
              %obj.Ur = factor.U;
              obj.Dr = factor.D;
              obj.invUr = inv(factor.U);
          end
      end

      % scaling parameter determinantion
      if any(strcmp(obj.localFilterFilt,{'ukf','dd1','dd2','sukf','sdd1','sdd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      aEpsilon = getOption(obj,'epsilon');
      aMaxIter = getOption(obj,'maxIter');
      meanDiff = aEpsilon + 1; % to make at least one iteration
      iterCount = 1;
      filtMeanLastIter = predMean;
      while (meanDiff > aEpsilon) && (iterCount < aMaxIter)
        switch obj.localFilterFilt
          case 'kalman'
            % matrices H for state and Delta for noise in state equation
            H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
            Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
            [filtMean,filtVar] = nefKalman.KalmanMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,H,Delta);
          case 'ukf'
            [filtMean,filtVar] = nefUKF.UKFMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,sParameter,obj.matrixDecompositionFunction); % scalingParameter = kappa
          case 'dd1'
            [filtMean,filtVar] = nefDD1.DD1MeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,sParameter); % scalingParameter = hpar
          case 'dd2'
            [filtMean,filtVar] = nefDD2.DD2MeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,sParameter); % scalingParameter = hpar
          case 'sukf'
            [filtMean,filtVar] = nefSUKF.SUKFMeasurementUpdate(predMean,evaluateS(predVar,[],[],[],[]),obj.sys.h,Input,vMean,obj.Sr,Time,Measurement,sParameter); % scalingParameter = kappa
          case 'sdd1'
            [filtMean,filtVar] = nefSDD1.SDD1MeasurementUpdate(predMean,evaluateS(predVar,[],[],[],[]),obj.sys.h,Input,vMean,obj.Sr,Time,Measurement,sParameter); % scalingParameter = hpar
          case 'sdd2'
            [filtMean,filtVar] = nefSDD2.SDD2MeasurementUpdate(predMean,evaluateS(predVar,[],[],[],[]),obj.sys.h,Input,vMean,obj.Sr,Time,Measurement,sParameter); % scalingParameter = hpar
          case 'skalman'
            % matrices H for state and Delta for noise in state equation
            H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
            Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
            [filtMean,filtVar] = nefSKalman.SKalmanMeasurementUpdate(predMean,evaluateS(predVar,[],[],[],[]),obj.sys.h,Input,vMean,obj.Sr,Time,Measurement,H,Delta);
          case 'udkalman'
              H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
             [filtMean,filtVar] = nefUDKalman.UDKalmanMeasurementUpdate(predMean,predVar.U,predVar.D,obj.sys.h,Input,vMean,obj.Dr,obj.invUr,Time,Measurement,H);
        end
        meanDiff = sqrt(sum((filtMeanLastIter-filtMean).^2));
        filtMeanLastIter = filtMean;
        iterCount = iterCount+1;
      end
      % returning nefGaussianRV
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
    end % function measurementUpdate
  end %methods
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ScalingParameterCheck
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sParameter = ScalingParameterCheck(obj)
      % SCALINGPARAMETERCHECK Checks the chosen scaling parameter against
      % the system description.

      if any(strcmp(getOption(obj,'localFilter'),{'kalman','skalman','udkalman'}))
        % no parameter is required
        if ~isnan(getOption(obj,'scalingParameter'))
          disp('**************')
          disp('nefItFilter(nefGSM) WARNING: Scaling parameter is redundant for Kalman-filter-type filters. See documentation.')
          disp('**************')
        end
        sParameter = getOption(obj,'scalingParameter');
      elseif any(strcmp(getOption(obj,'localFilter'),{'dd1','dd2','sdd1','sdd2'}))
        if isnan(getOption(obj,'scalingParameter')) % parameter wasn't entered
          sParameter = sqrt(3);
          disp(['nefItFilter(nefGSM): Option constant scaling parameter ("scalingParameter") unspecified, using default value: "',num2str(sParameter),'"'])
        elseif any(strcmp(getOption(obj,'localFilter'),{'dd1','sdd1'})) && (getOption(obj,'scalingParameter')<=0) % was entered but wrong
          sParameter = sqrt(3);
          disp('**************')
          disp('nefItFilter(nefGSM) WARNING: Scaling parameter for DD1-type filter must be positive. Default value set. See documentation.')
          disp('**************')
        elseif any(strcmp(getOption(obj,'localFilter'),{'dd2','sdd2'})) && (getOption(obj,'scalingParameter')<=1) % was entered but wrong
          sParameter = sqrt(3);
          disp('**************')
          disp('nefItFilter(nefGSM) WARNING: Scaling parameter for DD2-type filter must be greater then 1. Default value set. See documentation.')
          disp('**************')
        else % was entered and correct
          sParameter = getOption(obj,'scalingParameter');
        end
      elseif any(strcmp(getOption(obj,'localFilter'),{'ukf','sukf'}))
        if isnan(getOption(obj,'scalingParameter')) % parameter wasn't entered
          sParameter = max(3-obj.sys.dimState,0);
          disp(['nefItFilter(nefGSM): Option constant scaling parameter ("scalingParameter") unspecified, using default value: "',num2str(sParameter),'"'])
        elseif (getOption(obj,'scalingParameter')<=0)
          sParameter = max(3-obj.sys.dimState,0);
          disp('**************')
          disp('nefItFilter(nefGSM) WARNING: Scaling parameter for UKF-type filter must be positive. Default value set. See documentation.')
          disp('**************')
        else
          sParameter = getOption(obj,'scalingParameter');
        end
      end
    end % function ScalingParameterCheck
  end % methods static
end % the class
