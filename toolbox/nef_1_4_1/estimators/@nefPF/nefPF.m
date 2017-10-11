classdef nefPF < nefEstimator
  % file @nefPF/nefPF.m
  % nefPF Properties:
  %   nefPF_optsDef - definition of the nefPF options
  %   isPDFSystem - true for system of nefPDFSystem class
  % nefPF Methods:
  %   nefPF - class constructor (For list of options type 'help nefPF/nefPF')
  %   disp - display class info
  %   INTERNAL METHODS
  %   checkSystemAgainstSD - check system against requested sampling density (internal)
  %   ekf_SD - EKF based sampling density (internal)
  %   logFunctionalPrimaryWeights - logarithm of primary weights for functional auxiliary sampling density (internal)
  %   logPointPrimaryWeights - logarithm of primary weights for point auxiliary sampling density (internal)
  %   auxiliary_SD - sample and weights for axutiliary sampling density (internal)
  %   optimal_SD - sample and weights for optimal sampling density (internal)
  %   prior_SD - sample and weights for prior sampling density (internal)
  %   resampling - resampled samples and weights (internal)
  %   initialMeasurementUpdate - measurement update for initial time instant (internal)
  %   timeMeasurementUpdate - com,bined time and measurement update (internal)
  %   timeUpdate - time update (internal)
  %   prediction - redefinition of nefEstimator prediction (internal)
  %   filtering - redefinition of nefEstimator filtering (internal)

  % References:
  % A. Doucet, N. de Freitas and N. Gordon (2001):
  % Sequential Monte Carlo Methods in Practice, Springer

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % Definition of nefPF options
    nefPF_optsDef = {2,'','sampleSize',100,'sample size', @(x) isfloat(x) && (floor(x) == x) && (x > 0);
    2,'','samplingDensity','prior','sampling density',@(x) any(strcmpi(x,{'prior','optimal','functionalAuxiliary','pointAuxiliary','ekf'}));
    2,'samplingDensity:pointAuxiliary','PASDPointEstimate','mean','point estimate in point auxiliary sampling density',@(x) any(strcmpi(x,{'mean','mode','median','sample'}));
    2,'samplingDensity:functionalAuxiliary','logFASDPrimaryWeights' ,[],'primary weights in functional euxiliary sampling density',@(x) isa(x,'function_handle');
    2,'','resamplingSched','static','resampling schedule',@(x) any(strcmpi(x,{'static','dynamic'}));
    2,'resamplingSched:dynamic','resamplingSchedDynPar' ,0.7,'parameter for dynamic resampling schedule',@(x) isfloat(x) && (0 <= x) && (x <= 1);
    2,'resamplingSched:static','resamplingSchedStatPar',1,'parameter for static resampling schedule',@(x) isfloat(x) && (floor(x) == x) && (x > 0);
    2,'','resamplingAlg','systematic','resampling algorithm',@(x) any(strcmpi(x,{'multinomial','residual','stratified','systematic'}))};
  end % constant properties
  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    isPDFSystem % Indicator of class of the system
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefPF(system,varargin)
      % NEFPF Creates NEFPF object.
      %
      %   OBJ = NEFPF(system,varargin) creates a NEFPF object OBJ 
      %   representing particle filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed 
      %   via the standard Parameter-Value MATLAB mechanism
      %  PARAMETER              DEFAULT VAL  DESCRIPTION                  VALID VALUES
      %  =============================================================================
      %  sampleSize             100          sample size                  positive integer
      %  samplingDensity        prior        sampling density             'prior','optimal','functionalAuxiliary','pointAuxiliary','ekf'
      %  PASDPointEstimate      mean         point estimate in PASD       'mean','mode','median','sample'
      %  logFASDPrimaryWeights  []           primary weights in FASD      function_handle
      %  resamplingSched        static       resampling schedule          'static','dynamic'
      %  resamplingSchedDynPar  0.7          par. for dynamic res. sched. 0 <= x <= 1
      %  resamplingSchedStatPar 1            par. for static res. sched.  positive integer
      %  resamplingAlg          systematic   resampling algorithm         'multinomial','residual','stratified','systematic'
      %
      %   See also ESTIMATE.

      %
      % since there may not be a one-step prediction at disposal  in PF,
      % each field in estQueue represents last filtering, i.e. it has to
      % contain a last input as well in a structure with last filtering
      %
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFPF';
      p.addRequired('system',@(x) isa(x, 'nefSystem'));
      for  i = 1:size(nefPF.nefPF_optsDef,1)
        p.addParamValue(nefPF.nefPF_optsDef{i,3},[],nefPF.nefPF_optsDef{i,6});
      end
      p.parse(system, varargin{:});

      % calling nefEstimator constructor to set up tasType, taskPar and estQueue
      % estQueue is filled but not used at time=0
      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      obj.isPDFSystem = isa(system,'nefPDFSystem');
      processOptions(obj,nefPF.nefPF_optsDef,p.Results)
      checkSystemAgainstSD(obj)
      data.RV = obj.x0;
      obj.estQueue(end+1) = data;
    end % pf constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function disp(obj)
      % DISP Displays nefPF object parameters and options.
      fprintf('A nefPF object with parameters\n')
      showOptions(obj)
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDEFINED NEFESTIMATOR FILTERING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = filtering(obj,Input,Measurement)
      % FILTERING Performs filtering process based on measurement and input.

      % expand empty Input to proper length
      if isempty(Input)
        Input = zeros(1,size(Measurement,2));
      end
      val = cell(1,size(Measurement,2));
      % get last Estimate from estQueue
      filtEstimate = obj.estQueue(end);      
      % do for each measurement
      for k = 1:size(Measurement,2)
        if obj.time == 0 
          filtEstimate = initialMeasurementUpdate(obj,Input(:,k),Measurement(:,k));
        else
          % resampling step
          resEstimate = resampling(obj,filtEstimate);
          % filtering step (time and measurement update)
          filtEstimate = timeMeasurementUpdate(obj,resEstimate,Input(:,k),Measurement(:,k));
        end
        % return RV element of filtEstimate'
        % at each time instant for sequential dataProcessing
        % at last time instant for batch dataProcessing
        if strcmp(getOption(obj,'dataProcessing'),'sequential')
          val{k} = filtEstimate.RV;
        else % batch dataProcessing
          if k == size(Measurement,2)
            val{1} = filtEstimate.RV;
          end
        end
        val{k} = filtEstimate.RV;
        % increase time
        obj.time = obj.time + 1;
      end
      % store filtering estimate
      obj.estQueue(end+1) = filtEstimate;
    end % function filtering

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDEFINED NEFESTIMATOR PREDICTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = prediction(obj,Input,Measurement)
      % PREDICTION Performs prediction process based on measurement and input.

      % expand empty Input to proper length
      if isempty(Input)
        Input = zeros(1,size(Measurement,2)+obj.taskPar-1);
      end
      % prepare empty cell array for results
      % for time = 0 , return also prediction p(x_lag|z^-1)
      if strcmp(getOption(obj,'dataProcessing'),'sequential')
        if obj.time == 0
          val = cell(1,size(Measurement,2)+1);
        else
          val = cell(1,size(Measurement,2));
        end
      else % batch dataProcessing
        val = cell(1,1);
      end
      j = 0; % pointer to val
      % before estimation
      % obtain filtering estimate from queue
      % for time = 0 it contains X0
      filtEstimate = obj.estQueue(end);

      % for time = 0 , return also prediction p(x_lag|z^-1)
      if obj.time == 0 % use x0 instead of lastEstimate
        %----------------------------------------------------------
        % return prediction estimate p(x_lag|z^-1)
        % at each time instant for sequential dataProcessing
        % at last time instant for batch dataProcessing (in this case (k=0) if no measurement is given i.e. k=0 is the last time instant to be processed)
        if strcmp(getOption(obj,'dataProcessing'),'sequential') || size(Measurement,2) == 0
          predEstimate = obj.x0;
          for i = 1:obj.taskPar-1
            predEstimate = timeUpdate(obj,predEstimate,Input(:,i),obj.time+i-1);
          end
          j = j+1;
          val{j} = predEstimate.RV;
        end
      end
      %----------------------------------------------------------
      for k = 1:size(Measurement,2)
        % compute filtering estimate p(x_time|z^time)
        %----------------------------------------------------------
        %----------------------------------------------------------
        %----------------------------------------------------------
        if obj.time == 0 % different MeasurementUpdate
          % filtering step
          filtEstimate = initialMeasurementUpdate(obj,Input(:,k),Measurement(:,k));
          %----------------------------------------------------------
          %----------------------------------------------------------
          %----------------------------------------------------------
        else
          % resampling step
          resEstimate = resampling(obj,filtEstimate);
          % filtering step
          filtEstimate = timeMeasurementUpdate(obj,resEstimate,Input(:,k),Measurement(:,k));
        end
        %----------------------------------------------------------
        % compute and return prediction estimate p(x_time+lag|z^time)
        % at each time instant for sequential dataProcessing or at last time instant  for batch dataProcessing
        if strcmp(getOption(obj,'dataProcessing'),'sequential') || (k == size(Measurement,2))
          predEstimate = filtEstimate;
          for i = 1:obj.taskPar
            predEstimate = timeUpdate(obj,predEstimate,Input(:,i+k),obj.time+i+k-1);
            j = j+1;
            val{j} = predEstimate.RV;
          end
        end
        %----------------------------------------------------------
        % increase time instant
        obj.time = obj.time + 1;
      end
      % after estimation
      % store filtering estimate
      obj.estQueue(end+1) = filtEstimate;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIME UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,est,Input,Time)
      % TIMEUPDATE Performs time update step.
      N = getOption(obj,'sampleSize');
      values = zeros(obj.sys.dimState,N);
      for i = 1:N
        values(:,i) = drawSample(obj.sys.x,N,est.RV.values(:,i),Input,[],Time);
      end
      predEstimate.RV = nefEmpiricalRV(values,est.RV.weights);
      predEstimate.Input = Input;
    end % function  timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIME MEASUREMENT TUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = timeMeasurementUpdate(obj,est,Input,Measurement)
      % TIMEMEASUREMENTUPDATE Performs both time and measurements steps.

      % TODO ADD ADAPTIVE SAMPLE SIZE
      % switch according to the specified samplingDensity
      switch getOption(obj,'samplingDensity')
        %-----------------------------------------------------------------
        %------------------------------------------------------------PRIOR
        %-----------------------------------------------------------------
        case 'prior'
          [values logWeights] = prior_SD(obj,est,Input,Measurement);
          %-------------------------------------------------------------------
          %------------------------------------------------------------OPTIMAL
          %-------------------------------------------------------------------
        case 'optimal'
          [values logWeights] = optimal_SD(obj,est,Input,Measurement);
          %---------------------------------------------------------------------------
          %----------------------------------------------------------------- AUXILIARY
          %---------------------------------------------------------------------------
        case {'functionalauxiliary','pointauxiliary'}
          [values logWeights] = auxiliary_SD(obj,est,Input,Measurement);
          %---------------------------------------------------------------------------
          %----------------------------------------------------------------- EXTENDED KALMAN
          %---------------------------------------------------------------------------
        case 'ekf'
          [values logWeights] = ekf_SD(obj,est,Input,Measurement);
      end % switch
      % store estimate and new Input and normalize weights
      filtEstimate.Input = Input;
      filtEstimate.RV = nefEmpiricalRV(values,nefRV.normalizeLogWeights(logWeights));
    end % function timeMeasurementUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL MEASUREMENT UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = initialMeasurementUpdate(obj,Input,Measurement)
      % INITIALMEASUREMENTUPDATE Performs initial mesurement update step.

      logWeights = zeros(1,getOption(obj,'sampleSize')); % for numerical stability
      values = drawSample(obj.x0,getOption(obj,'sampleSize'));
      if obj.isPDFSystem
        % calculate weights according to the measurement pdf in obj.sys.z
        for i = 1:getOption(obj,'sampleSize')
          logWeights(i) = evalLogPDF(obj.sys.z,Measurement,values(:,i),Input,[],obj.time);
        end
      else % Eqsystem
        %        for i = 1:getOption(obj,'sampleSize')
        % calculate corresponding weights using the new Input
        %          hx(:,i) = obj.sys.h(values(:,i),Input,zeros(size(obj.sys.v)),obj.time);
        %          logWeights(i) = evalLogPDF(obj.sys.v,Measurement-hx(:,i),[],Input,[],obj.time);
        % calculate corresponding weights using the new Input
        %        end
        % USE LOGLIKELIHOOD INSTEAD OF COMPUTING IT AGAIN
        logWeights = logLikelihood(obj.sys,repmat(Measurement,1,getOption(obj,'sampleSize')),values,repmat(Input,1,getOption(obj,'sampleSize')),repmat(obj.time,1,getOption(obj,'sampleSize')));
      end
      % store Input for the next time instant
      filtEstimate.RV = nefEmpiricalRV(values,nefRV.normalizeLogWeights(logWeights));
      filtEstimate.Input = Input;
    end % function initialMeasurementUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function resEstimate = resampling(obj,filtEstimate)
      % RESAMPLING Resamples samples and weights.

      proceedResampling = 0;
      % assign all values in filtEstimate and change
      % appropriate fields if resampling proceeds
      resEstimate = filtEstimate;
      % static resampling check
      if strcmp(getOption(obj,'resamplingSched'),'static')
        if mod(obj.time,getOption(obj,'resamplingSchedStatPar')) == 0
          proceedResampling = 1;
        end
        % dynamic resampling check
      elseif strcmp(getOption(obj,'resamplingSched'),'dynamic')
        if efficientSampleSize(filtEstimate.RV) < getOption(obj,'resamplingSchedDynPar')
          proceedResampling = 1;
        end
      end
      if proceedResampling
        switch getOption(obj,'resamplingAlg')
          case 'multinomial'
            resEstimate.RV = multinomialResampling(filtEstimate.RV);
          case 'residual'
            resEstimate.RV = residualResampling(filtEstimate.RV);
          case 'stratified'
            resEstimate.RV = stratifiedResampling(filtEstimate.RV);
          case 'systematic'
            resEstimate.RV = systematicResampling(filtEstimate.RV);
        end % switch
      end % if
    end % function resampling

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAMPLING DENSITIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIOR SAMPLING DENSITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values logWeights] = prior_SD(obj,est,Input,Measurement)
      % PRIOR_SD Samples from prior sampling density.

      N = getOption(obj,'sampleSize');
      values = zeros(obj.sys.dimState,N);
      logWeights = zeros(1,N);
      if obj.isPDFSystem
        for i = 1:N
          % draw new samples ( use old Input stored in est
          values(:,i) = drawSample(obj.sys.x,1,est.RV.values(:,i),est.Input,[],obj.time-1);
          % calculate corresponding weights using the new Input
          logWeights(i) = evalLogPDF(obj.sys.z,Measurement,values(:,i),Input,[],obj.time) + log(est.RV.weights(i));
        end
      else % EqSystem
        % draw new samples ( use old Input stored in est
        wSample = drawSample(obj.sys.w,N,[],est.Input,[],obj.time-1);
        values = obj.sys.f(est.RV.values,repmat(est.Input,1,N),wSample,repmat(obj.time-1,1,N));
        % calculate corresponding weights using the new Input
        logWeights = logLikelihood(obj.sys,repmat(Measurement,1,N),values,repmat(Input,1,N),repmat(obj.time,1,N)) + log(est.RV.weights);
      end
    end % function prior_SD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIMAL SAMPLING DENSITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values logWeights] = optimal_SD(obj,est,Input,Measurement)
      % OPTIMAL_SD Samples from optimal sampling density.

      N = getOption(obj,'sampleSize');
      values = zeros(obj.sys.dimState,N);
      logWeights = zeros(1,N);
      % TODO: extension for GaussianSum
      if obj.isPDFSystem
        % variance of the measurement noise v_k
        R = evalVariance(obj.sys.z,[],Input,[],obj.time);
        % Jacobian of h w.r.t noise
        Delta = eye(obj.sys.z.dimRV);
        for i = 1:N
          % predictive means of state x_k
          xPredMean = evalMean(obj.sys.x,est.RV.values(:,i),est.Input,[],obj.time-1);
          % variance of the state noise w_k-1
          Q = evalVariance(obj.sys.x,est.RV.values(:,i),est.Input,[],obj.time-1);
          % Jacobian of h
          H = evalDiff1State(obj.sys.z.Mean,xPredMean,Input,[],obj.time);
          % compute mean and covariance matrix
          [fM,fV] = nefKalman.KalmanMeasurementUpdate(xPredMean,Q,obj.sys.z.Mean,Input,[],R,obj.time,Measurement,H,Delta);
          % draw samples
          values(:,i) = fM + chol(fV)'*randn(obj.sys.dimState,1);
          % predictive means of measurement z_k
          zPredMean = evalMean(obj.sys.z,xPredMean,Input,[],obj.time);
          % compute weights
          logWeights(i) = nefGaussianRV.evalGaussianPDF(Measurement,zPredMean,H*Q*H',1,[],[],[]) + log(est.RV.weights(i));
        end
      else % EqSystem
        % mean of the state noise w_k-1
        wMean = evalMean(obj.sys.w,[],est.Input,[],obj.time-1);
        % variance of the state noise w_k-1
        Q = evalVariance(obj.sys.w,[],est.Input,[],obj.time-1);
        % mean of the measurement noise v_k
        vMean = evalVariance(obj.sys.v,[],Input,[],obj.time);
        % variance of the measurement noise v_k
        R = evalVariance(obj.sys.v,[],Input,[],obj.time);
        for i = 1:N
          % predictive mean of state x_k
          xPredMean = obj.sys.f(est.RV.values(:,i),est.Input,wMean,obj.time-1);
          % Jacobian of h w.r.t state
          H = evalDiff1State(obj.sys.h,xPredMean,Input,[],obj.time);
          % Jacobian of h w.r.t noise
          Delta = evalDiff1Noise(obj.sys.h,xPredMean,Input,[],obj.time);
          % compute mean and covariance matrix
          [fM,fV] = nefKalman.KalmanMeasurementUpdate(xPredMean,Q,obj.sys.h,Input,vMean,R,obj.time,Measurement,H,Delta);
          % draw samples
          values(:,i) = fM + chol(fV)'*randn(obj.sys.dimState,1);
          % predictive means of measurement z_k
          zPredMeans = obj.sys.h(xPredMean,Input,vMean,obj.time);
          % compute weights
          logWeights(i) = nefGaussianRV.evalGaussianPDF(Measurement,zPredMean,H*Q*H',1,[],[],[]) + log(est.RV.weights(i));
        end
      end
    end % function optimal_SD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIARY SAMPLING DENSITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values,logWeights] = auxiliary_SD(obj,est,Input,Measurement)
      % AUXILIARY_SD Samples from auxiliary sampling density.

      N = getOption(obj,'sampleSize');
      values = zeros(obj.sys.dimState,N);
      logWeights = zeros(1,N);
      % compute primary weights according to sampling density subtype
      switch getOption(obj,'samplingDensity')
        case 'functionalauxiliary'
          [logPW] = logFunctionalPrimaryWeights(obj,est,Input,Measurement);
        case 'pointauxiliary'
          [logPW] = logPointPrimaryWeights(obj,est,Input,Measurement);
      end % switch
      % compute indices for selecting samples
      [idx] = nefRV.rndMultinomial(nefRV.normalizeLogWeights(logPW),N);
      if obj.isPDFSystem
        for i = 1:N
          % draw new samples ( use old Input stored in est)
          values(:,i) = drawSample(obj.sys.x,1,est.RV.values(:,idx(i)),est.Input,[],obj.time-1);
          % calculate corresponding weights using the new Inout
          logWeights(i) = evalLogPDF(obj.sys.z,Measurement,values(:,i),Input,[],obj.time) + log(est.RV.weights(i)) - logPW(idx(i));
        end % for
      else % nefEqSystem
        % draw new samples ( use old Input stored in est
        wSamples = drawSample(obj.sys.w,N,[],est.Input,[],obj.time-1);
        values = evaluateOverAllStates(obj.sys.f,RV.values(:,idx),est.Input,wSample,obj.time-1);
        % calculate corresponding weights using the new Inout
        measPrediction = evaluateOverAllStates(obj.sys.h,values(:,i),Input,zeros(size(obj.sys.v)),obj.time);
        for i = 1:N
          logWeights(i) = evalLogPDF(obj.sys.v,Measurement-measPrediction(:,i),[],Input,[],obj.time) + log(est.RV.weights(i)) - logPW(idx(i));
        end % for
      end
    end % function auxiliary_SD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOG POINT AUXILIARY SAMPLING DENSITY WEIGHTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [primaryWeights] = logPointPrimaryWeights(obj,est,Input,Measurement)
      % LOGPOINTPRIMARYWEIGHTS Calculates primary weights for point auxiliary sampling density.

      N = getOption(obj,'sampleSize');
      pointEstimate = zeros(obj.sys.dimState,N);
      primaryWeights = zeros(1,N);
      if obj.isPDFSystem
        switch getOption(obj,'PASDPointEstimate')
          case 'mean'
            pointEstimate = evalMean(obj.sys.x,est.RV.values,est.Input(:,ones(1,N)),[],(obj.time-1)*ones(1,N));
          case 'mode'
            pointEstimate = evalMode(obj.sys.x,est.RV.values,est.Input(:,ones(1,N)),[],(obj.time-1)*ones(1,N));
          case 'median'
            pointEstimate = evalMedian(obj.sys.x,est.RV.values,est.Input(:,ones(1,N)),[],(obj.time-1)*ones(1,N));
          case 'sample'
            for i = 1:N
              pointEstimate(:,i) = drawSample(obj.sys.x,1,est.RV.values(:,i),est.Input,[],obj.time-1);
            end
        end % switch
        % use predicted point estimate to compute primary weights
        for i = 1:N
          primaryWeights(i) = evalLogPDF(obj.sys.z,Measurement,pointEstimate(:,i),Input,[],obj.time);
        end
      else % nefEgSystem
        wPointEstimate = zeros(obj.sys.dimState,N);
        switch getOption(obj,'PASDPointEstimate')
          case 'mean'
            wPointEstimate = repmat(evalMean(obj.sys.w,[],est.Input,[],obj.time-1),1,N);
          case 'mode'
            wPointEstimate = repmat(evalMode(obj.sys.w,[],est.Input,[],obj.time-1),1,N);
          case 'median'
            wPointEstimate = repmat(evalMedian(obj.sys.w,[],est.Input,[],obj.time-1),1,N);
          case 'sample'
            wPointEstimate = drawSample(obj.sys.w,N,[],est.Input,[],obj.time-1);
        end % switch
        % use predicted point estimate to compute primary weights
        pointEstimate = evaluate(obj.sys.f,est.RV.values,est.Input(:,ones(1,N)),wPointEstimate(:,i),obj.(time-1)*ones(1,N));
        for i = 1:N
          primaryWeights(i) = evalLogPDF(obj.sys.z,Measurement,pointEstimate(:,i),Input,[],obj.time);
        end
      end
    end %% function logPointPrimaryWeight

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOG FUNCTIONAL AUXILIARY SAMPLING DENSITY WEIGHTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [primary_weights] = logFunctionalPrimaryWeights(obj,est,Input,Measurement)
      % LOGFUNCTIONALPRIMARYWEIGHTS Calculates primary weights for functional auxiliary sampling density.

      if isempty(getOption(obj,'FASDPrimaryWeights'))
        error('NEF:PF:functionalSDPWUndef','Primary weights for functional sampling density undefined!');
      end
      fcn = getOption(obj,'logFASDPrimaryWeights');
      primary_weights = fcn(obj,est.RV.values,est.Input,obj.time-1,Measurement,Input,obj.time);
    end % function logFunctionalPrimaryWeights

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EKF SAMPLING DENSITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values logWeights] = ekf_SD(obj,est,Input,Measurement)
      % EKF_SD Samples from the EKF sampling density.

      if obj.isPDFSystem
        error('NEF:PF:EKFSD_PDFSYSTEM','For EKF SD, system must be a nefEqSystem object');
      else % EqSystem
        N = getOption(obj,'sampleSize'); % to shorten code
        logwpi = zeros(1,N);
        values = zeros(obj.sys.dimState,N);
        % mean of the state noise w_k-1
        wMean = evalMean(obj.sys.w,[],est.Input,[],obj.time-1);
        % variance of the state noise w_k-1
        predVar = repmat(evalVariance(obj.sys.w,[],est.Input,[],obj.time-1),[1,1,N]);
        % mean of the measurement noise v_k
        aux = evalMean(obj.sys.v,[],Input,[],obj.time);
        vMean = aux(:,ones(1,N));
        % variance of the measurement noise v_k
        vVar = repmat(evalVariance(obj.sys.v,[],Input,[],obj.time),[1,1,N]);

        % predictive mean of state x_k
        predMean = evaluateOverAllStates(obj.sys.f,est.RV.values,est.Input,wMean,obj.time-1);
        % matrices H for state and Delta for noise in state equation
        H = evalDiff1State(obj.sys.h,predMean,Input(:,ones(1,N)),vMean,obj.time(:,ones(1,N)));
        Delta = evalDiff1State(obj.sys.h,predMean,Input(:,ones(1,N)),vMean,obj.time(:,ones(1,N)));
        % evaluating filtering mean and variance

        [filtMean,filtVar] = nefKalman.KalmanMeasurementUpdate(predMean,predVar,obj.sys.h,Input(:,ones(1,N)),vMean,vVar,obj.time(:,ones(1,N)),Measurement(:,ones(1,N)),H,Delta);
        for i = 1:N
          % draw samples
          values(:,i) = filtMean(:,i) + chol(filtVar(:,:,i))'*randn(obj.sys.dimState,1);
          % compute weights
          % sampling pdf
          logwpi(i) = nefGaussianRV.evalGaussianPDF(values(:,i),filtMean(:,i),filtVar(:,:,i),1,[],[],[]) + log(est.RV.weights(i));
        end
        % transition pdf
        logwx = logTransitionPDF(obj.sys,values,est.RV.values,repmat(est.Input,1,N),repmat(obj.time-1,1,N));
        % measurement pdf
        logwz = logLikelihood(obj.sys,repmat(Measurement,1,N),values,repmat(Input,1,N),repmat(obj.time,1,N));

        logWeights = logwx + logwz - logwpi + log(est.RV.weights);
      end
    end % function ekf_SD

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END SAMPLING DENESITIES END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  end % methods
  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECKSYSTEMAGAINSTSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkSystemAgainstSD(obj)
      % CHECKSYSTEMAGAINSTSD Checks system and sampling density for compatibility.

      switch getOption(obj,'samplingDensity')
        %-----------------------------------------------------------------
        %------------------------------------------------------------PRIOR
        %-----------------------------------------------------------------
        case 'prior'
          if obj.isPDFSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % no restriction on system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          else %nefEqSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.sys must contain loglikelihood specification
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.sys.log_likelihood_handle must be specified
            if isempty (obj.sys.log_likelihood_handle)
              error('NEF:PF:LogLikelihoodUnspecified','For prior SD, loglikelihood of the system must be specified');
            end
          end
          %-------------------------------------------------------------------
          %------------------------------------------------------------OPTIMAL
          %-------------------------------------------------------------------
        case 'optimal'
          % TODO: extension for GaussianSum

          if obj.isPDFSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % transition pdf must be Gaussian with variance independent of
            % state
            % measurement pdf must be Gaussian with variance independent of
            % state and mean must be a linear function of the state
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % obj.sys.x must be nefGaussianRV
            if ~isa(obj.sys.x,'nefGaussianRV')
              error('NEF:PF:xNotGaussian','For optimal SD, state of the system must be a Gaussian random variable');
            end
            % variance of obj.sys.x can not depend on previous state
            %if ~obj.sys.x.VarNumeric & obj.sys.x.Var.dimState > 0
            %  error('NEF:PF:QdependsOnState','For optimal SD, variance of the state RV can not depend on previous state')
            %end
            % obj.sys.z must be nefGaussianRV
            if ~isa(obj.sys.z,'nefGaussianRV')
              error('NEF:PF:zNotGaussian','For optimal SD, measurement of the system must be a Gaussian RV');
            end
            % variance of obj.sys.z can not depend on previous state
            if ~obj.sys.z.VarNumeric && obj.sys.z.Var.dimState > 0
              error('NEF:PF:RdependsOnState','For optimal SD, variance of the measurement RV can not depend on state')
            end
            % mean of obj.sys.z must be a linear function
            if ~obj.sys.z.Mean.isLinear
              error('NEF:PF:hIsNotLinear','For optimal SD, mean of the measurement RV must be a linear function')
            end
          else %nefEqSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % state equation must be additive and state noise Gaussian
            % measurement equation must be linear, additive and measurement
            % noise Gaussian
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.sys.w must be nefGaussianRV
            if ~isa(obj.sys.w,'nefGaussianRV')
              error('NEF:PF:wNotGaussian','For optimal SD, state noise must be a Gaussian RV');
            end
            % obj.sys.f must be additive
            if obj.sys.f.isAdditive
              error('NEF:PF:fNotAdditive','For optimal SD, function f must be additive');
            end
            % obj.sys.v must be nefGaussianRV
            if ~isa(obj.sys.v,'nefGaussianRV')
              error('NEF:PF:vNotGaussian','For optimal SD, measurement noise must be a Gaussian RV');
            end
            % obj.sys.h must be linear and additive
            if obj.sys.h.isLinear && obj.sys.h.isAdditive
              error('NEF:PF:hNotLinearAdditive','For optimal SD, function h must be linear and additive3');
            end
          end
          %---------------------------------------------------------------------------
          %----------------------------------------------------------------- AUXILIARY
          %---------------------------------------------------------------------------
        case {'functionalAuxiliary','pointAuxiliary'}
          if obj.isPDFSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % no restriction on system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          else %nefEqSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % no restriction on system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          end
          %---------------------------------------------------------------------------
          %----------------------------------------------------------------- EXTENDED KALMAN
          %---------------------------------------------------------------------------
        case 'ekf'
          if obj.isPDFSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % unable to use ekf SD for nefPDFSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            error('NEF:PF:ekfSDforPDFSystem','Can not compute ekf SD for PDF system')
          else %nefEqSystem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.sys must contain log_likelihood and transitionPDF specification
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.sys.log_likelihood_handle must be specified
            if isempty (obj.sys.log_likelihood_handle)
              error('NEF:PF:LogLikelihoodUnspecified','For EKF SD, loglikelihood of the system must be specified');
            end
            if isempty (obj.sys.log_transitionPDF_handle)
              error('NEF:PF:LogTransitionPDFUnspecified','For EKF SD, logtransitionPDF of the system must be specified');
            end
          end
      end % switch
    end % function checkSystemAgainstSD
  end %methods private
end % the class
