classdef nefPerformanceEvaluator < handle
  %file @nefPerformanceEvaluator/nefPerformanceEvaluator.m

  % ABSOLUTE ERROR MEASURES
  % MSEM - mean squared error matrix (state error)
  % RMSE - root mean squared error (state error)
  % AEE  - average Euklidean error (state error)
  % HAE  - harmonic average error (state error)
  % GAE  - geometric average error (state error)
  % MEDE - median error (state error)
  % MODE - mode error (state error)
  % RELATIVE ERROR MEASURES
  % RMSRE - root mean squared relative error (state error)
  % ARE   - average Euclidean relative error (state error)
  % BEEQ  - Bayesian estimation error quotient (state error, predicted state error)
  % EMER  - estimation error relative to measurement error (state estimate, predicted state estimate, measurement, function h) 
  % PERFORMANCE MEASURES
  % NCI - noncredibility index
  % ANEES - average normalized estimation error squared
  %

  % References:
  % X. R. Li and Z.-L. Zhao,
  % “Evaluation of Estimation Algorithms—Part I: Incomprehensive Performance Measures,” 
  % IEEE Transactions on Aerospace and Electronic Systems,
  %
  % X. R. Li, Z.-L. Zhao and V. P. Jilkov
  % Practical Measures and Test for Credibility of an Estimator
  % Proc. Workshop on Estimation, Tracking, and Fusion - A Tribute
  % to Yaakov Bar-Shalom, Monterey, CA, May 2001, pp. 481-495
  %
  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ASSIGNABLE VALUES AND DEFAULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    method@char % evaluation method
    dimState % state dimension (or parts to be measured)
    idxState % index of the state components to be measured
    dimMeasurement % measurement dimension
    timeSteps % #time steps
    MCRuns % #Monte Carlo runs
    storedMCRuns = 0; % #stored Monte Carlo runs
    system % system object
    nEstimators % # estimators to evaluate
    data % data
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefPerformanceEvaluator(system,timeSteps,MCRuns,nEstimators,varargin)
      % NEFPERFORMANCEEVALUATOR Creates NEFPERFORMANCEEVALUATOR object.
      %
      %   OBJ = NEFPERFORMANCEEVALUATOR(SYSTEM,TIMESTEPS,MCRUNS,NESTIMATORS,VARARGIN)
      %   creates a NEFPERFORMANCEEVALUATOR object OBJ
      %   representing a performance evaluator for SYSTEM and reserving storage space
      %   for TIMESTEPS of time steps, MCRUNS of Monte Carlo runs, NESTIMATORS of number
      %   of estimators. The performance measure can be chosen by specifying a property 'method'
      %   with values:
      %   ABSOLUTE ERROR MEASURES:
      %   'MSEM'  - mean squared error matrix (state error)
      %   'RMSE'  - (default) root mean squared error (state error)
      %   'AEE'   - average Euklidean error (state error)
      %   'HAE'   - harmonic average error (state error)
      %   'GAE'   - geometric average error (state error)
      %   'MEDE'  - median error (state error)
      %   'MODE'  - mode error (state error)
      %   RELATIVE ERROR MEASURES:
      %   'RMSRE' - root mean squared relative error (state error)
      %   'ARE'   - average Euclidean relative error (state error)
      %   'BEEQ'  - Bayesian estimation error quotient (state error, predicted state error)
      %   'EMER'  - estimation error relative to measurement error (state estimate, predicted state estimate, measurement, function h) 
      %   PERFORMANCE MEASURES:
      %   'NCI'   - noncredibility index
      %   'ANEES' - average normalized estimation error squared
      %   The index of the state elements to be used for evaluation can be specified via property
      %   'idxState' and value consisting of a vector of state indices
      %
      %   See also PROCESSDATA, PERFORMANCEVALUE


      methodValues = {'MSEM','RMSE','AEE','HAE','GAE','MEDE','MODE','RMSRE','ARE','BEEQ','EMER','NCI','ANEES'};

      p = inputParser;
      p.FunctionName = 'NEFPERFORMANCEEVALUATOR';
      p.addRequired('system',@(x) isa(x, 'nefSystem'));
      p.addRequired('timeSteps',@(x) isfloat(x) && (floor(x) == x) && (x > 0));
      p.addRequired('MCRuns',@(x) isfloat(x) && (floor(x) == x) && (x > 0));
      p.addRequired('nEstimators',@(x) isfloat(x) && (floor(x) == x) && (x > 0));
      p.addParamValue('method','RMSE',@(x) any(strcmpi(x,methodValues)));
      p.addParamValue('idxState',[1:system.dimState],@(x) all(isfloat(x)) && all(floor(x) == x) &&...
        all(x > 0) && all(x<=system.dimState) && (length(x)<=system.dimState));
      p.parse(system,timeSteps,MCRuns,nEstimators,varargin{:});

      obj.system = p.Results.system;
      obj.timeSteps = p.Results.timeSteps;
      obj.MCRuns = p.Results.MCRuns;
      obj.nEstimators = p.Results.nEstimators;
      obj.method = p.Results.method;
      obj.idxState = p.Results.idxState;


      obj.dimState = length(obj.idxState);
      obj.dimMeasurement = obj.system.dimMeasurement;
      switch obj.method
        case {'MSEM'}
          for h = 1:obj.nEstimators
            obj.data.MSEM{h} = zeros(obj.timeSteps,obj.dimState,obj.dimState);
          end
        case {'RMSE','AEE','HAE','GAE','MEDE','MODE'}
          for h = 1:obj.nEstimators
            obj.data.EucNorm{h} = zeros(obj.MCRuns,obj.timeSteps);
          end
        case {'RMSRE','ARE'}
          for h = 1:obj.nEstimators
            obj.data.EucNorm{h} = zeros(obj.MCRuns,obj.timeSteps);
            obj.data.reference{h} = zeros(obj.MCRuns,obj.timeSteps);
          end
        case {'BEEQ'}
          for h = 1:obj.nEstimators
            obj.data.EucNormFilt{h} = zeros(obj.MCRuns,obj.timeSteps);
            obj.data.EucNormPred{h} = zeros(obj.MCRuns,obj.timeSteps);
          end
        case {'EMER'}
          error('sorry, not yet fully implemented')
          if any(obj.idxState ~= [1:obj.sys.dimState])
            error('for EMER full state vector must be utilized')
          end
          for h = 1:obj.nEstimators
            obj.data.EucNormNum{h} = zeros(obj.MCRuns,obj.timeSteps);
            obj.data.EucNormDen{h} = zeros(obj.MCRuns,obj.timeSteps);
          end
        case {'NCI'}
          if MCRuns < 2
            error('At least 2 MCRUNS required for NCI')
          end
          for h = 1:obj.nEstimators
            obj.data.NEES{h} = zeros(obj.MCRuns,obj.timeSteps);
            obj.data.stateErrors{h} = zeros(obj.MCRuns,obj.timeSteps,obj.dimState);
          end
        case {'ANEES'}
          if MCRuns < 2
            error('At least 2 MCRUNS required for ANEES')
          end
          for h = 1:obj.nEstimators
            obj.data.NEES{h} = zeros(obj.MCRuns,obj.timeSteps);
          end
      end %switch
    end % nefPerformanceEvaluator

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % processData
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function processData(obj,trueData,estData)
      % PROCESSDATA Processes data obtained from Monte Carlo simulations
      %
      %   PROCESSDATA(OBJ,TRUEDATA,ESTDATA)
      %   processes data in the performance evaluator object OBJ.
      %   Suppose M simulations for K timesteps are available for N estimators
      %   Then
      %   ESTDATA is a three-dimensional NxMxK cell array of either:
      %   * nefRV objects provided by the estimators
      %   * structure with fields
      %         - stateMean
      %         - stateVariance (necessary only for NCI and ANEES)
      %   FOR  BEEQ:
      %   * structure with fields
      %         - filtering
      %         - prediction
      %     each containing either nefRV or a structure with a field stateMean
      %   TRUEDATA is a structure with fields 'state' (true state), 'measurement'
      %   (measurement) and 'input' (input)
      %   each being a two-dimensional MxK cell array of column
      %   vectors of proper dimensions.
      %   Note that the 'measurement' field is required for 'EMER' method only.
      %
      %   See also NEFPERFORMANCEEVALUATOR, PERFORMANCEVALUE

      % checking #trueData
      if ~isfield(trueData,'state')
        error('trueData must contain "state" field')
      end

      if size(trueData.state,1) > obj.MCRuns - obj.storedMCRuns
        error('first dimension of trueData.state is higher than storage space available')
      end

      stateMCRuns = size(trueData.state,1);

      if size(trueData.state,2) ~= obj.timeSteps
        size(estData)
        error('second dimension of trueData.state does not match the number of timeSteps')
      end

      if strcmp(obj.method,'EMER')
        if ~isfield(trueData,'measurement')
          error('for method EMER trueData must contain field "measurement"')
        end

        if size(trueData.measurement,1) ~= stateMCRuns
          error('inconsistent number of MCRuns for trueData.state and trueData.measurement')
        end

        if size(trueData.measurement,2) ~= obj.timeSteps
          error('second dimension of trueData.measurement does not match the number of timeSteps')
        end

        if ~isfield(trueData,'input')
          error('for method EMER trueData must contain field "input"')
        end

        if size(trueData.input,1) ~= stateMCRuns
          error('inconsistent number of MCRuns for trueData.state and trueData.measurement')
        end

        if size(trueData.input,2) ~= obj.timeSteps
          error('second dimension of trueData.measurement does not match the number of timeSteps')
        end
      end

      % checking estData
      if size(estData,1) ~= obj.nEstimators
        error('first dimension of estData does not match the number of estimators')
      end
      if size(estData,2) ~= stateMCRuns
        error('inconsistent number of MCRuns for trueData.state and estData')
      end
      if size(estData,3) ~= obj.timeSteps
        size(estData)
        error('third dimension of estData does not match the number of timeSteps')
      end

      newMCRuns = size(estData,2);

      %process estData
      switch obj.method
        case {'MSEM'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                stateError = trueData.state{i,j}(obj.idxState) - stateMean(obj.idxState); 
                obj.data.MSEM{h}(j,:,:) = obj.data.MSEM{h}(j,:,:) +...
                  reshape(stateError*stateError',1,obj.dimState,obj.dimState);
              end
            end
          end
        case {'RMSE','AEE','HAE','GAE','MEDE','MODE'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                stateError = trueData.state{i,j}(obj.idxState) - stateMean(obj.idxState);
                obj.data.EucNorm{h}(obj.storedMCRuns+i,j) = sqrt(stateError'*stateError);
              end
            end
          end
        case {'RMSRE','ARE'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                stateError = trueData.state{i,j}(obj.idxState) - stateMean(obj.idxState);
                obj.data.EucNorm{h}(obj.storedMCRuns+i,j) = sqrt(stateError'*stateError);
                obj.data.reference{h}(obj.storedMCRuns+i,j) = sqrt(trueData.state{i,j}(obj.idxState)'*...
                  trueData.state{i,j}(obj.idxState));
              end
            end
          end
        case {'BEEQ'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                % filtering error
                stateFiltMean = getMean(obj,estData{h,i,j}.filtering);
                stateFiltError = trueData.state{i,j}(obj.idxState) - stateFiltMean(obj.idxState);
                obj.data.EucNormFilt{h}(obj.storedMCRuns+i,j) = sqrt(stateFiltError'*stateFiltError);
                % prediction error
                statePredMean = getMean(obj,estData{h,i,j}.prediction);
                statePredError = trueData.state{i,j}(obj.idxState) - statePredMean(obj.idxState);
                obj.data.EucNormPred{h}(obj.storedMCRuns+i,j) = sqrt(statePredError'*statePredError);
              end
            end
          end
        case {'EMER'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                hx = sys.h(trueData.state{i,j},trueData.input{i,j},[],j-1);
                d_hx = hx - sys.h(stateMean,trueData.input{i,j},[],j-1);
                d_z = hx - trueData.measurement{i,j};
                obj.data.EucNormNum{h}(obj.storedMCRuns+i,j) = sqrt(d_hx'*d_hx);
                obj.data.EucNormDen{h}(obj.storedMCRuns+i,j) = sqrt(d_z'*d_z);
              end
            end
          end
        case {'NCI'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                stateError = trueData.state{i,j}(obj.idxState) - stateMean(obj.idxState);
                stateVariance = getVariance(obj,estData{h,i,j});
                obj.data.NEES{h}(obj.storedMCRuns+i,j) = stateError'/stateVariance(obj.idxState,obj.idxState)*stateError;
                obj.data.stateErrors{h}(obj.storedMCRuns+i,j,:) = stateError;
              end
            end
          end
        case {'ANEES'}
          for h = 1:obj.nEstimators
            for i = 1:newMCRuns
              for j = 1:obj.timeSteps
                stateMean = getMean(obj,estData{h,i,j});
                stateError = trueData.state{i,j}(obj.idxState) - stateMean(obj.idxState);
                stateVariance = getVariance(obj,estData{h,i,j});
                obj.data.NEES{h}(obj.storedMCRuns+i,j) = stateError'/stateVariance(obj.idxState,obj.idxState)*stateError;
              end
            end
          end
      end % switch

      % update MCRuns counter
      obj.storedMCRuns = obj.storedMCRuns + newMCRuns;
    end % function processData

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % performanceValue
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = performanceValue(obj,varargin)
      % PERFORMANCEVALUE Returns value of performance measure.
      %
      %   [VAL] = PERFORMANCEVALUE(OBJ,VARARGIN)
      %   Calculate value of the performance measure for object OBJ.
      %   The index of the estimators for which the value is requested
      %   can be specified by property 'estimatorID' with value given
      %   by a vector of estimator indices. By default performance 
      %   measure of all estimators is returned.
      %   
      %   VAL is then a NxK matrix of values with N being the number
      %   of estimators and K being the number of time steps
      %
      %   See also NEFPERFORMANCEEVALUATOR, PROCESSDATA

      estimatorID = getEstimatorID(obj,varargin{:});
      switch obj.method
        case {'MSEM'}
          for h = 1:length(estimatorID)
            for j = 1:obj.timeSteps
              val{h,j} = squeeze(obj.data.MSEM{h}(j,:,:))/obj.storedMCRuns;
            end
          end
        case {'RMSE'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:).^2).^0.5;
          end
        case {'AEE'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:));
          end
        case {'HAE'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:).^-1).^-1;
          end
        case {'GAE'}
          for h = 1:length(estimatorID)
            val(h,:) = exp(mean(log(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:))));
          end
        case {'MEDE'}
          for h = 1:length(estimatorID)
            val(h,:) = median(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:));
          end
        case {'MODE'}
        case {'RMSRE'}
          for h = 1:length(estimatorID)
            val(h,:) = mean((obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:).^2)...
              ./(obj.data.reference{estimatorID(h)}(1:obj.storedMCRuns,:).^2)).^0.5;
          end
        case {'ARE'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.EucNorm{estimatorID(h)}(1:obj.storedMCRuns,:)...
              ./obj.data.reference{estimatorID(h)}(1:obj.storedMCRuns,:));
          end
        case {'BEEQ'}
          for h = 1:length(estimatorID)
            % TODO check
            val(h,:) = mean(obj.data.EucNormPred{estimatorID(h)}(1:obj.storedMCRuns,:))./...
                       mean(obj.data.EucNormFilt{estimatorID(h)}(1:obj.storedMCRuns,:));
          end
        case {'EMER'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.EucNormNum{estimatorID(h)}(1:obj.storedMCRuns,:))...
              ./mean(obj.data.EucNormDen{estimatorID(h)}(1:obj.storedMCRuns,:));
          end
        case {'NCI'}
          for h = 1:length(estimatorID)
            eID = estimatorID(h); % for shorter notation
            for l = 1:obj.timeSteps
              % first calculating MSEM
              MSEM = zeros(obj.dimState);
              for i = 1:obj.storedMCRuns
                MSEM = MSEM + squeeze(obj.data.stateErrors{eID}(i,l,:))*squeeze(obj.data.stateErrors{eID}(i,l,:))';
              end
              % calculate MSEM inversion
              %iMSEM = inv(MSEM/obj.storedMCRuns);
              % calculate denominator
              for i = 1:obj.storedMCRuns
                NCIden(i,l) = obj.storedMCRuns*squeeze(obj.data.stateErrors{eID}(i,l,:))'/MSEM*squeeze(obj.data.stateErrors{eID}(i,l,:));
              end
            end
            val(h,:) = 10*mean(log10(obj.data.NEES{eID}(1:obj.storedMCRuns,:)))-...
              10*mean(log10(NCIden(1:obj.storedMCRuns,:)));
          end
        case {'ANEES'}
          for h = 1:length(estimatorID)
            val(h,:) = mean(obj.data.NEES{estimatorID(h)}(1:obj.storedMCRuns,:))./obj.dimState;
          end
      end % switch
    end % performanceValue

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getMean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [stateMean] = getMean(obj,estData)
      if isstruct(estData)
        stateMean = estData.stateMean;
      elseif isa(estData,'nefRV')
        stateMean = evalMean(estData);
      else
        error('estData element is neither a structure nor an nefRV object')
      end
    end % function getMean

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getVariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [stateVariance] = getVariance(obj,estData)
      if isstruct(estData)
        stateVariance = estData.stateVariance;
      elseif isa(estData,'nefRV')
        stateVariance = evalVariance(estData);
      else
        error('estData element is neither a structure nor an nefRV object')
      end
    end % function getVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getEstimatorID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = getEstimatorID(obj,varargin)
      p = inputParser;
      p.FunctionName = 'GETESTIMATORIDX';
      p.addParamValue('estimatorID',[1:obj.nEstimators],@(x) all(isfloat(x)) && all(floor(x) == x) && all(x > 0) && all(x<=obj.nEstimator));
      p.parse(varargin{:});

      val = p.Results.estimatorID;
    end % function getEstimatorID

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      % DISP Display object.
      fprintf('A nefPerformanceEvaluator object\n'); 
    end % function disp
  end % methods
end %class nefPerformanceEvaluator
