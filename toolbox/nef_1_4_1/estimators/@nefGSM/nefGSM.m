classdef nefGSM < nefEstimator
  %file @nefGSM/nefGSM.m (Gaussian Sum Filter)
  % nefGSM Properties:
  %   nefGSM_optsDef - definition of the nefGSM options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefGSM Methods:
  %   nefGSM - class constructor (For list of options type 'help nefGSM/nefGSM')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step distinguishing system description by equations or PDF's
  %   timeUpdateEqSystem - time update for system described by equations
  %   timeUpdatePDFSystem - time update for system described by PDF's
  %   measurementUpdate - measurement update step distinguishing system description by equations or PDF's
  %   measurementUpdateEqSystem - measurement update for system described by equations
  %   measurementUpdatePDFSystem - measurement update for system described by PDF's
  %   GSRVwithCholFactor - creation of Gaussian sum variable with covariances in form of factors
  %   GSRVwithUDFactor - creation of Gaussian sum variable with covariances in form of UD factors
  %
  % References:
  % Daniel L. Alspach, Harold W. Sorenson (1972):
  %   Nonlinear Bayesian estimation using Gaussian sum approximations,
  %   IEEE Trans. Automat. Contr., vol. 17, pp. 439 - 448

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % definition of estimator options
    nefGSM_optsDef = {2,'','localFilter','kalman','type of local filter',@(x) any(strcmpi(x,{'kalman','udkalman','skalman','ukf','dd1','dd2','sukf','sdd1','sdd2'}));
                      2,'','scalingParameter',NaN,'constant scaling parameter', @(x) isnumeric(x) && length(x)==1;
                      2,'','pruningMethod','none','mixtures pruning method',@(x) any(strcmpi(x,{'none','Threshold','MaxPosteriorProbability','CumulativeWeights','GenPseudoBayesian','DistanceBased'}));
                      2,'pruningMethod:Threshold','percentOfMaxWeight',NaN,'Preserve weights higher than specified percents of maximal weight value', @(x) isnumeric(x) && length(x)==1;
                      2,'pruningMethod:MaxPosteriorProbability','numberOfHighestWeight',NaN,'Number of preserved terms in the mixture with highest posterior probability)', @(x) isnumeric(x) && length(x)==1;
                      2,'pruningMethod:CumulativeWeights','cumulativeWeights',NaN,'Preserve the terms with highest weights which cumulative sum is less then or equal cumulativeWeights', @(x) isnumeric(x) && length(x)==1;
                      2,'pruningMethod:DistanceBased','qualityThreshold',NaN,'Sets the acceptance ratio for preserving the terms in the mixture', @(x) isnumeric(x) && length(x)==1;
                      2,'','sqrt','chol','matrix square-root used',@(x) any(strcmpi(x,{'chol','svd'}));
                      };
  end % constant properties


  properties (SetAccess = 'protected') % protected properties
    %localFilter@char; % type of local filter
    %scalingParameter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % true for PDF systems, false for Eq systems
    isPDFSystem
    % algoritms used for filtering and prediction
    localFilterFilt
    localFilterPred
    % following options apply only for nef EqSystem
    % the time-invariant state and measurement noise
    isStateNoiseVarConst
    isMeasNoiseVarConst
    % U and D factors of measurement noise variance for UD type filters
    Dr
    invUr
    % auxiliary variables for noises (e.g. for factorised versions, if necessary)
    wSVars

    vSVars
    vUVars
    vDVars
    % decomposition function for covariance matrix in UKF
    matrixDecompositionFunction;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefGSM(system,varargin)
      % NEFGSM Creates NEFGSM object.
      %
      %   OBJ = NEFGSM(system,varargin) creates a NEFGSM object OBJ
      %   representing Gaussian sum filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %
      %  PARAMETER                    DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =====================================================================================================================================
      %  nefEqSystem:localFilter      kalman           chosen local filter                  kalman,udkalman,skalman,ukf,dd1,dd2,sukf,sdd1,sdd2
      %  nefPDFSystem:localFilter     kalman           chosen local filter                  kalman,ukf,dd1,dd2
      %  scalingParameter             depends on       constant scaling parameter           if any then real, depends on local filter
      %                               localFilter
      %  pruningMethod                none             mixtures pruning method              MaxPosteriorProbability, Threshold, CumulativeWeights,
      %                                                                                     GenPseudoBayesian, Distancebased
      %  numberOfHighestWeight        NaN              Preserve this number of the          real number
      %                                                terms with highest weights
      %  percentOfMaxWeight           NaN              Preserve terms with weights higher   real number
      %                                                than given percents of maximal weight
      %  cumulativeWeights            NaN              Preserve the terms with highest      real number
      %                                                weights which cumulative sum is
      %                                                less then or equal cumulativeWeights
      %  qualityThreshold             NaN              The acceptance ratio for preserving  real number
      %                                                the terms in the mixture
      %
      %   See also ESTIMATE.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFGSM';
      p.addRequired('system',@(x) isa(x, 'nefSystem'));
      for i = 1:size(nefGSM.nefGSM_optsDef,1)
        p.addParamValue(nefGSM.nefGSM_optsDef{i,3},[],nefGSM.nefGSM_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefGSM.nefGSM_optsDef,p.Results);

      obj.isPDFSystem = isa(system,'nefPDFSystem'); % determine isPDFSystem attribute

      if obj.isPDFSystem % checks for nefPDFSystem
        if ~isa(system.x,'nefGaussianSumRV')
          error('NEF:nefGSM:xnotGaussianSum','State must be nefGaussianSumRV object')
        end
        if ~isa(system.z,'nefGaussianSumRV')
          error('NEF:nefGSM:znotGaussianSum','Measurement must be nefGaussianSumRV object')
        end

        %if ~any(strcmpi(obj.localFilter,{'kalman','ukf','dd1','dd2'}))
        if ~any(strcmpi(getOption(obj,'localFilter'),{'kalman','ukf','dd1','dd2'}))
          error('NEF:nefGSM:invalidFilterforPDF','Cannot use this local filter for nefPDFSystem model')
        end

      else % checks for nefEqSystem
        if ~isa(system.w,'nefGaussianSumRV')
          error('NEF:nefGSM:wnotGaussianSum','State noise must be nefGaussianSumRV object')
        end
        if ~isa(system.v,'nefGaussianSumRV')
          error('NEF:nefGSM:vnotGaussianSum','Measurement noise must be nefGaussianSumRV object')
        end

      end
      % common checks
      if ~isa(obj.x0,'nefGaussianSumRV')
        error('NEF:nefGSM:x0notGaussianSum','Initial condition must be nefGaussianSumRV object')
      end

      % setting matrix Decomposition function according to the options
      switch getOption(obj,'sqrt')
        case 'chol',
          obj.matrixDecompositionFunction = @(x) chol(x,'lower');
        case 'svd',
          obj.matrixDecompositionFunction = @(x) nefUKF.SVDDecomp(x);
      end;

      % scaling parameter check
      setOption(obj,'scalingParameter',nefItFilter.ScalingParameterCheck(obj))

      % prepare parameters according to the system class
      if obj.isPDFSystem
      else % the following lines are applicable for nefEqSystem only
        % computing square roots or UD of initial condition
        if any(strcmp(getOption(obj,'localFilter'),{'sukf','sdd1','skalman','sdd2'}))
          [obj.x0] = nefGSM.GSRVwithCholFactor(obj.x0,[],[],[],[]);
        elseif strcmp(getOption(obj,'localFilter'), 'udkalman')
          [obj.x0] = nefGSM.GSRVwithUDFactor(obj.x0,[],[],[],[]);
        end

        % check of linearity (if linear measurement or state function then appropriate part of Kalman filter used instead of chosen one)
        % prediciton part
        if  obj.sys.f.isLinear && ~any(strcmp(getOption(obj,'localFilter'),{'kalman','skalman','udkalman'}))
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
        if  obj.sys.h.isLinear && ~any(strcmp(getOption(obj,'localFilter'),{'kalman','skalman','udkalman'}))
          switch getOption(obj,'localFilter')
            case {'ukf', 'dd1', 'dd2'}
              obj.localFilterFilt = 'kalman';
            case {'sukf', 'sdd1', 'sdd2'}
              obj.localFilterFilt = 'skalman';
          end
        else
          obj.localFilterFilt = getOption(obj,'localFilter');
        end

        % if diff1Noise is Constant and measurement noise covariances are constant then Vtilde is also constant
        obj.isMeasNoiseVarConst =  obj.sys.h.isDiff1NoiseConst & obj.sys.v.VNumeric;
        % if diff1Noise is Constant and state noise covariances are constant then Vtilde is also constant
        obj.isStateNoiseVarConst =  obj.sys.f.isDiff1NoiseConst & obj.sys.w.VNumeric;
      end

      % predictors and, filters require prediction only in the queue
      data.RV = obj.x0;
      obj.estQueue(end+1) = data;
    end % nefGSM constructor

    function disp(obj)
      % DISP Displays nefGSM object parameters and options.

      fprintf('A nefGSM object with parameters\n')
      showOptions(obj)
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time)
      % TIMEUPDATE Distinguishes time update step according to system
      % description by means of equations or PDF's.

      if obj.isPDFSystem
        [predEstimate] = timeUpdatePDFSystem(obj,filtEst,Input,Time);
      else
        [predEstimate] = timeUpdateEqSystem(obj,filtEst,Input,Time);
      end
    end %timeUpdate function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATEEQSYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdateEqSystem(obj,filtEst,Input,Time)
      % TIMEUPDATEEQSYSTEM Performs time update step for systems described
      % by equations. Time update can be performed by means of any local
      % filter.

      % parameters of filtering or previous prediction estimate
      % parameters of prediction estimate
      % parameters of state noise
      switch obj.localFilterPred
        case {'kalman','dd1','dd2','ukf'}
          [wW,wM,wV] = evalParameters(obj.sys.w,[],Input,[],Time);
          [filtW,filtM,filtV] = evalParameters(filtEst.RV);
          filtN = filtEst.RV.N;
        case {'sukf','sdd1','skalman','sdd2'}
          if ((obj.isStateNoiseVarConst == 1) && (Time == 0)) || (obj.isStateNoiseVarConst == 0)
            obj.wSVars = evalSVars(obj.sys.w,[1:obj.sys.w.N],[],Input,[],Time);
          end
          wW = evalWeights(obj.sys.w,[1:obj.sys.w.N],[],Input,[],Time);
          wM = evalMeans(obj.sys.w,[1:obj.sys.w.N],[],Input,[],Time);

          filtN = filtEst.RV.N;
          filtW = evalWeights(filtEst.RV,[1:filtN]);
          filtM = evalMeans(filtEst.RV,[1:filtN]);
          filtSV = evalSVars(filtEst.RV,[1:filtN]);
        case 'udkalman'
          [wW,wM,wV] = evalParameters(obj.sys.w,[],Input,[],Time);

          filtN = filtEst.RV.N;
          filtW = evalWeights(filtEst.RV,[1:filtN]);
          filtM = evalMeans(filtEst.RV,[1:filtN]);
          [filtUV filtDV] = evalUDVars(filtEst.RV,[1:filtN]);
      end

      % scaling parameter determinantion
      if any(strcmp(obj.localFilterPred,{'ukf','dd1','dd2','sukf','sdd1','sdd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      pN = filtN*obj.sys.w.N;   % #filtering terms
      pW = zeros(1,pN);      % create new weights cell array
      pM = cell(1,pN);      % create new means cell array
      pV = cell(1,pN);      % create new variances cell array

      for i=1:pN,
        j=i-fix((i-1)/filtN)*filtN;
        n=1+fix((i-1)/filtN);

        pW(i) = filtW(j)*wW(n);


        % evaluating prediction mean and variance
        switch obj.localFilterPred
          case 'kalman'
            % matrix F for state in state equation
            F = evalDiff1State(obj.sys.f,filtM(:,j),Input,wM(:,n),Time);
            % matrix Gamma for noise in state equation
            Gamma = evalDiff1Noise(obj.sys.f,filtM(:,j),Input,wM(:,n),Time);
            [pM{i}, pV{i}] = nefKalman.KalmanTimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.f,Input,wM(:,n),wV(:,:,n),Time,F,Gamma);
          case 'ukf'
            [pM{i},pV{i}] = nefUKF.UKFTimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.f,Input,wM(:,n),wV(:,:,n),Time,sParameter,obj.matrixDecompositionFunction);
          case 'dd1'
            [pM{i},pV{i}] = nefDD1.DD1TimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.f,Input,wM(:,n),wV(:,:,n),Time,sParameter);
          case 'dd2'
            [pM{i},pV{i}] = nefDD2.DD2TimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.f,Input,wM(:,n),wV(:,:,n),Time,sParameter);
          case 'sukf'
            [pM{i},pV{i}] = nefSUKF.SUKFTimeUpdate(filtM(:,j),filtSV(:,:,j),obj.sys.f,Input,wM(:,n),obj.wSVars(:,:,n),Time,sParameter);
          case 'sdd1'
            [pM{i},pV{i}] = nefSDD1.SDD1TimeUpdate(filtM(:,j),filtSV(:,:,j),obj.sys.f,Input,wM(:,n),obj.wSVars(:,:,n),Time,sParameter);
          case 'sdd2'
            [pM{i},pV{i}] = nefSDD2.SDD2TimeUpdate(filtM(:,j),filtSV(:,:,j),obj.sys.f,Input,wM(:,n),obj.wSVars(:,:,n),Time,sParameter);
          case 'skalman'
            % matrix F for state in state equation
            F = evalDiff1State(obj.sys.f,filtM(:,j),Input,obj.sys.w.Means{n},Time);
            % matrix Gamma for noise in state equation
            Gamma = evalDiff1Noise(obj.sys.f,filtM(:,j),Input,obj.sys.w.Means{n},Time);
            [pM{i},pV{i}] = nefSKalman.SKalmanTimeUpdate(filtM(:,j),filtSV(:,:,j),obj.sys.f,Input,wM(:,n),obj.wSVars(:,:,n),Time,F,Gamma);
          case 'udkalman'
            % matrix F for state in state equation
            F = evalDiff1State(obj.sys.f,filtM(:,j),Input,obj.sys.w.Means{n},Time);
            % matrix Gamma for noise in state equation
            Gamma = evalDiff1Noise(obj.sys.f,filtM(:,j),Input,obj.sys.w.Means{n},Time);
            [pM{i},pV{i}] = nefUDKalman.UDKalmanTimeUpdate(filtM(:,j),filtUV(:,:,j),filtDV(:,:,j),obj.sys.f,Input,wM(:,n),Gamma*wV(:,:,n)*Gamma',Time,F);
        end
      end

      % normalize weights
      pW = nefRV.normalizeWeights(pW);
      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianSumRV(num2cell(pW),pM,pV,'check',0,'parameters','wmv');
    end % function timeUpdateEqSystem

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATEPDFSYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdatePDFSystem(obj,filtEst,Input,Time)
      % TIMEUPDATEEQSYSTEM Performs time update step for systems described
      % by PDF's. Time update can be performed by means of the UKF, DD1, and DD2 only.


      % parameters of filtering or previous prediction estimate
      [filtW,filtM,filtV] = evalParameters(filtEst.RV);
      filtN = filtEst.RV.N;

      pN = filtN*obj.sys.x.N;   % #filtering terms
      pW = zeros(1,pN);      % create new weights cell array
      pM = cell(1,pN);      % create new means cell array
      pV = cell(1,pN);      % create new variances cell array

      if any(strcmp(getOption(obj,'localFilter'),{'ukf','dd1','dd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      for i=1:pN,
        j=i-fix((i-1)/filtN)*filtN;
        n=1+fix((i-1)/filtN);
        xW = evalWeights(obj.sys.x,n,filtM(:,j),Input,[],Time);
        pW(i) = filtW(j)*xW;

        % evaluating prediction mean and variance
        xV = evalVars(obj.sys.x,n,filtM(:,j),Input,[],Time);
        wM = zeros(obj.sys.x.dimRV,1);
        switch getOption(obj,'localFilter')
          case 'kalman'
            % matrix F for state in state equation
            if obj.sys.x.MNumeric
              F = zeros(obj.sys.x.dimRV);
            else
              F = evalDiff1State(obj.sys.x.Means{n},filtM(:,j),Input,[],Time);
            end
            % matrix Gamma for noise in state equation
            Gamma = eye(obj.sys.x.dimRV);
            [pM{i}, pV{i}] = nefKalman.KalmanTimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.x.Means{n},Input,[],xV,Time,F,Gamma);
          case 'ukf'
            [pM{i},pV{i}] = nefUKF.UKFTimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.x.Means{n},Input,wM,xV,Time,sParameter,obj.matrixDecompositionFunction);
          case 'dd1'
            [pM{i},pV{i}] = nefDD1.DD1TimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.x.Means{n},Input,wM,xV,Time,sParameter);
          case 'dd2'
            [pM{i},pV{i}] = nefDD2.DD2TimeUpdate(filtM(:,j),filtV(:,:,j),obj.sys.x.Means{n},Input,wM,xV,Time,sParameter);
        end
      end

      % normalize weights
      pW = nefRV.normalizeWeights(pW);
      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianSumRV(num2cell(pW),pM,pV,'check',0,'parameters','wmv');
    end % function timeUpdatePDFSystem

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
      % MEASUREMENTUPDATE Distinguishes measurement update step according to system
      % description by means of equations or PDF's.

      if obj.isPDFSystem
        [filtEstimate] = measurementUpdatePDFSystem(obj,predEst,Input,Measurement,Time);
      else
        [filtEstimate] = measurementUpdateEqSystem(obj,predEst,Input,Measurement,Time);
      end
    end %measurementUpdate function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATEEQSYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdateEqSystem(obj,predEst,Input,Measurement,Time)
      % MEASUREMENTUPDATEEQSYSTEM Performs measurement update step for systems described
      % by equations. Measurement update can be performed by means of any local
      % filter.

      % parameters of prediction estimate
      % parameters of measurement noise
      switch obj.localFilterFilt
        case{'kalman','dd1','dd2','ukf'}
          [vW,vM,vV] = evalParameters(obj.sys.v,[],Input,[],Time);
          [predW,predM,predV] = evalParameters(predEst.RV);
          predN = predEst.RV.N;
        case {'sukf','sdd1','skalman','sdd2'}
          if ((obj.isMeasNoiseVarConst == 1) && (Time == 0)) || (obj.isMeasNoiseVarConst == 0)
            obj.vSVars = evalSVars(obj.sys.v,[1:obj.sys.v.N],[],Input,[],Time);
          end
          vW = evalWeights(obj.sys.v,[1:obj.sys.v.N],[],Input,[],Time);
          vM = evalMeans(obj.sys.v,[1:obj.sys.v.N],[],Input,[],Time);

          predN = predEst.RV.N;
          predW = evalWeights(predEst.RV,[1:predN]);
          predM = evalMeans(predEst.RV,[1:predN]);
          predSV = evalSVars(predEst.RV,[1:predN]);
        case 'udkalman'
          % UD factors are computed later (because of structure of input of nefUDKalman)
          [vW,vM,vV] = evalParameters(obj.sys.v,[],Input,[],Time);
          
          predN = predEst.RV.N;
          predW = evalWeights(predEst.RV,[1:predN]);
          predM = evalMeans(predEst.RV,[1:predN]);
          [predUV predDV] = evalUDVars(predEst.RV,[1:predN]);
      end

      if any(strcmp(obj.localFilterFilt,{'ukf','dd1','dd2','sukf','sdd1','sdd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      fN = predN*obj.sys.v.N;   % #filtering terms
      fW = zeros(1,fN);     % create new weights vector
      fM = cell(1,fN);      % create new means cell array
      fV = cell(1,fN);      % create new variances cell array

      for j = 1:fN
        i=j-fix((j-1)/predN)*predN;
        m=1+fix((j-1)/predN);


        switch obj.localFilterFilt
          case 'kalman'
            % matrix H for measurement in output equation
            H = evalDiff1State(obj.sys.h,predM(:,i),Input,vM(:,m),Time);
            % matrix Delta for noise in output equation
            Delta = evalDiff1Noise(obj.sys.h,predM(:,i),Input,vM(:,m),Time);
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefKalman.KalmanMeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.h,Input,vM(:,m),vV(:,:,m),Time,Measurement,H,Delta);
          case 'ukf'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefUKF.UKFMeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.h,Input,vM(:,m),vV(:,:,m),Time,Measurement,sParameter,obj.matrixDecompositionFunction);
          case 'dd1'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefDD1.DD1MeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.h,Input,vM(:,m),vV(:,:,m),Time,Measurement,sParameter);
          case 'dd2'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefDD2.DD2MeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.h,Input,vM(:,m),vV(:,:,m),Time,Measurement,sParameter);
          case 'sukf'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefSUKF.SUKFMeasurementUpdate(predM(:,i),predSV(:,:,i),obj.sys.h,Input,vM(:,m),obj.vSVars(:,:,m),Time,Measurement,sParameter);
          case 'sdd1'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefSDD1.SDD1MeasurementUpdate(predM(:,i),predSV(:,:,i),obj.sys.h,Input,vM(:,m),obj.vSVars(:,:,m),Time,Measurement,sParameter);
          case 'sdd2'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefSDD2.SDD2MeasurementUpdate(predM(:,i),predSV(:,:,i),obj.sys.h,Input,vM(:,m),obj.vSVars(:,:,m),Time,Measurement,sParameter);
          case 'skalman'
            % matrix H for measurement in output equation
            H = evalDiff1State(obj.sys.h,predM(:,i),Input,obj.sys.v.Means{m},Time);
            % matrix Delta for noise in output equation
            Delta = evalDiff1Noise(obj.sys.h,predM(:,i),Input,obj.sys.v.Means{m},Time);
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefSKalman.SKalmanMeasurementUpdate(predM(:,i),predSV(:,:,i),obj.sys.h,Input,vM(:,m),obj.vSVars(:,:,m),Time,Measurement,H,Delta);
          case 'udkalman'
            % matrix H for measurement in output equation
            H = evalDiff1State(obj.sys.h,predM(:,i),Input,obj.sys.v.Means{m},Time);
            % matrix Delta for noise in output equation
            Delta = evalDiff1Noise(obj.sys.h,predM(:,i),Input,obj.sys.v.Means{m},Time);
            if ((obj.isMeasNoiseVarConst==1) && (Time ==0)) || (obj.isMeasNoiseVarConst == 0)
                factored = nefUDFactorFunction(Delta*vV(:,:,m)*Delta');
                obj.Dr = factored.D;
                obj.invUr = inv(factored.U);
            end;
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefUDKalman.UDKalmanMeasurementUpdate(predM(:,i),predUV(:,:,i),predDV(:,:,i),obj.sys.h,Input,vM(:,m),obj.Dr,obj.invUr,Time,Measurement,H);
        end
        % compute weights
        dz = (Measurement-zPredMean);
        exponent = -0.5*dz'*invzPredVar*dz;
        zeta  = 1/(2*pi)^(obj.sys.dimMeasurement/2)/sqrt(detzPredVar)*exp(exponent);
        fW(j) = predW(i)*vW(m)*zeta;
        %logW(j) = exponent - 0.5*log(detzPredVar) ;
      end

      % normalize weights
      fW = nefRV.normalizeWeights(fW);
      % remove weight on the verge of numerical accuracy
      [cleanedWeights,cleanedMeans,cleanedVariances] = nefGaussianSumRV.pruningT(num2cell(fW),fM,fV,eps);

      % apply chosen pruning method
      [prunedWeights,prunedMeans,prunedVariances] = pruneMixturePDF(obj,cleanedWeights,cleanedMeans,cleanedVariances);

      % returning nefGaussianRV with symmetrized variance
      filtEstimate.RV = nefGaussianSumRV(prunedWeights,prunedMeans,prunedVariances,'check',0,'parameters','wmv');
    end % function measurementUpdateEqSystem

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATEPDFSYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdatePDFSystem(obj,predEst,Input,Measurement,Time)
      % MEASUREMENTUPDATEPDFSYSTEM Performs measurement update step for systems described
      % by PDF's. Measurement update can be performed by means of the UKF, DD1, and DD2 only.

      % parameters of prediction estimate
      [predW,predM,predV] = evalParameters(predEst.RV);
      predN = predEst.RV.N;

      fN = predN*obj.sys.z.N;   % #filtering terms
      fW = zeros(1,fN);      % create new weights vector
      fM = cell(1,fN);      % create new means cell array
      fV = cell(1,fN);      % create new variances cell array

      if any(strcmp(getOption(obj,'localFilter'),{'ukf','dd1','dd2'}))
        sParameter = getOption(obj,'scalingParameter');
      end

      for j = 1:fN
        i=j-fix((j-1)/predN)*predN;
        m=1+fix((j-1)/predN);

        zV = evalVars(obj.sys.z,m,predM(:,i),Input,[],Time);
        vM = zeros(obj.sys.z.dimRV,1);
        switch getOption(obj,'localFilter')
          case 'kalman'
            if obj.sys.z.MNumeric
              H = zeros(obj.sys.z.dimRV);
            else
              H = evalDiff1State(obj.sys.z.Means{m},predM(:,i),Input,[],Time);
            end
            Delta = eye(obj.sys.z.dimRV);
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefKalman.KalmanMeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.z.Means{m},Input,[],zV,Time,Measurement,H,Delta);
          case 'ukf'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefUKF.UKFMeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.z.Means{m},Input,vM,zV,Time,Measurement,sParameter,obj.matrixDecompositionFunction);
          case 'dd1'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefDD1.DD1MeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.z.Means{m},Input,vM,zV,Time,Measurement,sParameter);
          case 'dd2'
            [fM{j},fV{j},zPredMean,invzPredVar,detzPredVar] = nefDD2.DD2MeasurementUpdate(predM(:,i),predV(:,:,i),obj.sys.z.Means{m},Input,vM,zV,Time,Measurement,sParameter);
        end
      end

      % compute weights
      dz = (Measurement-zPredMean);
      exponent = -0.5*dz'*invzPredVar*dz;
      zeta  = 1/(2*pi)^(obj.sys.dimMeasurement/2)/sqrt(detzPredVar)*exp(exponent);
      zW = evalWeights(obj.sys.z,m,predM(:,i),Input,[],Time);
      fW(j) = predW(i)*zW*zeta;
      % normalize weights
      fW = nefRV.normalizeWeights(fW);
      % remove weight on the verge of numerical accuracy
      [cleanedWeights,cleanedMeans,cleanedVariances] = nefGaussianSumRV.pruningT(num2cell(fW),fM,fV,eps);

      % apply chosen pruning method
      [prunedWeights,prunedMeans,prunedVariances] = pruneMixturePDF(obj,cleanedWeights,cleanedMeans,cleanedVariances);

      % returning nefGaussianRV with symmetrized variance
      filtEstimate.RV = nefGaussianSumRV(prunedWeights,prunedMeans,prunedVariances,'check',0,'parameters','wmv');
    end % function measurementUpdatePDFSystem
  end %methods

  methods  (Access = 'private')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Prune the terms in the mixture given by GaussianSumRV
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [prunedWeights,prunedMeans,prunedVariances] = pruneMixturePDF(obj,Weights,Means,Vars)
          %pruneMixturePDF prunes the terms in the mixture given by
          %GaussianSumRV using selected method
          switch getOption(obj,'pruningMethod')
              case 'threshold'
                  [prunedWeights,prunedMeans,prunedVariances] = nefGaussianSumRV.pruningT(Weights,Means,Vars,getOption(obj,'percentOfMaxWeight'));
              case 'maxposteriorprobability'
                  [prunedWeights,prunedMeans,prunedVariances] = nefGaussianSumRV.pruningMPP(Weights,Means,Vars,getOption(obj,'numberOfHighestWeight'));
              case 'cumulativeweights'
                  [prunedWeights,prunedMeans,prunedVariances] = nefGaussianSumRV.pruningCW(Weights,Means,Vars,getOption(obj,'cumulativeWeights'));
              case 'genpseudobayesian'
                  [prunedWeights,prunedMeans,prunedVariances] = nefGaussianSumRV.pruningGPB1(Weights,Means,Vars);
              case 'distancebased'
                  [prunedWeights,prunedMeans,prunedVariances] = nefGaussianSumRV.pruningDBP(Weights,Means,Vars,getOption(obj,'qualityThreshold'));
              otherwise
                  prunedWeights = Weights;
                  prunedMeans = Means;
                  prunedVariances = Vars;
          end % switch
      end % pruneMixturePDF
  end % Private methods

  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GaussianSumRV to GaussianSumRV with CholFactorFunction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [y] = GSRVwithCholFactor(x,varargin)
    % GSRVWITHCHOLFACTOR creates Gaussian sum random variable with
    % covariances in form of their Cholesky factors.

        [Wx,Mx,Vx] = evalParameters(x,varargin{:});
        Vy = cell(1,x.N);
        % creating objects of NefCholFactorFunctions
        for i=1:x.N
          Vy{i} = nefCholFactorFunction(Vx(:,:,i),'matrixType','cov');
        end
        y = nefGaussianSumRV(num2cell(Wx),num2cell(Mx,1),Vy,'parameters','wmv');
    end % function GSRVwithCholFactor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GaussianSumRV to GaussianSumRV with UDFactorFunction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [y] = GSRVwithUDFactor(x,varargin)
    % GSRVWITHUDFACTOR creates Gaussian sum random variable with
    % covariances in form of their UD factors.

        [Wx,Mx,Vx] = evalParameters(x,varargin{:});
        Vy = cell(1,x.N);
        % creating objects of NefUDFactorFunctions
        for i=1:x.N
          Vy{i} = nefUDFactorFunction(Vx(:,:,i));
        end
        y = nefGaussianSumRV(num2cell(Wx),num2cell(Mx,1),Vy,'parameters','wmv');
    end % function GSRVwithCholFactor

  end % Static methods
end % the class
