classdef nefEnKF < nefEstimator
  % file @nefEnKF/nefEnKF.m
  % nefEnKF Properties:
  %   optsDef - first property
  %   opts - second property
  % nefEnKF Methods:
  %   nefEnKF - first method
  %

  % References:
  %
  %

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % definition of estimator options
    nefEnKF_optsDef = {2,'','ensembleSize',100,'ensemble size', @(x) isfloat(x) && (floor(x) == x) && (x > 0)};
  end % constant properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefEnKF(system,varargin)
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFENKF';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for  i = 1:size(nefEnKF.nefEnKF_optsDef,1)
        p.addParamValue(nefEnKF.nefEnKF_optsDef{i,3},[],nefEnKF.nefEnKF_optsDef{i,6});
      end
      p.parse(system, varargin{:});
      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % TODO test for additive measurement equation ??

      processOptions(obj,nefEnKF.nefEnKF_optsDef,p.Results)

      eS = getOption(obj,'ensembleSize');
      values = drawSample(obj.x0,eS);
      x0Ensemble.RV = nefEmpiricalRV(values,ones(1,eS)/eS);

      obj.estQueue(end+1) = x0Ensemble;
    end % nefEnKF constructor

    function disp(obj)
      fprintf('A nefEnKF object with parameters\n')
      showOptions(obj)
    end % function disp


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIME UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,est,Input,Time)
      eS = getOption(obj,'ensembleSize');
      wSample = drawSample(obj.sys.w,eS,[],est.Input,[],obj.time-1);
      predEnsemble = evaluate(obj.sys.f,est.RV.values,est.Input(:,ones(1,eS)),wSample,(obj.time-1)*ones(1,eS));
      predEstimate.RV = nefEmpiricalRV(predEnsemble,est.RV.weights);
      predEstimate.Input = Input;
    end % function  timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIME MEASUREMENT UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
      eS = getOption(obj,'ensembleSize');
      % parameters of prediction estimate
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = evalVariance(predEst.RV,[],[],[],[]);
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);
      % matrices H for state and Delta for noise in state equation
      H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
      Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);

      zPredMeans = evaluateOverAllStates(obj.sys.h,predEst.RV.values,Input,vMean,Time);
      perturbedMeasurement  = Measurement(:,ones(1,eS)) - Delta*drawSample(obj.sys.v,eS,[],Input,[],obj.time);

      % Kalman gain
      K = predVar*H'*inv(Delta*vVar*Delta' + H*predVar*H');
      filtEnsemble = predEst.RV.values + K*(perturbedMeasurement-zPredMeans);

      % store estimate and new Input and normalize weights
      filtEstimate.Input = Input;
      filtEstimate.RV = nefEmpiricalRV(filtEnsemble,predEst.RV.weights);
    end % function timeMeasurementUpdate
  end % methods
end % the class
