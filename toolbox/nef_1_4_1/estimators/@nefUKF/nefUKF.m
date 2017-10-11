classdef nefUKF < nefEstimator
  %file @nefUKF/nefUKF.m (Unscented Kalman Estimator)
  % nefUKF Properties:
  %   nefUKF_optsDef - definition of the nefUKF options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefUKF Methods:
  %   nefUKF - class constructor (For list of options type 'help nefUKF/nefUKF')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterCheck - check of chosen scaling parameter and warning spec.
  %   msp - computation of set of deterministically chosen weighted points
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   UKFTimeUpdate - time update of the UKF
  %   UKFMeasurementUpdate - measurement update of the UKF
  %   URTSSmoothUpdate - smoothing based on the UKF and Rauch-Tung-Striebel smoother
  %
  % References:
  % S. J. Julier, J. K. Uhlmann and H. F. Durrant-White (2000):
  % A new method for the nonlinear transformation of means and covariances in filters and estimators.
  % IEEE Trans. On AC 45(3), 477.482.
  %
  % M. Simandl and J. Dunik (2009):
  % Derivative-Free Estimation Methods: New Results and Performance Analysis
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
    nefUKF_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
                      2,'scalingParameterType:constant','parameterValue',@(x) max(3-x.sys.dimState,0),'constant scaling parameter', @(x) isfloat(x) && length(x)==1;
                      2,'scalingParameterType:adaptive','parameterDomain',[0 4 1],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
                      2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE','MMPE'}));
                      2,'','sqrt','chol','matrix square-root used',@(x) any(strcmpi(x,{'chol','svd'}))};
  end % constant properties

  properties (SetAccess = 'protected') % protected properties
    matrixDecompositionFunction;
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefUKF(system,varargin)
      % NEFUKF Creates NEFUKF object.
      %
      %   OBJ = NEFUKF(system,varargin) creates a NEFUKF object OBJ
      %   representing unscented Kalman filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  scalingParameterType        constant         scaling parameter can be             constant, adaptive
      %                                               constant or adaptively chosen
      %  parameterValue              3-dim(x),0       constant scaling parameter           real
      %  parameterDomain             [0 4 1]          adaptive grid search of scaling      nonnegative real, Kmin<Kmax
      %                                               parameter K within domain
      %                                               [Kmin, Kmax, Kstep] for
      %                                               adaptive setting
      %  parameterCriterion          MLE              selection of                         MLE, MMPE
      %                                               criterion for adaptive search
      %  sqrt                        chol             used matrix square-root              chol, ud, svd 
      %
      %
      %  See also ESTIMATE.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFUKF';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for  i = 1:size(nefUKF.nefUKF_optsDef,1)
        p.addParamValue(nefUKF.nefUKF_optsDef{i,3},[],nefUKF.nefUKF_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefUKF.nefUKF_optsDef,p.Results);
      % setting matrix Decomposition function according to the options
      switch getOption(obj,'sqrt')
        case 'chol',
          obj.matrixDecompositionFunction = @(x) nefCholFactorFunction.CholFactor(x);
        case 'svd',
          obj.matrixDecompositionFunction = @(x) nefUKF.SVDDecomp(x);
      end;
      % scaling parameter check
      obj = nefUKF.ScalingParameterCheck(obj);

      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end

      % Issue warnings
      if obj.sys.h.isLinear
        disp('**************')
        disp('nefUKF WARNING: Function in measurement equation is linear, hence the Kalman measurement update will be used instead')
        disp('**************')
      end
      if obj.sys.f.isLinear
        disp('**************')
        disp('nefUKF WARNING: Function in state equation is linear, hence the Kalman time update will be used instead')
        disp('**************')
      end
      obj.estQueue(end+1) = data;
    end % NEFUKF constructor

    function disp(obj)
      % DISP Displays nefUKF object parameters and options.

      fprintf('A nefUKF object with parameters\n')
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
        F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
        Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);
        [predMean,predVar] = nefKalman.KalmanTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,F,Gamma);
        %disp('linear f')
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen - value from the last filtering step
          aux_domain = getOption(obj,'parameterDomain');
          kappa = aux_domain(4);
        end
        [predMean, predVar] = nefUKF.UKFTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,kappa,obj.matrixDecompositionFunction);
      end

      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);

      % is timeUpdate is used for smoothing
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
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);

      % evaluating filtering mean and variance
      % - if measurement eq. is linear (with additive noise), then the predictive part of the KF is used
      if obj.sys.h.isLinear
        H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
        Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);
        [filtMean,filtVar] = nefKalman.KalmanMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,H,Delta);
        %disp('linear h')
      else
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen
          aux_domain = getOption(obj,'parameterDomain');
          aux_criterion = getOption(obj,'parameterCriterion');

          kappa = nefUKF.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,aux_domain(1:3),aux_criterion,obj.matrixDecompositionFunction);

          setOption(obj,'parameterDomain',[aux_domain(1:3), kappa]);
        end
        [filtMean,filtVar] = nefUKF.UKFMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,kappa,obj.matrixDecompositionFunction);
      end;

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
        F = predEstimate.auxData.F;
        [smoothM,smoothV] = nefKalman.RTSKalmanSmoothUpdate(initM,initV,filtM,filtV,predM,predV,F);
        %disp('linear f')
      else
        %kappa = getOption(obj,'scalingParameter');
        if strcmp(getOption(obj,'scalingParameterType'),'constant')
          % if kappa is constant during the experiment
          kappa = getOption(obj,'parameterValue');
        else
          % kappa is adaptively chosen - stored from the prediction
          kappa = predEstimate.auxData.kappa;
        end
        [smoothM,smoothV] = nefUKF.URTSSmoothUpdate(initM,initV,filtM,filtV,predM,predV,wMean,wVar,obj.sys.f,Time,Input,kappa,obj.matrixDecompositionFunction);
      end

      % returning nefGaussianRV with symmetrized variance
      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate

  end % methods
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SVDDecomp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function d = SVDDecomp(P)
      [u,s,v]=svd(P);
      d = u*s^0.5;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ScalingParameterCheck
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ScalingParameterCheck(obj)
    % SCALINGPARAMETERCHECK Checks the chosen scaling parameter against
    % the system description.
      if strcmp(getOption(obj,'scalingParameterType'),'constant')
        % check for negativity if kappa is constant
        if getOption(obj,'parameterValue')<0
          disp('**************')
          disp('nefUKF WARNING: This choice of scaling parameter (kappa) can cause the loss of positive definitness of state estimate covariance matrix.')
          disp('The scaling parameter kappa should be non-negative. Details can be found in documentation.')
          disp('**************')
        end
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
          disp('nefUKF WARNING: Function in measurement equation is linear and independent of the scaling parameter. Its optimisation is useless.')
          %disp('Better choice of the scaling parameter is "constant" (instead of adaptive).')
          disp('**************')
        end
      end
    end % function ScalingParameterCheck

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SigmaPointComputation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [chi,w] = msp(x, P, kappa,decompositionFunction)
      %   MSP Method msp for the calculation of sigma-points, called from methods for prediction, filtration and smoothing.
      %   For standard (non-square-root) estimation algorithms.
      %
      %   function [ch,w]=msp(x,P,kappa)
      %   input parameters                   - x - mean value,
      %                                      - P - covariance matrix,
      %                                      - kappa - parameter for distribution of sigma-points,
      %   output parameters                  - chi - calculated point
      %                                      - w -  weight of the point
      nx = size(x,1);
      decomp = sqrt(nx+kappa)*decompositionFunction(P);
      chi = [x (x(:,ones(1,2*nx)) + [decomp, -decomp])];
      w = [kappa 0.5*ones(1,2*nx)]/(nx+kappa);
    end % function msp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scaling Parameter Adaptive Choice
    % - MLE - maximum likelihood choice
    % - MMPE - minimum measurement prediciton error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function kappa_act = ScalingParameterAdaptiveChoice(predMean,predVar,h,Input,vMean,R,Time,Measurement,kdomain,kcriterion,matrixDecompositionFunction)
    % SCALINGPARAMETERADAPTIVECHOICE Performs adaptive choice of the
    % scaling parameter within predefined domain in filtering step.

    % number of sigma points and their computation
    nx = size(predMean,1);
    nv = size(vMean,1);
    if h.isAdditive
      Jns = 2*(nx) + 1;
    else
      Jns = 2*(nx+nv) + 1;
    end
    nz = size(Measurement,1);

    % grid points
    kappa = kdomain(1):kdomain(3):kdomain(2);

    if strcmpi(kcriterion,'MLE')
      % kappa with the MAXIMUM LIKELIHOOD is selected
      psi = zeros(1,length(kappa));
      for i=1:length(kappa)
        % predictive sigma-points
        dzetap = zeros(nz,Jns);
        if h.isAdditive
          [chip,w] = nefUKF.msp(predMean,predVar,kappa(i),matrixDecompositionFunction);
          % points for computation prediction of measurement
          dzetap = evaluateOverAllStates(h,chip,Input,vMean,Time);
        else
          [chip,w] = nefUKF.msp([predMean; vMean],[predVar, zeros(nx,nv); zeros(nv,nx), R],kappa(i),matrixDecompositionFunction);
          % points for computation prediction of measurement
          dzetap = evaluate(h,chip(1:nx,:),Input(:,ones(1,Jns)),chip(nx+1:nx+nv,:),Time(:,ones(1,Jns)));
        end

        % measurement prediction mean
        zPredMean = dzetap*w'; % zp

        % measurement prediction variance
        dzetap_ = bsxfun(@minus,dzetap,zPredMean);
        dzetap_w = dzetap_.*w(ones(1,nz),:);
        if h.isAdditive
          zPredVar = dzetap_w*dzetap_'+R; % Pzzp
        else
          zPredVar = dzetap_w*dzetap_'; % Pzzp
        end
        zPredVar = (zPredVar+zPredVar')/2;

        %psi(i) = mvnpdf(Measurement,zPredMean,zPredVar);
        %exponent = -0.5*(Measurement-zPredMean)'*inv(zPredVar)*(Measurement-zPredMean);
        exponent = -0.5*((Measurement-zPredMean)'/zPredVar)*(Measurement-zPredMean);
        psi(i)  = 1/(2*pi)^(nz/2)/sqrt(det(zPredVar))*exp(exponent);
      end % for
      kappa_act = kappa(logical(psi==max(psi)));
    elseif strcmpi(kcriterion,'MMPE')
      % kappa providing the MINIMUM MEASUREMENT PREDICTION ERROR is selected
      mpe = zeros(1,length(kappa));
      for i=1:length(kappa)
        % predictive sigma-points
        if h.isAdditive
          [chip,w] = nefUKF.msp(predMean,predVar,kappa(i),matrixDecompositionFunction);
          % points for computation prediction of measurement
          dzetap = evaluateOverAllStates(h,chip,Input,vMean,Time);
        else
          [chip,w] = nefUKF.msp([predMean; vMean],[predVar, zeros(nx,nv); zeros(nv,nx), R],kappa(i),matrixDecompositionFunction);
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

    % for the case of several values of kappa with same value of criterion
    kappa_act = kappa_act(1);
    end %function scaling parameter adaptive choice

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UKFTimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = UKFTimeUpdate(filtMean,filtVar,f,Input,wMean,Q,Time,kappa,decompositionFunction)
    % UKFTIMEUPDATE Implements time update of the UKF.

      numval = size(filtMean,2);

      %  number of sigma-points
      nx = size(filtMean,1);
      nw = size(wMean,1);
      if f.isAdditive
        Jns = 2*nx+ 1;
      else
        Jns = 2*(nx+nw)+ 1;
      end

      predMean = zeros(nx,numval);
      predVar = zeros(nx,nx,numval);

      for in = 1:numval
        % filtering sigma-points
        if f.isAdditive
          [chi,w] = nefUKF.msp(filtMean(:,in),filtVar(:,:,in),kappa,decompositionFunction);
        else
          [chi,w] = nefUKF.msp([filtMean(:,in); wMean(:,in)],[filtVar(:,:,in), zeros(nx,nw); zeros(nw,nx), Q(:,:,in)],kappa,decompositionFunction);
        end

        % computing predictive sigma-points
        %chip = evaluate(f,chi(1:nx,:),repmat(Input,1,Jns),chi(nx+1:nx+nw,:),repmat(Time,1,Jns));
        if f.isAdditive
          chip = evaluateOverAllStates(f,chi,Input,wMean(:,in),Time);
        else
          chip = evaluate(f,chi(1:nx,:),Input(:,ones(1,Jns)),chi(nx+1:nx+nw,:),Time(:,ones(1,Jns)));
        end

        % evaluating prediction mean
        predMean(:,in) = chip*w';

        % evaluating prediction variance
        chip_ = bsxfun(@minus,chip,predMean);
        chip_w = chip_.*w(ones(1,nx),:);
        if f.isAdditive
          predVar(:,:,in) = chip_w*chip_' + Q(:,:,in);
        else
          predVar(:,:,in) = chip_w*chip_';
        end
        predVar(:,:,in) = (predVar(:,:,in)+predVar(:,:,in)')/2;
      end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UKFMeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = UKFMeasurementUpdate(predMean,predVar,h,Input,vMean,R,Time,Measurement,kappa,decompositionFunction)
      % UKFMEASUREMENTUPDATE Implements measurement update of the UKF.

      numval = size(predMean,2);

      % number of sigma points
      nx = size(predMean,1);
      nv = size(vMean,1);
      if h.isAdditive
        Jns = 2*(nx) + 1;
      else
        Jns = 2*(nx+nv) + 1;
      end
      nz = size(Measurement,1);

      filtMean = zeros(nx,numval);
      zPredMean = zeros(nz,numval);
      zPredVar = zeros(nz,nz,numval);
      filtVar = zeros(nx,nx,numval);

      for in = 1:numval
        % predictive sigma-points
        if h.isAdditive
          [chip,w] = nefUKF.msp(predMean(:,in),predVar(:,:,in),kappa,decompositionFunction);
          % points for computation prediction of measurement
          dzetap = evaluateOverAllStates(h,chip,Input,vMean(:,in),Time);
        else
          [chip,w] = nefUKF.msp([predMean(:,in); vMean(:,in)],[predVar(:,:,in), zeros(nx,nv); zeros(nv,nx), R(:,:,in)],kappa,decompositionFunction);
          % points for computation prediction of measurement
          dzetap = evaluate(h,chip(1:nx,:),Input(:,ones(1,Jns)),chip(nx+1:nx+nv,:),Time(:,ones(1,Jns)));
        end


        % measurement prediction mean and prediction error
        zPredMean(:,in) = dzetap*w'; % zp
        dz = Measurement(:,in)-zPredMean(:,in);

        % measurement prediction variance
        dzetap_ = bsxfun(@minus,dzetap,zPredMean);
        dzetap_w = dzetap_.*w(ones(1,nz),:);
        if h.isAdditive
          zPredVar(:,:,in) = dzetap_w*dzetap_'+R(:,:,in); % Pzzp
        else
          zPredVar(:,:,in) = dzetap_w*dzetap_'; % Pzzp
        end
        zPredVar(:,:,in) = (zPredVar(:,:,in)+zPredVar(:,:,in)')/2;

        % cross-covariance predictive matrix of x and z
        chip_ = bsxfun(@minus,chip(1:nx,:),predMean);
        Pxzp = chip_*dzetap_w';

        % evaluating filtering mean and variance
        K = Pxzp/zPredVar(:,:,in);
        filtMean(:,in) = predMean(:,in) + K*dz;
        filtVar(:,:,in) = predVar(:,:,in) - K * zPredVar(:,:,in) * K';
        filtVar(:,:,in) = (filtVar(:,:,in)+filtVar(:,:,in)')/2;
      end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        % used for a particle and GSM filter
        invzPredVar = zeros(nz,nz,numval);
        detzPredVar = zeros(1,numval);
        for in = 1:numval
          invzPredVar(:,:,in) = inv(zPredVar(:,:,in));
          detzPredVar(:,in) = det(zPredVar(:,:,in));
        end
        varargout{3} = zPredMean;
        varargout{4} = invzPredVar;
        varargout{5} = detzPredVar;
      end
    end %function measurement update

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % URTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean,smoothVar] = URTSSmoothUpdate(initMean,initVar,filtMean,filtVar,predMean,predVar,wMean,Q,f,Time,Input,kappa,decompositionFunction)
    % URTSSMOOTHUPDATE Implements smoothing step based on the UKF and Rauch-Tung-Striebel smoother.

      nx = size(filtMean,1);
      nw = size(wMean,1);
      if f.isAdditive
        Jns = 2*(nx) + 1;
      else
        Jns = 2*(nx+nw) + 1;
      end

      % filtering sigma-points
      if f.isAdditive
        [chi_v,w] = nefUKF.msp(filtMean,filtVar,kappa,decompositionFunction);
        % computing predictive sigma-points
        chip_v = evaluateOverAllStates(f,chi_v,Input,wMean,Time);
      else
        [chi_v,w] = nefUKF.msp([filtMean; wMean],[filtVar, zeros(nx,nw); zeros(nw,nx), Q],kappa,decompositionFunction);
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
      chip_v_w = chip_v_.*w(ones(1,nx),:);

      % it should hold chip_v_*chip_v_w' =  predVar

      Pxx = chi_v_(1:nx,:)*chip_v_w';

      % smoother gain
      Kv = Pxx/predVar;

      % smoothed mean and covariance matrix
      smoothMean = filtMean - Kv * (predMean - initMean);
      smoothVar = filtVar - Kv * (predVar - initVar) * Kv';
      smoothVar = (smoothVar+smoothVar')/2;
    end

  end % methods
end % the class
