classdef nefSDD2 < nefEstimator
  %file @nefSDD2/nefSDD2.m (Square-Root Divided Difference Estimator 2nd Order)
  % nefSDD2 Properties:
  %   nefDD2_optsDef - definition of the nefSDD2 options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefSDD2 Methods:
  %   nefSDD2 - class constructor (For list of options type 'help nefSDD2/nefSDD2')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   SDD2TimeUpdate - time update of the SDD2
  %   SDD2MeasurementUpdate - measurement update of the SDD2
  %   SD1RTSSmoothUpdate - same as the SD1RTSSmoothUpdate
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
    nefSDD2_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
                      2,'scalingParameterType:constant','parameterValue',sqrt(3),'constant scaling parameter', @(x) isfloat(x) && length(x)==1 && x>1;
                      2,'scalingParameterType:adaptive','parameterDomain',[1.1 2.1 0.5],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && x(1)>1 && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
                      2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE','MMPE'}));};
  end % constant properties
  
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
    %hpar;
    isRConst;
    isQConst;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSDD2(system,varargin)
      % NEFSDD2 Creates NEFSDD2 object.
      %
      %   OBJ = NEFSDD2(system,varargin) creates a NEFSDD2 object OBJ 
      %   representing square-root divided difference filter 2nd order for model SYSTEM.
      %
      %   the following user-defined parameters can be changed 
      %   via the standard Parameter-Value MATLAB mechanism
      %  
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  scalingParameterType        constant         scaling parameter can be             constant, adaptive
      %                                               constant or adaptively chosen   
      %  parameterValue              sqrt(3)          constant scaling parameter           real, greater than 1
      %  parameterDomain             [1.1 2.1 0.5]    adaptive grid search of scaling      positive real, Kmin<Kmax
      %                                               parameter K within domain
      %                                               [Kmin, Kmax, Kstep] for 
      %                                               adaptive setting
      %  parameterCriterion          MLE              selection of criterion for           MLE, MMPE                       
      %                                               adaptive search
      %
      %  See also ESTIMATE.
      
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFSDD2';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for i = 1:size(nefSDD2.nefSDD2_optsDef,1)
        p.addParamValue(nefSDD2.nefSDD2_optsDef{i,3},[],nefSDD2.nefSDD2_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefSDD2.nefSDD2_optsDef,p.Results);
      
      % scaling parameter check
      obj = nefDD2.ScalingParameterCheck(obj);

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
    end % NEFSDD2 constructor
    
    function disp(obj)
      % DISP Displays nefSDD2 object parameters and options.
      
      fprintf('A nefSDD2 object with parameters\n')
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
        [predMean, predVar] = nefSDD2.SDD2TimeUpdate(filtMean,filtVar.S,obj.sys.f,Input,wMean,obj.wSVar,Time,hpar);
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
      if ((obj.isRConst == 1) && (isempty(obj.vSVar))) || (obj.isRConst == 0)
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
          hpar = nefSDD2.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,aux_domain(1:3),aux_criterion);                         
          setOption(obj,'parameterDomain',[aux_domain(1:3), hpar]);
        end
        [filtMean,filtVar] = nefSDD2.SDD2MeasurementUpdate(predMean,predVar.S,obj.sys.h,Input,vMean,obj.vSVar,Time,Measurement,hpar);
      end

      % returning nefGaussianRV with symmetrized variance
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
      % return also by-products if required
    end % function measurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SMOOTHUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Divided Difference RTS Smoother 2nd order is the same as the 1st order
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
    % included in nefDD2

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
          % differences of first order, Szx2, Szv2  containing divided
          % differences of second order, and auxiliary sums for mean
          % computation
          xd0 =evaluate(h,predMean,Input,vMean,Time);
          xd =evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar(i),predMean),Input,vMean,Time);
          xd_ =evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar(i),predMean),Input,vMean,Time);
          Szx1(:,idx) = (xd - xd_)/(2*hpar(i));
          Szx2(:,idx) = (xd + bsxfun(@minus,xd_,2*xd0))*sqrt(hpar(i)^2-1)/(2*hpar(i)^2);
          sumz = sum(xd + xd_,2);
          if h.isAdditive
            Szv1 = Sr;
            Szv2 = zeros(nz);
          else
            xd = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            xd_ = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            Szv1 = (xd - xd_)/(2*hpar(i));
            Szv2 = (xd + xd_ - 2*xd0)*sqrt(hpar(i)^2-1)/(2*hpar(i)^2);
            sumv =  sum(xd + xd_,2);
          end

          % measurement prediction mean and prediction error
          if additive
            zPredMean = ((hpar(i)^2-nx)/hpar(i)^2)*xd0+1/(2*hpar(i)^2)*sumz; % zp
          else
            zPredMean = ((hpar(i)^2-nx-nv)/hpar(i)^2)*xd0+1/(2*hpar(i)^2)*(sumz+sumv); % zp
          end
          dz = Measurement-zPredMean;

          % measurement prediction variance        
          zPredVar = [Szx1 Szv1 Szx2 Szv2]*[Szx1 Szv1 Szx2 Szv2]'; % Pzzp

          %psi(i) = mvnpdf(Measurement,zPredMean,zPredVar);        
          exponent = -0.5*(dz'/zPredVar)*dz;
          psi(i)  = 1/(2*pi)^(nz/2)/sqrt(det(zPredVar))*exp(exponent);
        end % for
        hpar_act = hpar(logical(psi==max(psi)));
      elseif strcmpi(kcriterion,'MMPE')
        % hpar providing the MINIMUM MEASUREMENT PREDICTION ERROR is selected
        mpe = zeros(1,length(hpar));      
        for i=1:length(hpar)        
          xd0 =evaluate(h,predMean,Input,vMean,Time);
          xd =evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar(i),predMean),Input,vMean,Time);
          xd_ =evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar(i),predMean),Input,vMean,Time);
          Szx1(:,idx) = (xd - xd_)/(2*hpar(i));
          Szx2(:,idx) = (xd + bsxfun(@minus,xd_,2*xd0))*sqrt(hpar(i)^2-1)/(2*hpar(i)^2);
          sumz = sum(xd + xd_,2);
          if h.isAdditive
          else
            xd = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            xd_ = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar(i),vMean),Time(:,ones(1,nv)));
            Szv1 = (xd - xd_)/(2*hpar(i));
            Szv2 = (xd + xd_ - 2*xd0)*sqrt(hpar(i)^2-1)/(2*hpar(i)^2);
            sumv =  sum(xd + xd_,2);
          end

          % measurement prediction mean and prediction error
          if h.isAdditive
            zPredMean = ((hpar(i)^2-nx)/hpar(i)^2)*xd0+1/(2*hpar(i)^2)*sumz; % zp
          else
            zPredMean = ((hpar(i)^2-nx-nv)/hpar(i)^2)*xd0+1/(2*hpar(i)^2)*(sumz+sumv); % zp
          end

          % measurement prediction error
          mpe(i)  = (Measurement-zPredMean)'*(Measurement-zPredMean);        
        end % for
        hpar_act = hpar(logical(mpe==min(mpe)));      
      end % if
    
      % for the case of several values of hpar with same likelihood
      hpar_act = hpar_act(1);
    end %function scaling parameter adaptive choice

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SDD2TimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = SDD2TimeUpdate(filtMean,filtSVar,f,Input,wMean,Sq,Time,hpar)
    % SDD2TIMEUPDATE Implements time update of the SDD2.

      nx = size(filtMean,1);
      nw = size(wMean,1);

      Sf = filtSVar;

      % Auxiliary matrices Sxx1, Sxw1  containing divided
      % differences of first order, Sxx2, Sxw2  containing divided
      % differences of second order, and auxiliary sums for mean
      % computation
        xd0 = evaluate(f,filtMean,Input,wMean,Time);
        xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean),Input,wMean,Time);
        xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean),Input,wMean,Time);
        Sxx1 = (xd - xd_)/(2*hpar);
        Sxx2 = (xd + bsxfun(@minus,xd_ ,2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
        sumx = sum(xd + xd_,2);
        if f.isAdditive
          Sxw1 = Sq;
          Sxw2 = zeros(nx);
        else
          xd = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,Sq*hpar,wMean),Time(:,ones(1,nw)));
          xd_ = evaluate(f,filtMean(:,ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,-Sq*hpar,wMean),Time(:,ones(1,nw)));
          Sxw1 = (xd - xd_)/(2*hpar);
          Sxw2 = (xd + bsxfun(@minus,xd_,2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
          sumw = sum(xd + xd_,2);
        end


      % evaluating prediction mean
      if f.isAdditive
        predMean = ((hpar^2-nx)/hpar^2)*xd0+1/(2*hpar^2)*sumx;
      else
        predMean = ((hpar^2-nx-nw)/hpar^2)*xd0+1/(2*hpar^2)*(sumx+sumw);
      end

      % evaluating prediction variance
      Sp = nefDD1.triag([Sxx1 Sxw1 Sxx2 Sxw2]);      
      predVar = nefCholFactorFunction(Sp,'matrixType','fact');
      
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SDD2MeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = SDD2MeasurementUpdate(predMean,predSVar,h,Input,vMean,Sr,Time,Measurement,hpar,additive)
    % SDD2MEASUREMENTUPDATE Implements measurement update of the SDD2.
      
      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);

      Sp = predSVar;

      % Auxiliary matrices Szx1, Szv1 containing divided
      % differences of first order, Szx2, Szv2  containing divided
      % differences of second order, and auxiliary sums for mean
      % computation
      xd0 =evaluate(h,predMean,Input,vMean,Time);
      xd =evaluateOverAllStates(h,bsxfun(@plus,Sp*hpar,predMean),Input,vMean,Time);
      xd_ =evaluateOverAllStates(h,bsxfun(@plus,-Sp*hpar,predMean),Input,vMean,Time);
      Szx1 = (xd - xd_)/(2*hpar);
      Szx2 = (xd + bsxfun(@minus,xd_,2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
      sumz = sum(xd + xd_,2);
      if h.isAdditive
        Szv1 = Sr;
        Szv2 = zeros(nz);
      else
        xd = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,Sr*hpar,vMean),Time(:,ones(1,nv)));
        xd_ = evaluate(h,predMean(:,ones(1,nv)),Input(:,ones(1,nv)),bsxfun(@plus,-Sr*hpar,vMean),Time(:,ones(1,nv)));
        Szv1 = (xd - xd_)/(2*hpar);
        Szv2 = (xd + bsxfun(@minus,xd_, 2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
        sumv =  sum(xd + xd_,2);
      end

      % measurement prediction mean and prediction error
      if h.isAdditive
        zPredMean = ((hpar^2-nx)/hpar^2)*xd0+1/(2*hpar^2)*sumz; % zp
      else
        zPredMean = ((hpar^2-nx-nv)/hpar^2)*xd0+1/(2*hpar^2)*(sumz+sumv); % zp
      end
      dz = Measurement-zPredMean;

      % measurement prediction variance
      Szp = [Szx1 Szv1 Szx2 Szv2];      

      % cross-covariance predictive matrix of x and z
      Pxzp = Sp*Szx1';

      % evaluating filtering mean and variance
      K = Pxzp/Szp'/Szp;
      filtMean = predMean + K*dz;
      Sf = nefDD1.triag([Sp-K*Szx1, K*Szv1, K*Szx2, K*Szv2]);
      filtVar = nefCholFactorFunction(Sf,'matrixType','fact');
      
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
    % D2RTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % same as for SDD1, included in nefSDD1
  end % methods
end % the class
