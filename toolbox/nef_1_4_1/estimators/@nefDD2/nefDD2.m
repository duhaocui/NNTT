classdef nefDD2 < nefEstimator
  %file @nefDD2/nefDD2.m (Divided Difference Estimator 2nd Order)
  % nefDD2 Properties:
  %   nefDD2_optsDef - definition of the nefDD2 options (For list of options type 'help nefEstimator/nefEstimator_optsDef')
  % nefDD2 Methods:
  %   nefDD2 - class constructor (For list of options type 'help nefDD2/nefDD2')
  %   disp - display class info
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   ScalingParameterCheck - check of chosen scaling parameter and warning spec.
  %   ScalingParameterAdaptiveChoice - adaptive computation of scaling parameter
  %   DD2TimeUpdate - time update of the DD2
  %   DD2MeasurementUpdate - measurement update of the DD2
  %   D2RTSSmoothUpdate - same as the D1RTSSmoothUpdate
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
    nefDD2_optsDef = {2,'','scalingParameterType','constant','scaling parameter (constant or adaptive setting)',@(x) any(strcmpi(x,{'constant','adaptive'}));
                      2,'scalingParameterType:constant','parameterValue',sqrt(3),'constant scaling parameter', @(x) isfloat(x) && length(x)==1 && x>1;
                      2,'scalingParameterType:adaptive','parameterDomain',[1.1 2.1 0.5],'domain for adaptive setting of parameter', @(x) all(isfloat(x)) && x(1)>1 && all(abs(x)==x) && size(x,1)==1 && size(x,2)==3 && x(1)<x(2);
                      2,'scalingParameterType:adaptive','parameterCriterion','MLE','selection of criterion for adaptive setting', @(x) any(strcmpi(x,{'MLE','MMPE'}));};
  end % constant properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefDD2(system,varargin)
      % NEFDD2 Creates NEFDD2 object.
      %
      %   OBJ = NEFDD2(system,varargin) creates a NEFDD2 object OBJ 
      %   representing divided difference filter 2nd order for model SYSTEM.
      %
      %   the following user-defined parameters can be changed 
      %   via the standard Parameter-Value MATLAB mechanism
      %  
      %  PARAMETER                   DEFAULT VAL      DESCRIPTION                          VALID VALUES
      %  =========================================================================================================================
      %  scalingParameterType        constant         scaling parameter can be             constant, adaptive
      %                                               constant or adaptively chosen   
      %  parameterValue              sqrt(3)          constant scaling parameter           real, greater than 1
      %  parameterDomain             [1.1 3.1 0.5]    adaptive grid search of scaling      positive real, Kmin<Kmax
      %                                               parameter K within domain 
      %                                               [Kmin, Kmax, Kstep] for 
      %                                               adaptive setting
      %  parameterCriterion          MLE              selection of criterion for           MLE, MMPE                       
      %                                               adaptive search
      %
      %  See also ESTIMATE.
      
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFDD2';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for i = 1:size(nefDD2.nefDD2_optsDef,1)
        p.addParamValue(nefDD2.nefDD2_optsDef{i,3},[],nefDD2.nefDD2_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefDD2.nefDD2_optsDef,p.Results);
      
      % scaling parameter check
      obj = nefDD2.ScalingParameterCheck(obj);

      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end
      obj.estQueue(end+1) = data;
    end % NEFDD2 constructor
    
    function disp(obj)
      fprintf('A nefDD2 object with parameters\n')
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
        [predMean, predVar] = nefDD2.DD2TimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,hpar);
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
          hpar = nefDD2.ScalingParameterAdaptiveChoice(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,aux_domain(1:3),aux_criterion);                         
          setOption(obj,'parameterDomain',[aux_domain(1:3), hpar]);
        end
        [filtMean,filtVar] = nefDD2.DD2MeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,hpar);
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
          disp('nef(S)DD2 WARNING: Function in measurement equation is linear and independent of the scaling parameter. Its optimisation is useless.')
          %disp('Better choice of the scaling parameter is "constant" (instead of adaptive).')
          disp('**************')
        end                
      end
    end % function ScalingParameterCheck


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MatrixTriangularization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % included in nefDD1
    
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
        if h.isAdditive
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
        % Auxiliary sums for mean computation
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
    % DD2TimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = DD2TimeUpdate(filtMean,filtVar,f,Input,wMean,Q,Time,hpar)
    % DD2TIMEUPDATE Implements time update of the DD2.
      
      numval = size(filtMean,2);

      nx = size(filtMean,1);
      nw = size(wMean,1);

      predMean = zeros(nx,numval);
      predVar = zeros(nx,nx,numval);

      for in = 1:numval
        Sf = nefCholFactorFunction.CholFactor(filtVar(:,:,in));
        Sq = nefCholFactorFunction.CholFactor(Q(:,:,in));

        % Auxiliary matrices Sxx1, Sxw1  containing divided
        % differences of first order, Sxx2, Sxw2  containing divided
        % differences of second order, and auxiliary sums for mean
        % computation
        xd0 = evaluate(f,filtMean(:,in),Input,wMean(:,in),Time);
        xd = evaluateOverAllStates(f,bsxfun(@plus,Sf*hpar,filtMean(:,in)),Input,wMean(:,in),Time);
        xd_ = evaluateOverAllStates(f,bsxfun(@plus,-Sf*hpar,filtMean(:,in)),Input,wMean(:,in),Time);
        Sxx1 = (xd - xd_)/(2*hpar);
        Sxx2 = (xd + bsxfun(@minus,xd_ ,2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
        sumx = sum(xd + xd_,2);
        if f.isAdditive
          Sxw1 = Sq;
          Sxw2 = zeros(nx);
        else
          xd = evaluate(f,filtMean(:,in*ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,Sq*hpar,wMean(:,in)),Time(:,ones(1,nw)));
          xd_ = evaluate(f,filtMean(:,in*ones(1,nw)),Input(:,ones(1,nw)),bsxfun(@plus,-Sq*hpar,wMean(:,in)),Time(:,ones(1,nw)));
          Sxw1 = (xd - xd_)/(2*hpar);
          Sxw2 = (xd + bsxfun(@minus,xd_,2*xd0))*sqrt(hpar^2-1)/(2*hpar^2);
          sumw = sum(xd + xd_,2);
        end

        % evaluating prediction mean
        if f.isAdditive
          predMean(:,in) = ((hpar^2-nx)/hpar^2)*xd0+1/(2*hpar^2)*sumx;
        else
          predMean(:,in) = ((hpar^2-nx-nw)/hpar^2)*xd0+1/(2*hpar^2)*(sumx+sumw);
        end

        % evaluating prediction variance
        %Sp = nefDD1.triag([Sxx1 Sxw1 Sxx2 Sxw2]);
        Sp = [Sxx1 Sxw1 Sxx2 Sxw2];
        predVar(:,:,in) = Sp*Sp';
      end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DD2MeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = DD2MeasurementUpdate(predMean,predVar,h,Input,vMean,R,Time,Measurement,hpar)
    % DD2MEASUREMENTUPDATE Implements measurement update of the DD2.  
      
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
          zPredMean(:,in) = ((hpar^2-nx)/hpar^2)*xd0+1/(2*hpar^2)*sumz; % zp
        else
          zPredMean(:,in) = ((hpar^2-nx-nv)/hpar^2)*xd0+1/(2*hpar^2)*(sumz+sumv); % zp
        end
        dz = Measurement(:,in)-zPredMean(:,in);

        % measurement prediction variance
        Szp = [Szx1 Szv1 Szx2 Szv2];
        zPredVar(:,:,in) = Szp*Szp'; % Pzzp

        % cross-covariance predictive matrix of x and z
        Pxzp = Sp*Szx1';

        % evaluating filtering mean and variance
        K = Pxzp/zPredVar(:,:,in);
        filtMean(:,in) = predMean(:,in) + K*dz;
        %Sf = nefDD1.triag([Sp-K*Szx1, K*Szv1, K*Szx2, K*Szv2]);
        Sf = [Sp-K*Szx1, K*Szv1, K*Szx2, K*Szv2];
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
    % D2RTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % same as for DD1, included in nefDD1
  end % methods
end % the class
