classdef nefSIF < nefEstimator
  %file @nefSIF/nefSIF.m Stochastic Integration Filter
  % nefSIF Methods:
  %   nefSIF - class constructor
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   smoothUpdate - smoothing step with respect to system description
  %   SIFTimeUpdate - implements time update of the SIF
  %   SIFMeasurementUpdate - implements measurement update of the SIF
  %   RURTSSmoothUpdate - implements smoothing step based on the SIF and Rauch-Tung-Striebel smoother
  %
  % References:
  % J. Dunik, O. Straka and M. Simandl (2013):
  %   Stochastic Integration Filter
  %   IEEE Trans. Automat. Contr., vol. 58(6), pp. 1561 - 1566
  %   DOI: 10.1109/TAC.2013.2258494

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (Constant)
    % definition of estimator options
    nefSIF_optsDef = {2,'','SIorder',3,'stochastic integration rule order', @(x) isscalar(x) && abs(x) ==x && (mod(x,1) ==0);
        2,'','Nmin',1,'minimum number of iterations', @(x) isscalar(x) && abs(x) ==x;
        2,'','Nmax',3,'maximum number of iterations', @(x) isscalar(x) && abs(x) ==x;
        2,'','Eps',1e-6,'epsilon stop condition', @(x) isscalar(x) && abs(x) ==x;};
  end % constant properties
  properties (SetAccess = 'protected') % protected properties
    optNmin;
    optNmax;
    optEps;
    optSIorder;
  end % properties
  methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSIF(system,varargin)
      % NEFSIF Creates NEFSIF object.
      %
      %   OBJ = NEFSIF(system,varargin) creates a NEFSIF object OBJ
      %   representing Stochastic Integration Filter for model SYSTEM.
      %
      %   the following user-defined parameters can be changed
      %   via the standard Parameter-Value MATLAB mechanism
      %
      %  PARAMETER   DEFAULT VAL    DESCRIPTION                          VALID VALUES
      %  ===========================================================================================================
      %  Nmin          1           minimum number of iterations         nonnegative real, scalar
      %
      %  Nmax          3           maximum number of iterations         nonnegative real, scalar
      %
      %  Eps          1e-6         epsilon stop condition               nonnegative real, scalar
      %
      %  SIorder      3            order of SIF (3 for RUKF)            nonnegative int, scalar
      %
      %  See also ESTIMATE.
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFSIF';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      for  i = 1:size(nefSIF.nefSIF_optsDef,1)
          p.addParamValue(nefSIF.nefSIF_optsDef{i,3},[],nefSIF.nefSIF_optsDef{i,6});
      end
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      % options processing
      processOptions(obj,nefSIF.nefSIF_optsDef,p.Results);
      obj.optNmin = getOption(obj,'Nmin');
      obj.optNmax = getOption(obj,'Nmax');
      obj.optEps = getOption(obj,'Eps');
      obj.optSIorder = getOption(obj,'SIorder');
      if system.f.isLinear
        disp('linear f -> switching to using Kalman Filter time update')
      end
      if strcmpi(getOption(obj,'taskType'),'fixedLagSmoothing')
        %  smoothers require a complex structure in the queue
        data.predEstimate.RV = obj.x0;
        data.time = -1; % for debugging purposes
        if system.f.isLinear
          disp('linear f -> switching to using Kalman Filter smooth update')
        end
      else
        % predictors and, filters require prediction only in the queue
        data.RV = obj.x0;
      end
      obj.estQueue(end+1) = data;
      if system.h.isLinear
        disp('linear h -> switching to using Kalman Filter measurement update')
      end

    end % NEFSIF constructor

    function disp(obj)
      % DISP Displays nefUKF object parameters and options.

      fprintf('A nefSIF object with parameters\n')
      showOptions(obj)
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time,varargin)
    % TIMEUPDATE Performs time update step with respect to system description.

      p = inputParser;
      p.addParamValue('smoothingPurpose',0,@(x) (x ==0) || (x ==1) );
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
      else
        [predMean, predVar] = nefSIF.SIFTimeUpdate(filtMean,filtVar,obj.sys.f,Input,wMean,wVar,Time,obj.optNmin,obj.optNmax,obj.optEps,obj.optSIorder);
      end

      % returning nefGaussianRV with symmetrized variance
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);

      % is timeUpdate is used for smoothing
      if p.Results.smoothingPurpose
          predEstimate.auxData.Input = Input;
          predEstimate.auxData.Time = Time;
          predEstimate.auxData.wMean = wMean;
          predEstimate.auxData.wVar = wVar;
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
      else
        [filtMean,filtVar] = nefSIF.SIFMeasurementUpdate(predMean,predVar,obj.sys.h,Input,vMean,vVar,Time,Measurement,obj.optNmin,obj.optNmax,obj.optEps,obj.optSIorder);
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
      else
        [smoothM,smoothV] = nefSIF.RURTSSmoothUpdate(initM,initV,filtM,filtV,predM,predV,wMean,wVar,obj.sys.f,Time,Input,obj.optNmin,obj.optNmax,obj.optEps,obj.SIForder);
      end

      % returning nefGaussianRV with symmetrized variance
      smoothEstimate.RV = nefGaussianRV(smoothM,smoothV,'check',0);
    end % function smoothUpdate

  end % methods

  methods (Static)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIFTimeUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean, predVar] = SIFTimeUpdate(filtMean,filtVar,f,Input,wMean,Q,Time,Nmin,Nmax,Eps,SIorder)
    % UKFTIMEUPDATE Implements time update of the SIF.
      numval = size(filtMean,2);
      nx = size(filtMean,1);
      nw = size(wMean,1);
      predMean = zeros(nx,numval);
      predVar = zeros(nx,nx,numval);
      if f.isAdditive
        enx = nx;
      else
        enx = nx+nw;
      end
      for in = 1:numval
        if f.isAdditive
          efMean = filtMean(:,in);
          efVar = filtVar(:,:,in);
        else
          efVar = [filtVar(:,:,in),zeros(nx,nx);zeros(nx,nx), Q(:,:,in)];
          efMean = [filtMean(:,in);wMean(:,in)];
        end
        Sp = nefCholFactorFunction.CholFactor(efVar);
        Ix  = zeros(nx,1);
        Vx  = zeros(nx,1);
        IPx  = zeros(nx);
        VPx  = zeros(nx);
        N = 0; % number of iterations
        while N<Nmin || all([N<Nmax, any([(norm(Vx)> Eps) (norm(VPx)> Eps)])])
          N = N+1;
          [SCRSigmaPoints,w] = nefSIF.stochasticCubatureRulePoints(enx,SIorder);
          Jns = length(w); % # sigmapoints
          xpoints = bsxfun(@plus,Sp*SCRSigmaPoints,efMean);
          fpoints = zeros(nx,Jns);
          if f.isAdditive
            fpoints = evaluateOverAllStates(f,xpoints,Input,wMean(:,in),Time);
          else
            fpoints = evaluate(f,xpoints(1:nx,:),Input(:,ones(1,Jns)),xpoints(nx+1:enx,:),Time(:,ones(1,Jns)));
          end
          SumRx = fpoints*w';
          SumRPx=(w(ones(nx,1),:).*fpoints)*fpoints';
          % update Ix
          Dx = (SumRx-Ix)/N;
          Ix = Ix+Dx;
          Vx = (N-2)*Vx/N+Dx.^2;
          % update IPx
          DPx = (SumRPx-IPx)/N;
          IPx = IPx+DPx;
          VPx = (N-2)*VPx/N+DPx.^2;
        end
        predMean(:,in) = Ix;
        predVar(:,:,in) = IPx - Ix*Ix';
        if f.isAdditive
          predVar(:,:,in) = predVar(:,:,in) + Q(:,:,in);
        end
        predVar(:,:,in) = (predVar(:,:,in)+predVar(:,:,in)')/2;
      end
    end % function timeUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIFMeasurementUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = SIFMeasurementUpdate(predMean,predVar,h,Input,vMean,R,Time,Measurement,Nmin,Nmax,Eps,SIorder)
      % UKFMEASUREMENTUPDATE Implements measurement update of the SIF.
      numval = size(predMean,2);
      nx = size(predMean,1);
      nv = size(vMean,1);
      nz = size(Measurement,1);
      filtMean = zeros(nx,numval);
      filtVar = zeros(nx,nx,numval);
      Pzp = zeros(nz,nz,numval);
      if h.isAdditive
        enx = nx;
      else
        enx = nx+nv;
      end
      for in = 1:numval
        if h.isAdditive
          epMean = predMean(:,in);
          epVar = predVar(:,:,in);
        else
          epVar = [predVar(:,:,in),zeros(nx,nv);zeros(nv,nx), R(:,:,in)];
          epMean = [predMean(:,in);vMean(:,in)];
        end
        Sp = nefCholFactorFunction.CholFactor(epVar);
        Iz  = zeros(nz,1);
        Vz  = zeros(nz,1);
        IPz  = zeros(nz);
        VPz  = zeros(nz);
        IPxz = zeros(nx,nz);
        VPxz  = zeros(nx,nz);
        N = 0; % number of iterations
        while N<Nmin || all([N<Nmax, any([(norm(Vz)> Eps) (norm(VPz)> Eps) (norm(VPxz)> Eps)])])
          N = N+1;
          [SCRSigmaPoints,w] = nefSIF.stochasticCubatureRulePoints(enx,SIorder);
          Jns = length(w); % # sigmapoints
          xpoints = bsxfun(@plus,Sp*SCRSigmaPoints,epMean);
          if h.isAdditive
            hpoints = evaluateOverAllStates(h,xpoints,Input,vMean(:,in),Time);
          else
            hpoints = evaluate(h,xpoints(1:nx,:),Input(:,ones(1,Jns)),xpoints(nx+1:enx,:),Time(:,ones(1,Jns)));
          end
          SumRz = hpoints*w';
          SumRPz = hpoints*(hpoints.*w(ones(nz,1),:))';
          SumRPxz = xpoints(1:nx,:)*(hpoints.*w(ones(nz,1),:))';
          % update Iz
          Dz = (SumRz-Iz)/N;
          Iz = Iz+Dz;
          Vz = (N-2)*Vz/N+Dz.^2;
          % update Pz
          DPz = (SumRPz-IPz)/N;
          IPz = IPz+DPz;
          VPz = (N-2)*VPz/N+DPz.^2;
          % update Pxz
          DPxz = (SumRPxz-IPxz)/N;
          IPxz = IPxz+DPxz;
          VPxz = (N-2)*VPxz/N+DPxz.^2;
        end
        zp = Iz;
        Pzp(:,:,in) = IPz - zp*zp';
        if h.isAdditive
          Pzp(:,:,in) = Pzp(:,:,in) + R(:,:,in);
        end
        Pxzp = IPxz - predMean(:,in)*zp';
        K = Pxzp/Pzp(:,:,in);
        filtMean(:,in) = predMean(:,in)+ K*(Measurement(:,in)-zp);
        filtVar(:,:,in) = predVar(:,:,in) - K*Pzp(:,:,in)*K';
        filtVar(:,:,in) = (filtVar(:,:,in)+filtVar(:,:,in)')/2;
        mieig = min(eig(filtVar(:,:,1)));
        if mieig < 0
          filtVar(:,:,1) =  filtVar(:,:,1) - mieig*eye(nx);
        end
      end
      varargout{1} = filtMean;
      varargout{2} = filtVar;
      if nargout == 5
        % used for a particle and GSM filter
        invzPredVar = zeros(nz,nz,numval);
        detzPredVar = zeros(1,numval);
        for in = 1:numval
          invzPredVar(:,:,in) = inv(Pzp(:,:,in));
          detzPredVar(:,in) = det(Pzp(:,:,in));
        end
        varargout{3} = Pzp(:,:,in);
        varargout{4} = invzPredVar;
        varargout{5} = detzPredVar;
      end
    end %function measurement update

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RURTSSmoothUpdate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [smoothMean, smoothVar] = RURTSSmoothUpdate(initMean,initVar,filtMean,filtVar,predMean,predVar,wMean,Q,f,Time,Input,Nmin,Nmax,epsilon,SIForder)
    % URTSSMOOTHUPDATE Implements smoothing step based on the SIF and Rauch-Tung-Striebel smoother.
      nx = size(filtMean,1);
      enx = nx;
      if f.isAdditive
      else
        enx = nx*2;
      end
      Jns = 1+2*enx;
      if f.isAdditive
        enVar = filtVar;
      else
        enVar = [filtVar,zeros(nx,nx);zeros(nx,nx), Q];
      end
      Sp = nefCholFactorFunction.CholFactor(enVar);
      Ix  = zeros(nx,1);
      IPx  = zeros(nx);
      if f.isAdditive
        efMean = filtMean;
      else
        efMean = [filtMean;wMean];
      end
      N = 0;
      while N<Nmin || all([N<Nmax, any([(norm(Ix)> epsilon) (norm(IPx)> epsilon)])])
        N = N+1;
        [sigmapoints,w] = nefSIF.stochasticCubatureRule(enx,SIForder);
        xpoints = bsxfun(@plus,Sp*sigmapoints,efMean);
        if f.isAdditive
          fpoints = evaluateOverAllStates(f,xpoints,Input,wMean,Time);
        else
          fpoints = evaluate(f,xpoints(1:nx,:),Input(:,ones(1,Jns)),xpoints((nx+1):enx,:),Time(:,ones(1,Jns)));
        end
        %w = [1-enx/rho^2,0.5*ones(1,Jns-1)/rho^2];
        SumRx = fpoints*w';
        SumRPx = bsxfun(@minus,xpoints(1:nx,:),filtMean)*(w(ones(1,nx),:).*fpoints)';
        Dx = (SumRx-Ix)/N;
        Ix = Ix+Dx;
        DPx = (SumRPx-IPx)/N;
        IPx = IPx+DPx;
      end
      Ppred = IPx - bsxfun(@minus,xpoints(1:nx,:),filtMean)*(w(ones(1,nx),:).*Ix(:,ones(1,Jns)))';
      Kv = Ppred/predVar;
      smoothMean = filtMean - Kv * (predMean - initMean);
      smoothVar  = filtVar  - Kv * (predVar  - initVar ) * Kv';
      smoothVar  = (smoothVar+smoothVar')/2;

    end % function timeUpdate

    function [SCRSigmaPoints,weights] = stochasticCubatureRulePoints(nx,order)
      switch order
        case 1
          % TODO
          X = randn(nx,1);
          SCRSigmaPoints = [X,-X];
          weights = [0.5 0.5];
        case 3
          CRSigmaPoints = nefSIF.cubatureRulePoints(nx,order);
          rho = sqrt(chi2rnd(nx+2));
          Q = nefSIF.RandOrthMat(nx);
          SCRSigmaPoints = Q*rho*CRSigmaPoints;
          weights=[1-nx/rho^2,0.5*ones(1,2*nx)/rho^2];
        case 5
          % generating random values
          r=sqrt(chi2rnd(2*nx+7));
          q=betarnd(nx+2,3/2);
          rho=r*sin(asin(q)/2);
          delta=r*cos(asin(q)/2);

          % calculating weights
          c1up=nx+2-delta^2;
          c1do=rho^2*(rho^2-delta^2);
          c2up=nx+2-rho^2;
          c2do=delta^2*(delta^2-rho^2);
          cdo=2*(nx+1)^2*(nx+2);
          c3=(7-nx)*nx^2;
          c4=4*(nx-1)^2;
          coef1=c1up*c3/cdo/c1do;
          coef2=c2up*c3/cdo/c2do;
          coef3=c1up*c4/cdo/c1do;
          coef4=c2up*c4/cdo/c2do;

          weights=[1-nx*(rho^2+delta^2-nx-2)/(rho^2*delta^2),ones(1,2*nx+2)*coef1,ones(1,2*nx+2)*coef2,ones(1,nx*(nx+1))*coef3,ones(1,nx*(nx+1))*coef4];

          %calculating sigma points
          Q = nefSIF.RandOrthMat(nx);
          v=Q*nefSIF.simplex(nx);
          y=zeros(nx,nx*(nx+1)/2);
          cnt=0;
          for j=1:nx+1
            for i=1:j-1
              cnt=cnt+1;
              y(:,cnt)=(v(:,j)+v(:,i))/norm(v(:,j)+v(:,i),2);
            end
          end

          SCRSigmaPoints=[zeros(nx,1),-rho*v,rho*v,-delta*v,+delta*v,-rho*y,rho*y,-delta*y,delta*y];
        otherwise
          disp('We are working hard, but...')
      end
    end

    function V=simplex(m)
      V=zeros(m,m+1);
      for i=1:m
        V(i,i)=sqrt((m+1)*(m-i+1)/m/(m-i+2));
        for j=i+1:m+1
          V(i,j)=-sqrt((m+1)/((m-i+1)*m*(m-i+2)));
        end
      end
    end %simplex

    function [CRSigmaPoints] = cubatureRulePoints(nx,order)
      switch order
        case 3
          CRSigmaPoints = [zeros(nx,1),eye(nx),-eye(nx)];
        otherwise
          disp('We are working hard, but...');
      end
    end

    function [M] = RandOrthMat(n, tol)
      % M = RANDORTHMAT(n)
      % generates a random n x n orthogonal real matrix.
      %
      % M = RANDORTHMAT(n,tol)
      % explicitly specifies a thresh value that measures linear dependence
      % of a newly formed column with the existing columns. Defaults to 1e-6.
      %
      % In this version the generated matrix distribution *is* uniform over the manifold
      % O(n) w.r.t. the induced R^(n^2) Lebesgue measure, at a slight computational
      % overhead (randn + normalization, as opposed to rand ).
      %
      % (c) Ofek Shilon , 2006.
      if nargin ==1
        tol = 1e-6;
      end
      M = zeros(n); % prealloc
      % gram-schmidt on random column vectors
      vi = randn(n,1);
      % the n-dimensional normal distribution has spherical symmetry, which implies
      % that after normalization the drawn vectors would be uniformly distributed on the
      % n-dimensional unit sphere.
      M(:,1) = vi ./ norm(vi);
      for i = 2:n
        nrm = 0;
        while nrm<tol
          vi = randn(n,1);
          vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
          nrm = norm(vi);
        end
        M(:,i) = vi ./ nrm;
      end %i
    end  % RandOrthMat

  end %static methods

end % classdef
