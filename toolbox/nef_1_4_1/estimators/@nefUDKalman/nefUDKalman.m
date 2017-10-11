classdef nefUDKalman < nefEstimator
  %file @nefUDKalman/nefUDKalman.m (Square-Root (Extended) Kalman Filter - based on UD factorization)
  % nefUDKalman Methods:
  %   nefUDKalman - class constructor
  %   INTERNAL METHODS
  %   timeUpdate - time update step with respect to system description
  %   measurementUpdate - measurement update step with respect to system description
  %   UDKalmanTimeUpdate - time update of the square-root (E)KF
  %   UDKalmanMeasurementUpdate - measurement update of the square-root (E)KF
  %
  % References:
  % D. Simon (2006):
  % Optimal State Estimation: Kalman, H infinity and Nonlinear Approaches.
  % John Wiley & Sons, Inc.

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties

    %%%%%%%%%%%%%%%%%%%%%
    % MATRICES
    %%%%%%%%%%%%%%%%%%%%%
    Ur@double = [];
    Dr@double = [];
    Rtilde@double = [];
    invUr@double = [];
    %%%%%%%%%%%%%%%%%%%%%%%
    %PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    isRtildeConst;
  end

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefUDKalman(system,varargin)
      % NEFUDKALMAN Creates NEFUDKALMAN object.
      %
      %   OBJ = NEFUDKALMAN(system,varargin) creates a NEFUDKALMAN object OBJ 
      %   representing UD factorized Kalman filter for model SYSTEM.
      
      % allow nfEqSystem only with additive and Gaussian noises

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFUDKALMAN';
      p.addRequired('system',@(x) isa(x, 'nefEqSystem'));
      p.parse(system,varargin{:});

      obj@nefEstimator(system,p.Unmatched) % call nefEstimator constructor

      if ~isa(obj.x0,'nefGaussianRV')
        error('NEF:nefKalman:wnotGaussian','Initial condition must be nefGaussianRV object')
      end
      if ~isa(system.w,'nefGaussianRV')
        error('NEF:nefKalman:wnotGaussian','State noise must be nefGaussianRV object')
      end
      if ~isa(system.v,'nefGaussianRV')
        error('NEF:nefKalman:vnotGaussian','Measurement noise must be nefGaussianRV object')
      end

      % change x0 to x0 with UD factored covariance matrix
      obj.x0 = nefGaussianRV(evalMean(obj.x0),nefUDFactorFunction(evalVariance(obj.x0)),'check',0);

      % prepare estQueue
      data.RV = obj.x0;
      obj.estQueue(end+1) = data;

      obj.Ur = [];
      obj.Dr = [];
      obj.Rtilde = [];

      % if diff1Noise is Constant and R is constant then Rtilde is also constant
      obj.isRtildeConst =  obj.sys.h.isDiff1NoiseConst & obj.sys.v.VarNumeric;
    end % UDkalman constructor
    
    function disp(obj)
      % DISP Displays nefUDKalman object parameters and options.
      
      fprintf('A nefUDKalman object with parameters\n')
      showOptions(obj)
    end % function disp
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMEUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predEstimate] = timeUpdate(obj,filtEst,Input,Time)
    % TIMEUPDATE Performs time update step with respect to system
    % description.
      
      % parameters of filtering or previous prediction estimate
      filtMean = evalMean(filtEst.RV,[],[],[],[]);
      filtVar = filtEst.RV.Var;
      %filtVar = evalVariance(filtEst.RV,[],[],[],[]);
      % mean and variance of state noise
      wMean = evalMean(obj.sys.w,[],Input,[],Time);
      wVar = evalVariance(obj.sys.w,[],Input,[],Time);
      % matrices F for state and Gamma for noise in state equation
      F = evalDiff1State(obj.sys.f,filtMean,Input,wMean,Time);
      Gamma = evalDiff1Noise(obj.sys.f,filtMean,Input,wMean,Time);

      % evaluating prediction mean and variance
      [predMean,predVar] = nefUDKalman.UDKalmanTimeUpdate(filtMean,filtVar.U,filtVar.D,obj.sys.f,Input,wMean,Gamma*wVar*Gamma',Time,F);

      % returning nefGaussianRV
      predEstimate.RV = nefGaussianRV(predMean,predVar,'check',0);
    end % function timeUpdate


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASUREMENTUPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtEstimate] = measurementUpdate(obj,predEst,Input,Measurement,Time)
    % MEASUREMENTUPDATE Performs measurement update step with respect to
    % system description.
      
      predMean = evalMean(predEst.RV,[],[],[],[]);
      predVar = predEst.RV.Var;
      %predVar = evalVariance(predEst.RV,[],[],[],[]);
      % mean and variance of measurement noise
      vMean = evalMean(obj.sys.v,[],Input,[],Time);
      vVar = evalVariance(obj.sys.v,[],Input,[],Time);                   % R
      % matrices H for state and Delta for noise in state equation
      H = evalDiff1State(obj.sys.h,predMean,Input,vMean,Time);
      Delta = evalDiff1Noise(obj.sys.h,predMean,Input,vMean,Time);

      % obtain UD factors of Rtilde = Ur*Dr*Ur'
      if obj.isRtildeConst == 1  % 1 = true, 0 = false
        % for constant Rtilde compute its factors at the first step
        if isempty(obj.Rtilde) % first time evaluation; probably Time==0
          obj.Rtilde = Delta*vVar*Delta';
          factored = nefUDFactorFunction(obj.Rtilde);
          obj.Ur = factored.U;
          obj.Dr = factored.D;
          obj.invUr = inv(factored.U);
        end
        % for non-constant Rtilde compute its factors at each step
      else
        obj.Rtilde = Delta*vVar*Delta';
        factored = nefUDFactorFunction(obj.Rtilde);
        obj.Ur = factored.U;
        obj.Dr = factored.D;
        obj.invUr = inv(factored.U);
      end

      [filtMean,filtVar] = nefUDKalman.UDKalmanMeasurementUpdate(predMean,predVar.U,predVar.D,obj.sys.h,Input,vMean,obj.Dr,obj.invUr,Time,Measurement,H);

      % returning nefGaussianRV
      filtEstimate.RV = nefGaussianRV(filtMean,filtVar,'check',0);
    end % function measurementUpdate

  end %methods
  methods(Static)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bierman UD Measurement Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [varargout] = UDKalmanMeasurementUpdate(predMean,predUVar,predDVar,h,Input,vMean,Dr,invUr,Time,Measurement,H)
      % UDKALMANMEASUREMENTUPDATE Implements measurement update of the
      % UD factorized (E)KF.
      
      % transormation of measurement and of matrix H
      Measurement = invUr*Measurement;
      H = invUr*H;

      % filtering mean can be computed at once
      % to omit inversion recursive (sequential) computation is preferred
      %zPredMean =  invUr*h(predMean,Input,vMean,Time);
      %dz = Measurement-zPredMean;
      %zPredVar = Dr + H*predVar*H';
      %invzPredVar = inv(zPredVar);
      %K = predVar*H'*invzPredVar;
      %filtMean = predMean + K*dz;

      % sequential computation of filtering mean and covariance matrix
      % - initialization
      xm = predMean;
      nz = size(Measurement,1);
      %Px = predVar;
      %factor = nefUDFactorFunction(Px);
      %factor = nefUDFactorFunction(predVar);
      %U = factor.U;
      %D = diag(factor.D);
      U = predUVar;
      D = diag(predDVar);
      % - sequential utilisition of measurement
      for i=1:1:nz
        % -- variance of particular measurement prediction
        % -- UDU' = Px
        alpha = H(i,:)*(U*D*U')*H(i,:)' + Dr(i);
        % % -- partial filter gain
        % K = U*(D*U'*H(i,:)')/alpha;
        % -- part of partial filter gain
        Kp = (D*U'*H(i,:)')/alpha;

        % -- partial measurement predicition and innovation
        zPredMean =  invUr*evaluate(h,xm,Input,vMean,Time);
        dz = Measurement-zPredMean;
        % -- partial mean
        % xm = xm + K*dz(i);
        xm = xm + U*Kp*dz(i);
        % -- partial covariance matrix
        % RHM = D-(alpha)^(-1)*(D*U'*H(i,:)')*(D*U'*H(i,:)')';
        RHM = D-Kp*(D*U'*H(i,:)')';
        factor = nefUDFactorFunction((RHM+RHM')/2);
        Ubar = factor.U;
        Dbar = diag(factor.D);
        U = U*Ubar;
        D = Dbar;
        % Px = U*D*U';
      end

      % resulting filtering mean and covariance matrix, storing
      %filtMean = xm;
      %filtVar = nefUDFactorFunction(U,diag(D));
      %filtVar = U*D*U';
      varargout{1} = xm;
      varargout{2} = nefUDFactorFunction(U,diag(D));
      if nargout == 5
        % used for GSM filter, not for a particle filter TODO
        zPredVar = diag(Dr) + H*(U*D*U')*H'; % Pzzp        
        %zPredVar = Dr + H*evaluate(predVar)*H'; % Pzzp    
        varargout{3} = zPredMean;
        varargout{4} = inv(zPredVar);
        varargout{5} = det(zPredVar);
      end;
    end % function UDKalmanMeasurementUpdate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bierman UD Time Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [predMean,predVar] = UDKalmanTimeUpdate(filtMean,filtUVar,filtDVar,f,Input,wMean,Qtilde,Time,F)
    % UDKALMANTIMEUPDATE Implements time update of the UD factorized (E)KF.
      
      % predictive mean computation
      predMean = evaluate(f,filtMean,Input,wMean,Time);

      % predictive covariance matrix computation
      % - initialization
      nx = size(filtMean,1);
      U = filtUVar;
      D = diag(filtDVar);
      W = [F*U, eye(nx)];
      Dhat = [D, zeros(nx); zeros(nx), Qtilde];
      % - computation
      V = zeros(nx,2*nx);
      V(nx,:) = W(nx,:);
      Upred = eye(nx);
      for i=nx-1:-1:1
        su = zeros(1,2*nx);
        for j=i+1:nx
          Upred(i,j) = (W(i,:)*Dhat*V(j,:)')/(V(j,:)*Dhat*V(j,:)');
          su = su + Upred(i,j)*V(j,:);
        end
        V(i,:) = W(i,:) - su;
      end
      Dpred = V*Dhat*V';
      %predVar = Upred*Dpred*Upred';
      predVar = nefUDFactorFunction(Upred,diag(Dpred));
    end; %function UDKalmanTimeUpdate
  end % static methods
end % the class
