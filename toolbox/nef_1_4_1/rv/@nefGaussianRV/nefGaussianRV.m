classdef nefGaussianRV < nefRV
  %file @nefGaussianRV/nefGaussianRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    Mean % mean
    Var  % covariance matrix
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ADDITIONAL PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    precVarEigVec = [] % eigenvectors of Variance
    precVarEigVal = [] % eigenvalues of Variance
    precVarInv = [] % inversion of Variance
    precVarLowTriag = [] % Lower Triangular Decomposition of Variance
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%
    MeanNumeric % is Mean numeric
    VarNumeric  % is Var numeric
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefGaussianRV(m,v,varargin)
      % NEFGAUSSIANRV Creates nefGAUSSIANRV object.
      %
      %   OBJ = NEFGAUSSIANRV(M,V) creates a NEFGAUSSIANRV object OBJ
      %   representing a Gaussian random variable with mean M and
      %   covariance matrix V
      %
      %   See also NEFRV, NEFUNIFORMRV, NEFGAUSSIANSUMRV, NEFPOINTMASSRV.

      % ugly hack to make things faster (skipping creating parser for constructors called with 'check',0)
      if (nargin-2 == 2) && isstr(varargin{1}) && strcmpi(varargin{1},'check') && isnumeric(varargin{2}) && (varargin{2} == 0)
        % skip creating parser
        obj.check = 0;
      else
        p = inputParser;
        p.KeepUnmatched = true;
        p.FunctionName = 'NEFGAUSSIANRV';
        p.addRequired('m',@(x) isnumeric(x) || isa(x,'nefFunction'));
        p.addRequired('v',@(x) isnumeric(x) || isa(x,'nefFunction'));
        p.addParamValue('check',1,@(x)x==0 || x==1);
        p.parse(m,v, varargin{:});
        obj.check = p.Results.check;
      end

      if obj.check
        if (size(m,2) ~= 1)
          error('NEF:nefGaussianRV:Mean','Mean must be a column vector')
        end
        if (size(v,1) ~= size(v,2))
          error('NEF:nefGaussianRV:Covariance','Covariance matrix must be square')
        end
        if (size(v,1) ~= size(m,1))
          error('NEF:nefGaussianRV:InconsMeanCov','Inconsistent mean and covariance matrix')
        end
        try obj = storeSizes(obj,{m,v});catch ME1,rethrow(ME1),end
        % if Variance is numeric, check it (symmetry, positive semidefinitness)
        if isnumeric(v)
          try nefRV.checkCovarianceMatrix(v);catch ME1,rethrow(ME1),end
        end
      end

      obj.MeanNumeric = isnumeric(m);
      obj.VarNumeric = isnumeric(v);

      obj.Mean = m;
      obj.Var = v;

      obj.dimRV = size(m,1);
      obj.isContinuous = 1;
    end % nefGaussianRV constructor


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varargin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,VARAGIN) returns value of the Gaussian pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(1,[],[]);
      %      OBJ = NEFGAUSSIANRV(M,1);
      %      POINT = 1;
      %      VAL = EVALPDF(OBJ,POINT,10,[],[],[]);
      %
      %   See also DRAWSAMPLE.

      [dim,N] = size(point);
      if obj.check
        if (dim ~= obj.dimRV)
          error('NEF:nefGaussianRV:PointDim','Wrong dimension of point')
        end
        checkVarargin(obj,varargin{:})
      end
      val = zeros(1,N);
      M = evalMean(obj,varargin{:});
      V = evalVariance(obj,varargin{:});
      detV = det(V);
      % if V is positive definite
      if detV > 0
        if obj.VarNumeric
          if isempty(obj.precVarInv)
            % compute inverse matrix and store it for later use
            obj.precVarInv = inv(V);
          end
          invV = obj.precVarInv; % invV is precomputed in this case
        else
          invV = inv(V);
        end
        [val] = nefGaussianRV.evalGaussianPDF(point,M,V,0,detV,invV,[]);
      else % if V is singular
        if obj.VarNumeric % v and d are precomputed in this case
          if isempty(obj.precVarEigVec)
            % compute positive eigenvectors and eigenvalues and store for later use
            [obj.precVarEigVec,obj.precVarEigVal] = eig(V);
          end
          v = obj.precVarEigVec;
          d = obj.precVarEigVal;
        else
          [v,d] = eig(V);
        end
        [val] = nefGaussianRV.evalGaussianPDF(point,M,V,0,detV,[],[v d]);
      end % if
    end % funtion evalPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE NATURAL LOGARITHM OF PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalLogPDF(obj,point,varargin)
      % EVALLOGPDF Evaluates natural logarithm of probability density function at a point.
      %
      %   VAL = EVALLOGPDF(OBJ,POINT,VARAGIN) returns logarithm of the Gaussian pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(1,[],[]);
      %      OBJ = NEFGAUSSIANRV(M,1);
      %      POINT = 1;
      %      VAL = EVALLOGPDF(OBJ,POINT,10,[],[],[]);
      %
      %   See also DRAWSAMPLE.

      [dim,N] = size(point);
      if obj.check 
        if (dim ~= obj.dimRV)
          error('NEF:nefGaussianRV:PointDim','Wrong dimension of point')
        end
        checkVarargin(obj,varargin{:})
      end
      val = zeros(1,N);
      M = evalMean(obj,varargin{:});
      V = evalVariance(obj,varargin{:});
      detV = det(V);
      % if V is positive definite
      if detV > 0
        if obj.VarNumeric
          if isempty(obj.precVarInv)
            % compute inverse matrix and store it for later use
            obj.precVarInv = inv(V);
          end
          invV = obj.precVarInv; % invV is precomputed in this case
        else
          invV = inv(V);
        end
        [val] = nefGaussianRV.evalGaussianPDF(point,M,V,1,detV,invV,[]);
      else % if V is singular
        if obj.VarNumeric % v and d are precomputed in this case
          if isempty(obj.precVarEigVec)
            % compute positive eigenvectors and eigenvalues and store for later use
            [obj.precVarEigVec,obj.precVarEigVal] = eig(V);
          end
          v = obj.precVarEigVec;
          d = obj.precVarEigVal;
        else
          [v,d] = eig(V);
        end
        [val] = nefGaussianRV.evalGaussianPDF(point,M,V,1,detV,[],[v d]);
      end % if
    end % funtion evalLogPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MARGINAL PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMarginalPDF(obj,point,idx,varargin)
      % EVALPDF Evaluates marginal probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,IDX,VARAGIN) returns value of the Gaussian pdf of the random
      %   variable OBJ at a point given by POINT. IDX specifies components of the marginal.
      %   VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION([1;2],[],[]);
      %      OBJ = NEFGAUSSIANRV(M,eye(2));
      %      POINT = 1;
      %      IDX = [2];
      %      VAL = EVALMARGINALPDF(OBJ,POINT,IDX,10,[],[],[]);
      %
      %   See also.

      [dim,N] = size(point);
      if obj.check 
        if (dim ~= length(idx))
          error('NEF:nefGaussianRV:PointDim','Wrong dimension of point')
        end
        if (any(idx>obj.dimRV) | any(idx<=0))
          error('NEF:nefGaussianRV:WrongIdx','Wrong indices of the component')
        end
        checkVarargin(obj,varargin{:})
      end
      val = zeros(1,N);
      M = evalMean(obj,varargin{:});
      V = evalVariance(obj,varargin{:});
      % get components belonging to the marginal
      M = M(idx);
      V = V(idx,idx);
      % OK, now continue as usual, but this time nothing is precomputed
      detV = det(V);
      % if V is positive definite
      if detV > 0
        invV = inv(V);
        val = nefGaussianRV.evalGaussianPDF(point,M,V,0,detV,invV,[]);
      else % if V is singular
        [v,d] = eig(V);
        val = nefGaussianRV.evalGaussianPDF(point,M,V,0,detV,[],[v d]);
      end % if
    end % funtion evalMarginalPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW SAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varargin)
      % DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples
      %   of the Gaussian random variable OBJ. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(0.9,1,[]);
      %      OBJ = NEFGAUSSIANRV(M,1);
      %      COUNT = 5;
      %      VAL = DRAWSAMPLE(OBJ,COUNT,10,1,[],[]);
      %
      %   See also EVALPDF.
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      M = evalMean(obj,varargin{:});
      V = evalVariance(obj,varargin{:});
      if obj.VarNumeric
        if isempty(obj.precVarLowTriag)
          % compute Lower Triangular Decomposition of Variance and store for later use
          obj.precVarLowTriag = nefCholFactorFunction.CholFactor(V);
        end
        L = obj.precVarLowTriag;
      else
        L = nefCholFactorFunction.CholFactor(V);
      end
      val = (M*ones(1,count)) + L*randn(obj.dimRV,count);
    end % funtion drawSample

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MEAN (1st central moment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMean(obj,varargin)
      % EVALMEAN Evaluates mean of the random variable.
      %
      %   VAL = EVALMEAN(OBJ,VARARGIN) returns numerical value of the
      %   mean given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALVARIANCE.

      if obj.MeanNumeric
        % replicate obj.Mean by 1 if no varargin is given or by maximum number columns in varargins
        if (nargin-1) == 0 % no varargins
          val = obj.Mean;
        elseif (nargin-1) == 4 % all varargins
          numval = max([size(varargin{1},2) size(varargin{2},2) size(varargin{3},2) size(varargin{4},2)]);
          if numval == 0 % four empty arguments
            val = obj.Mean;
          else
            val = obj.Mean(:,ones(1,numval));
          end
        else
          error('NEF:nefGaussianRV:WrongNumberOfVararginsForMean','Either none or four arguments for mean evaluation allowed');
        end
      else
        val = evaluate(obj.Mean,varargin{:});
      end
    end % funtion evalMean

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MAP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMAP(obj,varargin)
      % EVALMEAN Evaluates maximum aposteriori probability of the random variable.
      %
      %   VAL = EVALMAP(OBJ,VARARGIN) returns numerical value of the
      %   mode given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALMEAN EVALMEDIAN EVALVARIANCE.
      [val] = evalMean(obj,varargin);
    end % funtion evalMAP

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MEDIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMedian(obj,varargin)
      % EVALMEAN Evaluates mean of the random variable.
      %
      %   VAL = EVALMEDIAN(OBJ,VARARGIN) returns numerical value of the
      %   median given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALMEAN EVALMAP EVALVARIANCE.
      [val] = evalMean(obj,varargin);
    end % funtion evalMedian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE VARIANCE (2nd central moment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalVariance(obj,varargin)
      % EVALVARIANCE Evaluates covariance matrix of the random variable.
      %
      %   VAL = EVALVARIANCE(OBJ,VARARGIN) returns numerical value of the
      %   covariance matrix given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALMEAN.

      if obj.VarNumeric
        % replicate obj.Var by 1 if no varargin is given or by maximum number columns in varargins
        if nargin-1 == 0 % no varargins
          val = obj.Var;
        elseif nargin-1 == 4 % all varargins
          numval = max([size(varargin{1},2) size(varargin{2},2) size(varargin{3},2) size(varargin{4},2)]);
          if numval == 0 % four empty arguments
            val = obj.Var;
          else
            val = obj.Var(:,:,ones(1,numval));
          end
        else
          error('NEF:nefGaussianRV:WrongNumberOfVararginsForVar','Either none or four arguments for variance evaluation allowed');
        end
      else
        val = evaluate(obj.Var,varargin{:});
        % check covariance matrix (symmetry, positive semidefinitness)
        if obj.check
          try nefRV.checkCovarianceMatrix(val);catch ME1,rethrow(ME1),end
        end
      end
    end % funtion evalVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISABLECHECKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disableChecks(obj)
      obj.check = 0;

      if ~obj.MeanNumeric
        disableChecks(obj.Mean);
      end
      if ~obj.VarNumeric
        disableChecks(obj.Var);
      end
    end % function disableChecks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefGaussianRV object N{M,V}\n');
      M = obj.Mean
      V = obj.Var
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check wheter the Covariance
    % is not given as function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function isNumeric = isVarNumeric(obj)
      % isVarNumeric Checks wheter the Covariance is numeric
      isNumeric = obj.VarNumeric;
    end % function isVarNumeric

  end %methods
  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALPARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A,B] = evalParameters(obj,varargin)
      % EVALPARAMETERS Evaluates parameters (including necessary checks)
      A = evalMean(obj,varargin{:});
      B = evalVariance(obj,varargin{:});
    end % function evalParameters
  end %methods private

  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evalGaussianPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalGaussianPDF(point,M,V,logPDF,detV,invV,VEigValVec)
      % EVALGAUSSIANPDF Evaluates Gaussian probability density function at a point.
      %
      %   VAL = EVALGAUSSIANPDF(POINT,M,V,LOGPDF,DETV,INVV,VEIGVALVEC))
      %   POINT - points which the pdf should be evaluated at
      %   M - mean
      %   V - covariance matrix
      %   LOGPDF = 1  the function provides the pdf value
      %   LOGPDF = 0  the function provides logarithm of the pdf
      %   DETV Determinant of the covariance matrix if it is known or empty otherwise
      %   INVV Inverse of the covariance matrix if it is known or empty otherwise
      %   VEIGVALVEC Eiggenvalues and vectors of the covariance matrix if they are known or empty otherwise

      [dim,N] = size(point);
      val = zeros(1,N);

      % get det(V)
      if numel(detV) == 0
        detV = det(V);
      else
        detV = detV;
      end

      % if v is positive definite
      if detV > 0
        % get inv(V)
        if numel(invV) == 0
          invV = inv(V);
        else
          invV = invV;
        end
        xm = point-M;
        exponent = -0.5*sum(xm.*(invV*xm),1);
        if logPDF
          val = exponent - (dim/2)*log(2*pi) - 0.5*log(detV);
        else
          val = exp(exponent)/((2*pi)^(dim/2)*sqrt(detV));
        end
      else % if v is singular
        % get eigenvalues and eigenvectors
        if numel(VEigValVec) == 0
          [v,d] = eig(V);
        else
          v = VEigValVec(:,1:dim);
          d = VEigValVec(:,dim+1:2*dim);
        end
        D = diag(d);
        idxpos = D>0;
        idxzer = D==0;
        P = v(:,idxpos); % eigenvectors corresponding to positive eigenvalues
        Q = v(:,idxzer); % eigenvectors corresponding to zero eigenvalues (span null space)
        for i = 1:N
          y = P'*(point(:,i)-M); % center points and transform using eigenvectors corresponding to positive eigenvalues
          z = Q'*(point(:,i)-M); % center points and transform using eigenvectors corresponding to zero eigenvalues
          if z'*z > 0 % if x is not in the hyperplane => pdf=0
            if logPDF
              val(i) = -Inf;
            else
              val(i) = 0;
            end
          else % x is in the hyperplane => compute pdf
            % covariance matrix is given by P'AP, i.e. diagonal matrix of positive eigenvalues
            % if hyperplane is a point, value of PDF is infinite
            if numel(idxpos) == 0
              val(i) = Inf;
            else
              invVp = diag(D(idxpos).^-1);
              detVp = prod(D(idxpos));
              if logPDF
                val(i) = -0.5*(y'*invVp*y) - (length(idxpos)/2)*log(2*pi) - 0.5*log(detVp);
              else
                val(i) = exp(-0.5*y'*invVp*y)/(2*pi)^(length(idxpos)/2)/sqrt(detVp);
              end
            end
          end % if
        end% for
      end % if
    end % funtion evalGaussianPDF
  end % static methods
end % classdef
