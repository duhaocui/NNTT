classdef nefGaussianSumRV < nefRV
  %file @nefGaussianSumRV/nefGaussianSumRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIANSUM PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    N@double = 0; % number of terms
    Weights@cell = {}; % weights
    Means@cell = {}; % means
    Vars@cell = {}; % variances
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ADDITIONAL PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    precVarEigVec@cell = {} % eigenvectors of Variances
    precVarEigVal@cell = {} % eigenvalues of Variances
    precVarInv@cell = {} % inversion of Variances
    precVarLowTriag@cell = {} % Lower Triangular Decomposition of Variances
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%
    WNumeric % are weights numerical
    MNumeric % are means numerical
    VNumeric % are variances numerical
    %%%%%%%%%%%%%%%%%%%%%%%
    % type of parameters specification
    %%%%%%%%%%%%%%%%%%%%%%%
    parameters@char
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % nefGaussianSumRV CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefGaussianSumRV(varargin)
      % NEFGAUSSIANSUMRV Creates a nefGaussianSumRV object.
      %
      %   OBJ = NEFGAUSSIANSUMRV(VARARGIN) creates
      %   a NEFGAUSSIANSUMRV object OBJ representing random variable
      %   with Gaussian sum probability density function.
      %
      %   The parameters may be specified as
      %
      %   WEIGHT-MEAN-COVARIANCE MATRIX TRIPLETS
      %   OBJ = NEFGAUSSIANSUMRV(W1,M1,V1,W2,M2,V2,...) creates a NEFGAUSSIANSUMRV
      %   object OBJ representing random variable with Gaussian sum probability
      %   density function. The parameters of the Gaussian sum are given by tripplets:
      %   weight Wi, mean Mi, covariance matrix Vi. Thus the probability density function has
      %   the following form \sum_{i=1}^N Wi * N(Mi,Vi).
      %
      %   THREE CELL ARRAYS OF WEIGHTS, MEANS AND COVARIANCE MATRICES
      %   OBJ = NEFGAUSSIANSUMRV(W,M,V,...) creates a NEFGAUSSIANSUMRV
      %   object OBJ representing random variable with Gaussian sum probability
      %   density function. The parameters of the Gaussian sum are given by weights W,
      %   means M and covariance matrices V. The parameters W,M,V are cell arrays of proper
      %   dimmensions.
      %
      %   WEIGHT-NEFGAUSSIANRV PAIRS
      %   OBJ = NEFGAUSSIANSUMRV(W1,GRV1,W2,GRV2,...) creates a NEFGAUSSIANSUMRV
      %   object OBJ representing random variable with Gaussian sum probability
      %   density function. The parameters of the Gaussian sum are given by pairs
      %   weight Wi, nefGaussianRV object GRVi. Thus the probability density function has
      %   the following form \sum_{i=1}^N Wi * GRVi.
      %
      %   The way of specification is guessed by constructor or may be explicitly specified
      %   as a param-value pair at the end of arguments
      %
      %   See also NEFRV, NEFUNIFORMRV, NEFGAUSSIANRV, NEFPOINTMASSRV.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % SPLIT VARARGIN INTO PARAMETER PART
      % AND PARAM-VALUE PART
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % set pointer to the first string in param-value pairs
      pointer = nargin+1;
      continueSearching = 1;
      while continueSearching && (pointer-2) > 0
        continueSearching = 0;
        if ischar(varargin{pointer-2})
          pointer = pointer-2;
          continueSearching = 1;
        end
      end


      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFGAUSSIANSUMRV';
      p.addParamValue('parameters','undefined',@(x)any(strcmpi(x,{'triplets','wmv','wnefgaussianrv'})));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(varargin{pointer:end});
      obj.check = p.Results.check;
      obj.parameters = p.Results.parameters;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % SELECT PARAMETERS STRUCTURE AND
      % CONVERT THEM TO THE Weights,Means,Vars TRIPLET
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      switch lower(obj.parameters)
        case 'triplets'
          [Weights,Means,Vars] = tripplets2wmv(varargin{1:pointer-1});
        case 'wmv'
          [Weights,Means,Vars] = wmv2wmv(varargin{1:pointer-1});
        case 'wnefgaussianrv'
          [Weights,Means,Vars] = wnefGaussianRV2wmv(varargin{1:pointer-1});
        case 'undefined'
          warning('NEF:GaussianSumRV:undefinedSpecification','nefGaussianSumRV: Undefined specification of parameters - trying to guess')
          % OK, now the hard part, guessing
          nPar = pointer-1;
          % if nPar is 3 and first par is a vector -> try wmv structure
          if (nPar == 3) && isvector(varargin{1}) && iscell(varargin{1})
            [Weights,Means,Vars] = wmv2wmv(varargin{1:nPar});
            % if mod nPar,2 is 0 and second par is nefGaussianRV
          elseif (mod(nPar,2) == 0) && isa(varargin{2},'nefGaussianRV')
            [Weights,Means,Vars] = wnefGaussianRV2wmv(varargin{1:nPar});
          elseif mod(nPar,3)== 0
            [Weights,Means,Vars] = tripplets2wmv(varargin{1:nPar});
          else
            error('NEF:nefGaussianSumRV:GivingUpGuessing','Giving up guessing parameters specification, please use parameter "parameters"')
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%
      % STORE PARAMETERS
      %%%%%%%%%%%%%%%%%%%%%%%
      obj.Weights = Weights;
      obj.Means = Means;
      obj.Vars = Vars;
      obj.N = length(Weights);
      if obj.check
        %length(Weights)
        checkWMV(obj);
      end
      % set attributes
      % after checkWMV all Weights are either numeric of nefFunction
      % the isequal applies to Means and Vars
      obj.WNumeric = isnumeric(obj.Weights{1});
      obj.MNumeric = isnumeric(obj.Means{1});
      obj.VNumeric = isnumeric(obj.Vars{1});
      try obj = storeSizes(obj,{obj.Weights{:} obj.Means{:} obj.Vars{:}});catch ME1,rethrow(ME1),end

      % if variances are numeric, check them here
      if obj.check 
        if obj.VNumeric
          for i = 1:obj.N
            try nefRV.checkCovarianceMatrix(obj.Vars{i});catch ME1,rethrow(ME1),end
          end
        end
        % if weights are numeric, check them here
        if obj.WNumeric
          try nefRV.checkWeights(Weights);catch ME1,rethrow(ME1),end
        end
      end

      obj.dimRV = size(obj.Means{1},1);

      % initialize cell arrays for precomputation
      obj.precVarEigVec = cell(obj.N);
      obj.precVarEigVal = cell(obj.N);
      obj.precVarInv = cell(obj.N);
      obj.precVarLowTriag= cell(obj.N);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % CONVERSION FUNCTIONS DEFINITION
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [Weights,Means,Vars] = tripplets2wmv(varargin)
        % TRIPPLETS2WMV
        %
        %   Weights,Means,Vars = TRIPPLETS2WMV(VARARGIN)
        %
        %   See also.

        if obj.check && (mod(nargin,3) > 0)
          error('NEF:nefGaussianSumRV:TrippletsParNum','number of random variable parameters must be a multiple of three')
        end
        div3 = nargin/3;
        Weights = cell(1,div3);
        Means = cell(1,div3);
        Vars = cell(1,div3);
        for i = 1:div3
          Weights{i} = varargin{3*i-2};
          Means{i} = varargin{3*i-1};
          Vars{i} = varargin{3*i};
        end
      end % tripples2wmv constructor


      function [Weights,Means,Vars] = wmv2wmv(varargin)
        % WMV2WMV
        %
        %   Weights,Means,Vars = WMV2WMV(VARARGIN)
        %
        %   See also.

        if obj.check && (nargin ~=3)
          error('NEF:nefGaussianSumRV:WMVParNum','Three random variable parameters required')
        end
        tmpW = varargin{1};
        tmpM = varargin{2};
        tmpV = varargin{3};
        lengths(1)=length(tmpW);
        lengths(2)=length(tmpM);
        lengths(3)=length(tmpV);

        if obj.check && (sum(mean(lengths) ~= lengths))
          error('NEF:nefGaussianSumRV:WMVLengths','incompatible lengths of random variable parameters')
        end

        Weights = cell(1,lengths(1));
        Means = cell(1,lengths(2));
        Vars = cell(1,lengths(3));
        for i = 1:length(tmpW)
          Weights{i} = tmpW{i};
          Means{i} = tmpM{i};
          Vars{i} = tmpV{i};
        end
      end % wmv2wmv constructor

      function [Weights,Means,Vars] = wnefGaussianRV2wmv(varargin)
        % WNEFGAUSSIANRV2WMV
        %
        %   Weights,Means,Vars = WNEFGAUSSIANRV2WMV(VARARGIN)
        %
        %   See also.
        if obj.check && (mod(nargin,2) > 0)
          error('NEF:nefGaussianSumRV:WnefGaussParNum','Number of random variable parameters must be a multiple of two')
        end
        N = nargin/2;
        Weights = cell(1,N);
        Means = cell(1,N);
        Vars = cell(1,N);
        for i = 1:N
          Weights{i} = varargin{2*i-1};
          tmpobj = varargin{2*i};
          if obj.check && ~isa(tmpobj,'nefGaussianRV')
            error('NEF:nefGaussianSumRV:NotnefGauss','not a nefGaussianRV object')
          end
          Means{i} = tmpobj.Mean;
          Vars{i} = tmpobj.Var;
        end
      end % wnefGaussianRV2wmv constructor
    end % nefGaussianSumRV constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE PDF AT THE SPECIFIED POINT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varargin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,VARAGIN) returns value of the Gaussian sum pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      W1 = 0.5;
      %      M1 = NEFLINEARFUNCTION(1,[],[]);
      %      V1 = 1;
      %      W2 = 0.5;
      %      M2 = NEFLINEARFUNCTION(2,[],[]);
      %      V2 = 1;
      %      OBJ = NEFGAUSSIANSUMRV('tripplets',W1,M1,V1,W2,M2,V2);
      %      POINT = 0;
      %      VAL = EVALPDF(OBJ,POINT,10,[],[],[]);
      %
      %   See also DRAWSAMPLE.

      [dim,N] = size(point);
      if obj.check 
        if (dim ~= obj.dimRV)
          error('NEF:nefGaussianSumRV:PointDim','Wrong dimension of point')
        end
        checkVarargin(obj,varargin{:})
      end

      val = zeros(1,N);
      [Weights,Means,Vars] = evalParameters(obj,varargin{:});

      for i = 1:obj.N
        % if Vars is positive definite
        detV = det(Vars(:,:,i));
        if detV > 0
          if obj.VNumeric % in this case the inverse mey already be calculated
            if isempty(obj.precVarInv{i})
              % compute inverse matrix and store it for later use
              obj.precVarInv{i} = inv(Vars(:,:,i));
            end
            invV = obj.precVarInv{i}; % invV is precomputed in this case
          else % its no use to store it, as it is a function
            invV = inv(Vars(:,:,i));
          end
          [val] = val + Weights(i)*nefGaussianRV.evalGaussianPDF(point,Means(:,i),Vars(:,:,i),0,detV,invV,[]);
        else % if Vars is singular
          if obj.VarNumeric % v and d are precomputed in this case
            if isempty(obj.precVarEigVec{i})
              % compute positive eigenvectors and eigenvalues and store for later use
              [obj.precVarEigVec{i},obj.precVarEigVal{i}] = eig(Vars(:,:,i));
            end
            v = obj.precVarEigVec{i};
            d = obj.precVarEigVal{i};
          else
            [v,d] = eig(Vars(:,:,i));
          end
          [val] = val + Weights(i)*nefGaussianRV.evalGaussianPDF(point,Means(:,i),Vars(:,:,i),0,detV,[],[v d]);
        end % if
      end % for
    end % funtion evalPDF
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
          error('NEF:nefGaussianSumRV:PointDim','Wrong dimension of point')
        end
        if obj.check && (any(idx>obj.dimRV) | any(idx<=0))
          error('NEF:nefGaussianSumRV:WrongIdx','Wrong indices of the component')
        end
        checkVarargin(obj,varargin{:})
      end
      val = zeros(1,N);
      [Weights,Means,Vars] = evalParameters(obj,varargin{:});
      for i = 1:obj.N
        % get components belonging to the marginal
        [val] = val + Weights(i)*nefGaussianRV.evalGaussianPDF(point,Means(idx,i),Vars(idx,idx,i),0,[],[],[]);
      end % for
    end % funtion evalMarginalPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW SAMPLE(S) OF THE RV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      %      W1 = 0.5;
      %      M1 = NEFLINEARFUNCTION(0.9,1,[]);
      %      V1 = 1;
      %      W2 = 0.5;
      %      M2 = NEFLINEARFUNCTION(1.9,1,[]);
      %      V2 = 1;
      %      OBJ = NEFGAUSSIANSUMRV('tripplets',W1,M1,V1,W2,M2,V2);
      %      COUNT = 5;
      %      VAL = DRAWSAMPLE(OBJ,COUNT,10,1,[],[]);
      %
      %   See also EVALPDF.

      if obj.check
        checkVarargin(obj,varargin{:})
      end
      val = zeros(obj.dimRV,count);
      [Weights,Means,Vars] = evalParameters(obj,varargin{:});

      idx = nefRV.rndMultinomial(Weights,count);
      for i = 1:count
        if obj.VNumeric
          if isempty(obj.precVarLowTriag{idx(i)})
            % compute Lower Triangular Decomposition of Variance and store for later use
            obj.precVarLowTriag{idx(i)} = nefCholFactorFunction.CholFactor(Vars(:,:,idx(i)));
          end
          L = obj.precVarLowTriag{idx(i)};
        else
          L = nefCholFactorFunction.CholFactor(Vars(:,:,idx(i)));
        end
        val(:,i) = Means(:,idx(i)) + L*randn(obj.dimRV,1);
      end
    end % funtion drawSample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MEAN OF THE RV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMean(obj,varargin)
      % EVALMEAN Computes mean of the Gaussian sum RV.
      %
      %   VAL = EVALMEAN(OBJ,VARARGIN) returns numerical value of the
      %   mean given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALVARIANCE.

      [Weights,Means,Vars] = evalParameters(obj,varargin{:});
      val = Means * Weights';
    end % funtion evalMean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE COVARIANCE MATRIX OF THE RV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalVariance(obj,varargin)
      % EVALVARIANCE Evaluates covariance matrix of the random variable.
      %
      %   VAL = EVALVARIANCE(OBJ,VARARGIN) returns numerical value of the
      %   covariance matrix given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   See also EVALMEAN.

      [Weights,Means,Vars] = evalParameters(obj,varargin{:});
      % TOTAL MEAN
      TM = Means * Weights';
      % TOTAL COVARIANCE MATRIX
      val = -TM*TM';
      for i = 1:obj.N
        val = val + Weights(i)*(Vars(:,:,i)+Means(:,i)*Means(:,i)');
      end
    end % funtion variance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRINT OBJECT INFO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefGaussianSumRV object\n');
    end % function disp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISABLECHECKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disableChecks(obj)
      obj.check = 0;

      if ~obj.WNumeric
        for i = 1:obj.N
          disableChecks(obj.Weights{i});
        end
      end
      if ~obj.MNumeric
        for i = 1:obj.N
          disableChecks(obj.Means{i});
        end
      end
      if ~obj.VNumeric
        for i = 1:obj.N
          disableChecks(obj.Vars{i});
        end
      end
    end % function disableChecks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALWEIGHTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Weights] = evalWeights(obj,idx,varargin)
      % EVALWEIGHTS Evaluates weights
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      Weights = zeros(1,length(idx));
      if obj.WNumeric
        Weights = cell2mat({obj.Weights{idx}});
      else
        for i = 1:length(idx)
          Weights(i) = evaluate(obj.Weights{idx(i)},varargin{:});
        end
      end
      if obj.check && ~obj.WNumeric
        try nefRV.checkWeights(Weights);catch ME1,rethrow(ME1),end
      end
    end % function evalWeights
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALMEANS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Means] = evalMeans(obj,idx,varargin)
      % EVALMEANS Evaluates means
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      Means = zeros(obj.dimRV,length(idx));
      if obj.MNumeric
        Means = cell2mat({obj.Means{idx}});
      else
        for i = 1:length(idx)
          Means(:,i) = evaluate(obj.Means{idx(i)},varargin{:});
        end
      end
    end % function evalMeans
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALVARS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Vars] = evalVars(obj,idx,varargin)
      % EVALVARS Evaluates covariance matrices
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      Vars = zeros(obj.dimRV,obj.dimRV,length(idx));
      if obj.VNumeric
        for i = 1:length(idx)
          Vars(:,:,i) = cell2mat({obj.Vars{idx(i)}});
        end
      else
        for i = 1:length(idx)
          Vars(:,:,i) = evaluate(obj.Vars{idx(i)},varargin{:});
          if obj.check
            try nefRV.checkCovarianceMatrix(Vars(:,:,i));catch ME1,rethrow(ME1),end
          end
        end
      end
    end % function evalVars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evalsvars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [SVars] = evalSVars(obj,idx,varargin)
      % evalvars evaluates square root of covariance matrices
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      SVars = zeros(obj.dimRV,obj.dimRV,length(idx));
      if obj.VNumeric
        for i = 1:length(idx)
          SVars(:,:,i) = nefCholFactorFunction.CholFactor(cell2mat({obj.Vars{idx(i)}}));
        end
      else
        for i = 1:length(idx)
          if isa(obj.Vars{idx(i)},'nefCholFactorFunction') % if is is already factored
            SVars(:,:,i) = obj.Vars{idx(i)}.S;
          else % otherwise compute it
            SVars(:,:,i) = nefCholFactorFunction.CholFactor(evaluate(obj.Vars{idx(i)},varargin{:}));
          end
        end
      end
    end % function evalSVars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evalsvars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Us,Ds] = evalUDVars(obj,idx,varargin)
      % evalvars evaluates UD factor of covariance matrices
      if obj.check
        checkVarargin(obj,varargin{:})
      end
      Us = zeros(obj.dimRV,obj.dimRV,length(idx));
      Ds = zeros(obj.dimRV,1,length(idx));
      if obj.VNumeric
        for i = 1:length(idx)
          [Us(:,:,i),Ds(:,:,i)] = nefUDFactorFunction.UDFactor(cell2mat({obj.Vars{idx(i)}}));
        end
      else
        for i = 1:length(idx)
          if isa(obj.Vars{idx(i)},'nefUDFactorFunction') % if is is already factored
            Us(:,:,i) = obj.Vars{idx(i)}.U;
            Ds(:,:,i) = obj.Vars{idx(i)}.D;
          else % otherwise compute it
            [Us(:,:,i),Ds(:,:,i)] = nefUDFactorFunction.UDFactor(evaluate(obj.Vars{idx(i)},varargin{:}));
          end
        end
      end
    end % function evalUDVars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALPARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Weights,Means,Vars] = evalParameters(obj,varargin)
      % EVALPARAMETERS Evaluates parameters (including necessary checks)
      idx = 1:obj.N;
      [Weights] = evalWeights(obj,idx,varargin{:});
      [Means] = evalMeans(obj,idx,varargin{:});
      [Vars] = evalVars(obj,idx,varargin{:});
    end % function evalParameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check wheter the Covariance
    % is not given as function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function isNumeric = isVarNumeric(obj)
      % isVarNumeric Checks wheter the Covariance is not given as function
      isNumeric = (obj.WNumeric && obj.MNumeric && obj.VNumeric)
    end % function isVarNumeric
  end %methods

  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK CONSISTENCY OF THE TRIPLET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkWMV(obj)
      % CHECKWMV Check consistency of WMV parameters.
      %
      %   CHECKWM(OBJ) checks consistency and types of weights, means and
      %   covariance matrices.
      %
      %   See also.

      Weights = obj.Weights;
      Means = obj.Means;
      Vars = obj.Vars;
      Wdim = zeros(obj.N,2);
      Mdim = zeros(obj.N,2);
      Vdim = zeros(obj.N,2);
      WFun = zeros(obj.N,1);
      MFun = zeros(obj.N,1);
      VFun = zeros(obj.N,1);
      WNumeric = zeros(obj.N,1);
      MNumeric = zeros(obj.N,1);
      VNumeric = zeros(obj.N,1);
      % object type check
      for i = 1:obj.N
        Wdim(i,:) = size(Weights{i});
        Mdim(i,:) = size(Means{i});
        Vdim(i,:) = size(Vars{i});
        WNumeric(i) = isnumeric(Weights{i});
        WFun(i) = isa(Weights{i},'nefFunction');
        MNumeric(i) = isnumeric(Means{i});
        MFun(i) = isa(Means{i},'nefFunction');
        VNumeric(i) = isnumeric(Vars{i});
        VFun(i) = isa(Vars{i},'nefFunction');
      end

      % checking types
      if (sum(WNumeric) ~= obj.N) && (sum(WFun) ~= obj.N)
        error('NEF:nefGaussianSumRV:WeightsType','Weights must be either all numeric or all nefFunction')
      end
      if (sum(MNumeric) ~= obj.N) && (sum(MFun) ~= obj.N)
        error('NEF:nefGaussianSumRV:MeansType','Means must be either all numeric or all nefFunction')
      end
      if (sum(VNumeric) ~= obj.N) && (sum(VFun) ~= obj.N)
        error('NEF:nefGaussianSumRV:VariancesType','Variances must be either all numeric or all nefFunction')
      end

      % checking weights
      if ~isequal(Wdim,ones(obj.N,2)) % both size(Weights,1) and size(Weights,2) must be 1 for each weight
        error('NEF:nefGaussianSumRV:ScalarWeights','Weights must be scalar')
      end
      % checking means
      if ~isequal(Mdim(:,2),ones(obj.N,1))  % size(Means,2) must be 1 for each mean
        error('NEF:nefGaussianSumRV:VectorMean','Mean must be a column vector')
      end
      if any(Mdim(:,1) ~= Mdim(1,1)) % size(Means,1) must be isequal for all means
        error('NEF:nefGaussianSumRV:InconsMean','Inconsistent mean sizes')
      end
      % checking variances
      if ~isequal(Vdim(:,1),Vdim(:,2)) % size(Vars,1) and size(Vars,2) must be equal for each variance matrix
        error('NEF:nefGaussianSumRV:SquareCovar','Each covariance matrix must be square')
      end
      if any(Vdim(:,1) ~= Vdim(1,1)) % size(Vars,1) must be isequal for all variance matrix
        error('NEF:nefGaussianSumRV:InconsCovar','Incompatible covariance matrix sizes')
      end
      % checking means with respect to variances
      if Mdim(1,1) ~= Vdim(1,2) % means and weights must correcpond
        error('NEF:nefGaussianSumRV:InconsCovarMean','Incompatible covariance matrix sizes and mean sizes')
      end

    end % checkWMV function

  end %methods private
  methods(Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDUCTION TECHNIQUES - Threshold Pruning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [redWeights,redMeans,redVars] = pruningT(Weights,Means,Vars,percentOfMaxWeight)
      if (percentOfMaxWeight <= 0 ) || (percentOfMaxWeight > 100)
        error('NEF:nefGaussianSumRV:notvalidValue','The parameter has to be in the interval (0,100> to give percents!');
      end
      cWeights = cell2mat(Weights);
      threshold = percentOfMaxWeight * max(cWeights) / 100;
      idx = find(cWeights>=threshold);
      %                 if isempty(idx) % should not happen
      %                   error('NEF:nefGaussianSumRV:prunedAll','Pruned all terms')
      %                 end
      redWeights = num2cell(nefRV.normalizeWeights(cWeights(idx)));
      redMeans = Means(idx);
      redVars = Vars(idx);
    end % function pruningT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDUCTION TECHNIQUES - Maximum Posterior Probability Pruning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [redWeights,redMeans,redVars] = pruningMPP(Weights,Means,Vars,N)
      if N >= length(Weights)
        %warning('NEF:nefGaussianSumRV:lessWeightThanSpecified','Less or equal terms in mixture than the specified maximal number. No pruning done.')
        %fprintf('Less or equal terms in mixture than the specified maximal number. No pruning done.\n');
        redWeights = Weights;
        redMeans = Means;
        redVars = Vars;
      else
        [sWeights,idx] = sort(cell2mat(Weights),'descend');
        redWeights = num2cell(nefRV.normalizeWeights(sWeights(1:N)));
        redMeans = Means(idx(1:N));
        redVars = Vars(idx(1:N));
      end
    end % function pruningMPP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDUCTION TECHNIQUES - Cumulative Weights Pruning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [redWeights,redMeans,redVars] = pruningCW(Weights,Means,Vars,Cum)
      [sWeights,idx] = sort(cell2mat(Weights),'descend');
      csWeights = cumsum(sWeights);
      nWeights = max(length(find(csWeights<=Cum)),1); % at least one weight must be present
      redWeights = num2cell(nefRV.normalizeWeights(sWeights(1:nWeights)));
      redMeans = Means(idx(1:nWeights));
      redVars = Vars(idx(1:nWeights));
    end % function pruningCW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDUCTION TECHNIQUES - Generalized Pseudo Bayesian 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [redWeights,redMeans,redVars] = pruningGPB1(Weights,Means,Vars)
      % TOTAL MEAN
      TM = zeros(size(Means{1}));
      for i = 1:length(Weights);
        TM = TM + Weights{i}*Means{i};
      end
      % TOTAL COVARIANCE MATRIX
      val = -TM*TM';
      for i = 1:length(Weights)
        val = val + Weights{i}*(Vars{i}+Means{i}*Means{i}');
      end
      redWeights = {1};
      redMeans = {TM};
      redVars = {val};
    end % function reductionGPB1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REDUCTION TECHNIQUES - Distance based pruning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [redWeights,redMeans,redVars] = pruningDBP(Weights,Means,Vars,R)
      if (R < 0 ) || (R > 1)
        error('NEF:nefGaussianSumRV:notvalidValue','The decrease rate of distance based pruning has to be in the interval <0,1>!');
      end
      [sWeights,idx] = sort(cell2mat(Weights),'descend');
      sMeans = Means(idx);
      sVars = Vars(idx);
      nWeights = length(Weights);
      %----------------------------------------------------------
      % compute the distance between the individual mixture terms
      %----------------------------------------------------------
      D = zeros(nWeights,nWeights); % initialize distances storage
      % compute the higher triangular and diagonal part of the distances
      % matrix
      for i = 1:nWeights
        for j = 1:i
          %TODO the following line is strange
          D(i,j) = sWeights(i)*sWeights(j)*nefGaussianRV.evalGaussianPDF(sMeans{i},sMeans{j},sVars{i}+sVars{j},0,[],[],[]); nefGaussianRV.evalGaussianPDF(sMeans{i},sMeans{j},sVars{i}+sVars{j},0,[],[],[]);
          %D(i,j) = sWeights(i)*sWeights(j)/sqrt(2*pi*(sVars{i}+sVars{j}))*exp(-0.5*(sMeans{i}-sMeans{j}).^2/(sVars{i}+sVars{j}));
        end
      end
      % and copy it to the lower trianular part
      for i = 1:nWeights
        for j = i+1:nWeights
          D(i,j) = D(j,i);
        end
      end

      totalSum = sum(sum(D));
      idxAcc = 1;  % initialize list of accepted term with the higest weighted term
      idxRej = []; % initialize list of rejected term (its empty at start of the pruning)

      % get the distance of the term with the highest weight
      dMax = totalSum -2*sum(sum(D(idxAcc,:)))/sum(sWeights(idxAcc))...
        +sum(sum(D(idxAcc,idxAcc)))/sum(sWeights(idxAcc))^2;
      % test the remaining term if they contribute significantly to the pruned mixture            
      for i = 2:nWeights
        tempIdxAcc = [idxAcc i]; % temporalily add the term to the new mixture
        dCurrent = totalSum -2*sum(sum(D(tempIdxAcc,:)))/sum(sWeights(tempIdxAcc))...
          +sum(sum(D(tempIdxAcc,tempIdxAcc)))/sum(sWeights(tempIdxAcc)).^2;
        % test the condition whether the current term has
        % considerable influence on the reduction of the distance
        if dCurrent < dMax*exp(-R*length(tempIdxAcc)) %(1-0.05*length(tempIdxAcc));
          idxAcc = [idxAcc i];
        else
          idxRej = [idxRej i];
        end
      end

      redWeights = num2cell(sWeights(idxAcc)/sum(sWeights(idxAcc)));
      redMeans = sMeans(idxAcc);
      redVars = sVars(idxAcc);
    end % function reductionDBP
  end % methods static

end % classdef
