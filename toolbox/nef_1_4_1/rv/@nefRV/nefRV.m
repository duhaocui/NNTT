classdef nefRV < handle
  %file @nefRV/nefRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2010 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % DIMENSION OF THE RV
    %%%%%%%%%%%%%%%%%%%%%
    dimRV@double = 0;
    %%%%%%%%%%%%%%%%%%%%%
    % RV PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%
    isContinuous@double = 0;
    %%%%%%%%%%%%%%%%%%%%
    % SIZES OF STATE, INPUT
    % AND NOISE IN PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%
    stateSizes@double = [];
    inputSizes@double = [];
    noiseSizes@double = [];
    timeSizes@double = [];
    %%%%%%%%%%%%%%%%%%%%%
    % TOTAL VARIABLES DIM
    %%%%%%%%%%%%%%%%%%%%%
    dimState@double = [];
    dimInput@double = [];
    dimTime@double = [];
    %%%%%%%%%%%%%%%%%%%%%
    % PERFORM CHECKS
    %%%%%%%%%%%%%%%%%%%%%
    check@double = 1;
  end % properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefRV()
      % NEFRV Creates nefRV object.
      %
      %   OBJ = NEFRV() creates a NEFRV object OBJ
      %
      %   See also NEFUNIFORMRV, NEFGAUSSIANRV, NEFGAUSSIANSUMRV, NEFPOINTMASSRV.
    end % nefrv constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE NATURAL LOGARITHM OF PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalLogPDF(obj,point,varagin)
      % EVALPDF Evaluates natural logarithm of probability density function at a point.
      %
      %   VAL = EVALLOGPDF(OBJ,VARAGIN) returns logarithm of the pdf of the random
      %   variable OBJ at a point POINTVARARGIN.
      %
      %   See also drawSample.

      error('NEF:nefRV:evalLogPDFUndefined','Method evalLogPDF is not implemented for this random variable.')
    end % function evalLogPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varagin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,VARAGIN) returns value of the pdf of the random
      %   variable OBJ at a point POINTVARARGIN.
      %
      %   See also drawSample.

      error('NEF:nefRV:evalPDFUndefined','Method evalPDF is not implemented for this random variable.')
    end % function evalPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW SAMPLE OF THE RV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varagin)
      % DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples 
      %   of the random variable OBJ
      %
      %   See also evalPDF.

      error('NEF:nefRV:drawSampleUndefined','Method drawSample is not implemented for this random variable.')
    end % function drawSample

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
      %   See also EVALMAP EVALMEDIAN EVALVARIANCE.

      error('NEF:nefRV:evalMeanUndefined','Method evalMean is not implemented for this random variable.')
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

      error('NEF:nefRV:evalMAPUndefined','Method evalMAP is not implemented for this random variable.')
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

      error('NEF:nefRV:evalMedianUndefined','Method evalMedian is not implemented for this random variable.')
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
      %   See also EVALMEAN EVALMAP EVALMEDIAN.

      error('NEF:nefRV:evalVarianceUndefined','Method evalVariance is not implemented for this random variable.')
    end % funtion evalVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = plotPDF(obj,lo,up,varargin)
      % PLOTPDF plots pdf of the uniform RV.
      %
      %   VAL = PLOTPDF(OBJ,LO,UP,VARARGIN) plots probability density
      %   function of the uniform random variable OBJ within the interval
      %   given by LO and UP bounds. The random variable parameters are
      %   evaluated using state, input and time specified in VARARGIN.
      %
      %   See also.

      if obj.check && (size(lo,1) ~= obj.dimRV)
        error('NEF:nefRV:LoDim','Wrong dimension of lo')
      end

      if obj.check && (size(up,1) ~= obj.dimRV)
        error('NEF:nefRV:UpDim','Wrong dimension of up')
      end
      switch obj.dimRV
        case 1
          h = (up - lo)/1000;
          x = lo:h:up;
          y = evalPDF(obj,x,varargin{:});
          plot(x',y')
        case 2
          h1 = (up(1,1) - lo(1,1))/100;
          h2 = (up(2,1) - lo(2,1))/100;
          X = (lo(1,1)):h1:(up(1,1));
          Y = (lo(2,1)):h2:(up(2,1));
          Z = zeros(length(X),length(Y));
          for i = 1:length(X)
            for j = 1:length(Y)
              Z(i,j) = evalPDF(obj,[X(i) Y(j)]',varargin{:});
            end
          end
          mesh(X,Y,Z)
        otherwise 
          error('NEF:nefRV:Plot','Only one or two dimensional RV can be plotted')
      end
    end % function plotPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK AND SET DIMRV PROPERTY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = set.dimRV(obj,dimrv)
      if (floor(dimrv) ~= dimrv) || (dimrv < 0)
        error('NEF:nefRV:Dimension','Dimension must be integer greater than zero')
      else
        obj.dimRV = dimrv;
      end
    end % function set.dimRV

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK AND STORE SIZES OF STATE,
    % INPUT NOISE AND TIME OF EACH 
    % PARAMETER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = storeSizes(obj,parList)
      % STORESIZES Stores arguments sizes of all of the random variable parameters.
      %
      %   OBJ = STORESIZES(OBJ,PARLIST) determines arguments sizes of all parameters 
      %   specified in the PARLIST and stores them in parameters of the random variable OBJ
      %
      %   See also.

      for i = 1:length(parList)
        if isnumeric(parList{i})
          obj.stateSizes(i) = 0;
          obj.inputSizes(i) = 0;
          obj.noiseSizes(i) = 0;
          obj.timeSizes(i) = 0;
        else
          obj.stateSizes(i) = parList{i}.dimState;
          obj.inputSizes(i) = parList{i}.dimInput;
          obj.noiseSizes(i) = parList{i}.dimNoise;
          obj.timeSizes(i) = parList{i}.dimTime;
        end
      end
      if obj.check
        % check nonzero state sizes
        if ~isequal(obj.stateSizes(obj.stateSizes>0),circshift(obj.stateSizes(obj.stateSizes>0),[1 1]))
          error('NEF:nefRV:StateSizes','Inconsistent state sizes in parameters')
        end
        % check nonzero input sizes
        if ~isequal(obj.inputSizes(obj.inputSizes>0),circshift(obj.inputSizes(obj.inputSizes>0),[1 1]))
          error('NEF:nefRV:InputSizes','Inconsistent input sizes in parameters')
        end
        % check noise sizes
        if ~isempty(obj.noiseSizes(obj.noiseSizes>0))
          error('NEF:nefRV:DependOnNoise','Random variable parameters must not depend on another random variable')
        end
      end
      % eval variables dim
      obj.dimState = max(obj.stateSizes);
      obj.dimInput = max(obj.inputSizes);
      obj.dimTime = max(obj.timeSizes);

    end % function storeSizes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECIAL SUBSCRIPTED REFERENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        function val = subsref(obj,s)
    %        % Implement a special subscripted reference
    %        switch s(1).type
    %        case '()'
    %            val = evalPDF(obj,s.subs{:});
    %        case '.' % general reference without knowledge of
    %                 % specified property
    %                 % should be general for all child classes
    %            tmpobj = eval(['obj.' s(1).subs]);
    %            ['obj.' s(1).subs]
    %            val = subsref(tmpobj,s(2:end));
    %       % otherwise
    %       %     error(???Specify value for x as obj(x)???)
    %        end
    %        end % subsref

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE OF THE RV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = size(obj,varargin)
      % SIZE Returns size of the random variable.
      %
      %   VAL = SIZE(OBJ,VARARGIN) Size of random variable is returned. The input and output 
      %   arguments can be used in accord with MATLAB function SIZE for a matrix
      %
      %   See also.
      val = size(zeros(obj.dimRV,1),varargin{:});
    end % function size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefRV object\n');
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check wheter the Covariance 
    % is not given as function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function isNumeric = isVarNumeric(obj)
      % isVarNumeric Checks wheter the Covariance is not given as function
      error('NEF:nefRV:isVarNumericUndefined','Method isVarNumeric is not implemented for this random variable.')
    end % function isVarNumeric

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check Varargins passed to evalXXXX functions
    % either none or four arguments with a single value in each allowed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkVarargin(obj,varargin)
      if (nargin-1 == 0)
      elseif(nargin-1 == 4)
        numval = max([size(varargin{1},2) size(varargin{2},2) size(varargin{3},2) size(varargin{4},2)]);
        if numval > 1
          error('NEF:nefRV:TooManyArguments','Only single value for each argument allowed')
        end
      else
        error('NEF:nefRV:WrongNumberOfVarargins','Either none or four arguments evaluation allowed');
      end
    end % function checkVarargin
  end %methods

  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW AND INDEX FROM MULTINOMIAL 
    % DISTRIBUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [idx] = rndMultinomial(w,N)
      % RNDMULTINOMIAL Draws samples from the multinomial pdf
      %
      %   IDX = RNDMULTINOMIAL(W,N) returns N indices in IDX by sampling from 
      %   multinomial pdf with parameters W.
      %
      %   See also.

      idx = zeros(1,N);
      u = fliplr(cumprod(rand(1,N).^(1./(N:-1:1))));
      cumDist = cumsum(w);
      j = 1;
      for i=1:N
        while (u(i)>cumDist(j))
          j = j+1;
        end
        idx(i) = j;
      end;
    end % function rndMultinomial

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK COVARIANCE MATRIX FOR SYMMETRY
    % AND POSITIVE SEMIDEFINITNESS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkCovarianceMatrix(V)
      % CHECKCOVARIANCEMATRIX Checks covariance matrix for symmetry and semidefinite positivity.
      %
      if any(V ~= V')
        error('NEF:nefRV:CovMatSymmetric','Covariance matrix must be symmetric')
      end
      if any(eig(V) < 0)
        error('NEF:nefRV:CovMatPosSemidef','Covariance matrix must be positive semidefinite')
      end
    end % function CheckCovarianceMatrix

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK WEIGHTS IF THEY SUM TO ONE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkWeights(W)
      % CHECKWEIGHT Checks sum of the weights

      if ~(abs(1- sum(cell2mat(W))) < 3*eps)
        %sum(cell2mat(W)) ~= 1
        error('NEF:nefRV:WeightsSumToOne','Weights must sum to one')
      end
    end % function CheckWeights

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NORMALIZE WEIGHTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function wn = normalizeWeights(w)
      % NORMALIZEWEIGHTS Normalizes weights.
      %
      %   WN = NORMALIZEWEIGHTS(W) normalizes W.
      %
      %   See also.

      % implementation of the Kahan summation algorithm
      compensation = 0;
      sum=0;
      for step = 1:length(w)
        y = w(step) - compensation;
        t = sum + y;
        compensation = ( t - sum ) - y;
        sum = t;
      end

      if sum == 0
        warning('NEF:nefRV:zeroWeights','All weights are zero.')
        wn = ones(size(w))/length(w);
      else
        wn = w/sum;
      end
    end % function normalizeWeights

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NORMALIZE WEIGHTS GIVEN THEIR LOGARITHM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function wn = normalizeLogWeights(logw)
      % NORMALIZELOGWEIGTHS Normalizes weights given their logarithm.
      %
      %   WN = NORMALIZELOGWEIGTHS(LOGW) calculates normalized weights
      %   based on their logarithm.
      %
      %   See also.

      % scale the weights so that the greatest is one
      logw = logw - max(logw);
      % calculate actual weights
      wn = exp(logw);
      % normalize
      wn = wn/sum(wn);
    end % function normalizeLogWeights


  end % static methods
end % classdef
