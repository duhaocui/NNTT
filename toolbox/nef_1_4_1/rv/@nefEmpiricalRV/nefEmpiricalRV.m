classdef nefEmpiricalRV < nefRV
  %file @nefEmpiricalRV/nefEmpiricalRV.m
  
  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia
  
  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%
    % PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    values = [];
    weights = [];
    sampleSize= 0;
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefEmpiricalRV(values,weights,varargin)
      % NEFPOINTMASSRV Creates a NEFPOINTMASS object.
      %
      %   OBJ = NEFPOINTMASSRV(VALUES,WEIGHTS) creates a NEFPOINTMASS object OBJ
      %   representing a point mass random variable given by points VALUES and 
      %   masses WEIGHTS. VALUES must be a [D,N] matrix of N vectors of dimension D.
      %   WEIGHTS must be a row vector of length N.
      %
      %   See also NEFRV, NEFGAUSSIANRV, NEFGAUSSIANSUMRV, NEFUNIFORMRV.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFEMPIRICALRV';
      p.addRequired('values',@(x) isnumeric(x));
      p.addRequired('weights',@(x) isnumeric(x));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(values,weights,varargin{:});
      obj.check = p.Results.check;

      if obj.check
        if (size(values,2) ~= length(weights))
          error('NEF:nefEmpiricalRV:InconValWe','Inconsistent VALUES and WEIGHTS sizes.')
        end
        if (size(weights,1) ~= 1) 
          error('NEF:nefEmpiricalRV:Weights','WEIGHTS must be a row vector.')
        end
        if (weights < 0)
          error('NEF:nefEmpiricalRV:NegWeights','WEIGHTS must be positive.')
        end
      end
      obj.values = values;
      obj.weights = nefRV.normalizeWeights(weights);
      obj.sampleSize = length(weights);
      obj.dimRV = size(values,1);
    end % nefEmpiricalRV constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varargin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,VARAGIN) returns value of the point mass pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   NOTE: as nefEmpiricalRV consist of constant values only, values in VARARGIN will 
      %   not affect the evalPDF value
      %
      %   See also DRAWSAMPLE.
      [dim,N] = size(point);
      if obj.check && (dim ~= obj.dimRV)
        error('NEF:nefGaussianRV:PointDim','Wrong dimension of point')
      end
      val = zeros(1,N);
      for i = 1:N
        if any(all(obj.values == repmat(point(:,i),1,obj.sampleSize),1));
          val(i) = Inf;
        else
          val(i) = 0;
        end
      end
    end % function evalPDF

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
      if obj.check && (dim ~= length(idx))
        error('NEF:nefGaussianRV:PointDim','Wrong dimension of point')
      end
      if obj.check && (any(idx>obj.dimRV) | any(idx<=0))
        error('NEF:nefGaussianRV:WrongIdx','Wrong indices of the component')
      end
      val = zeros(1,N);
      for i = 1:N
        if any(all(obj.values(idx,:) == repmat(point(:,i),1,obj.sampleSize),1));
          val(i) = Inf;
        else
          val(i) = 0;
        end
      end
    end % funtion evalMarginalPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW SAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varargin)
      % DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples 
      %   of the point mass random variable OBJ. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   NOTE: as nefEmpiricalRV consist of constant values only, values in VARARGIN will 
      %   not affect the evalPDF value
      %
      %   See also EVALPDF.

      idx = nefRV.rndMultinomial(obj.weights,count);
      val = obj.values(:,idx);obj.dimRV,ones(1,count);
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
      %   NOTE: as nefEmpiricalRV consist of constant values only, values in VARARGIN will 
      %   not affect the evalPDF value
      %
      %   See also EVALVARIANCE.

      val = obj.values*obj.weights';
    end % function evalMean

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
      %   NOTE: as nefEmpiricalRV consist of constant values only, values in VARARGIN will 
      %   not affect the evalPDF value
      %
      %   See also EVALMEAN.

      M = obj.values*obj.weights';
      V = zeros(obj.dimRV);
      for i = 1:obj.sampleSize
        V  = V + (obj.values(:,i)*obj.values(:,i)')*obj.weights(i);
      end
      val = V - M*M';
    end % funtion evalVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EFFICIENT SAMPLE SIZE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ess = efficientSampleSize(obj)
      % EFFICIENTSAMPLESIZE Computes efficient sample size.
      %
      %   ESS = EFFICIENTSAMPLESIZE computes efficient sample size
      %   of the RV. It is given as N/(1+sum(weights^2))
      %
      % See also.

      ess = obj.sampleSize/(1+sum(obj.weights.^2));
    end % function efficientSampleSize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NORMALIZE WEIGHTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function nobj = normalize(obj)
      % NORMALIZE Normalizes weights.
      %
      %   NOBJ = NORMALIZE(OBJ) normalizes weights of the random variable.
      %
      %   See also.

      nobj = obj;
      nobj.weights = nefRV.normalizeWeights(obj.weights);
    end % function normalize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefEmpiricalRV object\n');
    end % function disp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MULTINOMIAL RESAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [resobj] = multinomialResampling(obj)
      values = obj.values;
      idx = nefRV.rndMultinomial(obj.weights,obj.sampleSize);
      resobj = nefEmpiricalRV(values(:,idx),ones(1,obj.sampleSize)/obj.sampleSize);
    end % function multinomialResampling

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESIDUAL RESAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [resobj] = residualResampling(obj)
      idx = [];
      values_count = floor(obj.sampleSize*obj.weights);
      for i = 1:obj.sampleSize
        idx = [idx i*ones(1,values_count(i))];
      end
      if length(idx)<obj.sampleSize
        residual_weights = obj.weights - values_count/obj.sampleSize;
        residual_weights = residual_weights/sum(residual_weights);
        idx = [idx nefRV.rndMultinomial(residual_weights,obj.sampleSize-length(idx))];
      end
      resobj = nefEmpiricalRV(values(:,idx),ones(1,obj.sampleSize)/obj.sampleSize);
    end % function residualResampling

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STRATIFIED RESAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [resobj] = stratifiedResampling(obj)
      values = obj.values;
      u = ((0:(obj.sampleSize-1)) + rand(1,obj.sampleSize))/obj.sampleSize;
      idx = zeros(1,obj.sampleSize);
      cumDist = cumsum(obj.weights);
      j = 1;
      for i=1:obj.sampleSize
        while (u(i)>cumDist(j))
          j = j+1;
        end
        idx(i) = j;
      end;
      resobj = nefEmpiricalRV(values(:,idx),ones(1,obj.sampleSize)/obj.sampleSize);
    end % function stratifiedResampling

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SYSTEMATIC RESAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [resobj] = systematicResampling(obj)
      values = obj.values;
      u = ((0:(obj.sampleSize-1)) + rand*ones(1,obj.sampleSize))/obj.sampleSize;
      idx = zeros(1,obj.sampleSize);
      cumDist = cumsum(obj.weights);
      j = 1;
      for i=1:obj.sampleSize
        while (u(i)>cumDist(j))
          j = j+1;
        end
        idx(i) = j;
      end;
      resobj = nefEmpiricalRV(values(:,idx),ones(1,obj.sampleSize)/obj.sampleSize);
    end % function systematicResampling

  end %methods
end % classdef
