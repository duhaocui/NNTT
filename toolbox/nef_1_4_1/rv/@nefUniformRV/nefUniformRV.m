classdef nefUniformRV < nefRV
  %file @nefUniformRV/nefUniformRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVUNIFORM PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    A % lower bound
    B % upper bound
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%
    ANumeric % is A numeric
    BNumeric % is B numeric
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefUniformRV(A,B,varargin)
      % NEFUNIFORMRV Creates a NEFUNIFORMRV object.
      %
      %   OBJ = NEFUNIFORMRV(A,B) creates a NEFUNIFORM object OBJ representing
      %   a uniform random variable with parameters A and B.
      %
      %   See also NEFRV, NEFGAUSSIANRV, NEFGAUSSIANSUMRV, NEFPOINTMASSRV.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFUNIFORMRV';
      p.addRequired('A',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addRequired('B',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(A,B, varargin{:});
      obj.check = p.Results.check;
      % convert possible number to nefconstFunction

      if obj.check
        if (size(A,2) ~= 1)
          error('NEF:nefUniformRV:ColumnPar','A must be column vector')
        end
        if (size(B,2) ~= 1)
          error('NEF:nefUniformRV:ColumnPar','B must be column vector')
        end
        if (size(A,1) ~= size(B,1))
          error('NEF:nefUniformRV:InconsPar','Inconsistent parameters')
        end
        try obj = storeSizes(obj,{A,B});catch ME1 ,rethrow(ME1),end
        % if both parameters are numeric, execute check
        if isnumeric(A) && isnumeric(B)
          checkPar(obj,A,B);
        end
      end

      obj.ANumeric = isnumeric(A);
      obj.BNumeric = isnumeric(B);

      obj.A = A;
      obj.B = B;

      obj.dimRV = size(A,1);
      obj.isContinuous = 1;
    end % nefUniformRV constructor


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
      %      A = NEFLINEARFUNCTION(0.9,1,0);
      %      B = NEFLINEARFUNCTION(1.9,1,0);
      %      OBJ = NEFUNIFORMRV(A,B);
      %      POINT = 1;
      %      VAL = EVALPDF(OBJ,POINT,0.7,0.1,[],[]);
      %
      %   See also DRAWSAMPLE.

      if obj.check 
        if (size(point,1) ~= obj.dimRV)
          error('NEF:nefUniformRV:PointDim','Wrong dimension of point')
        end
        checkVarargin(obj,varargin{:})
      end
      [A,B] = evalParameters(obj,varargin{:});
      val = zeros(1,size(point,2));
      for i = 1:size(point,2)
        if ((point(:,i)) >= A ) && ((point(:,i)) <= B)
          val(i) = 1/(prod(B-A));
        else
          val(i) = 0;
        end
      end
    end % funtion evalPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRAW SAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varargin)
      % DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples
      %   of the uniform random variable OBJ. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be
      %   specified as empty.
      %
      %   Example
      %      A = NEFLINEARFUNCTION(0.9,1,0);
      %      B = NEFLINEARFUNCTION(1.9,1,0);
      %      OBJ = NEFUNIFORMRV(A,B);
      %      COUNT = 5;
      %      VAL = DRAWSAMPLE(OBJ,COUNT,10,1,[],[]);
      %
      %   See also EVALPDF.

      [A,B] = evalParameters(obj,varargin);
      val = A*ones(1,count) + ((B - A)*ones(1,count)).*rand(obj.dimRV,count);
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

      [A,B] = evalParameters(obj,varargin{:});
      val = [(A + B)/2];
    end % funtion evalMean

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
      [A,B] = evalParameters(obj,varargin{:});
      val = diag(((B - A).^2)/12);
    end % function evalVariance

    function disp(obj)
      fprintf('A nefUniformRV object U{A,B}\n');
    end % function disp
  end %methods

  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALPARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A,B] = evalParameters(obj,varargin)
      % EVALPARAMETERS Evaluates parameters (including necessary checks)
      if obj.ANumeric
        A = obj.A;
      else
        A = evaluate(obj.A,varargin{:});
      end

      if obj.BNumeric
        B = obj.B;
      else
        B = evaluate(obj.B,varargin{:});
      end
      if obj.check && (~obj.ANumeric  || ~obj.BNumeric)
        checkPar(obj,A,B)
      end
    end % function evalParameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [v] = checkPar(obj,A,B)
      % CHECKPAR Checks parameters of the uniform random variable.
      %
      %   CHECKPAR(OBJ,A,B) check parameters A and B of the uniform
      %   random variable OBJ.
      %
      %   See also.

      if any((A >= B))
        error('NEF:nefUniformRV:Pars','A must be lower than B')
      end
    end % function checkPar
  end %methods private
end % classdef
