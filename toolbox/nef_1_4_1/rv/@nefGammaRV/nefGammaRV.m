classdef nefGammaRV < nefRV
  %file @nefGammaRV/nefGammaRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    C  % mean
    Lambda % covariance matrix
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%
    CNumeric % is numeric
    LambdaNumeric % is numeric

  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefGammaRV(C,Lambda,varargin)
      % NEFGamma Creates nefGammaRV object.
      %
      %   OBJ = NEFGammaRV(M,V) creates a NEFGammaRV object OBJ 
      %   representing a Gamma random variable with  C and 
      %   Gamma parameters
      %
      %   See also NEFRV,NEFBETARV, NEFUNIFORMRV,NEFGAUSSIANRV, NEFGAUSSIANSUMRV, NEFPOINTMASSRV.


      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFGAMMARV';
      p.addRequired('C',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addRequired('Lambda',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(C,Lambda, varargin{:});
      obj.check = p.Results.check;

      if obj.check
        if (size(obj.C,1) ~= 1)
          error('NEF:nefGammaRV:C','This is scalar version only')
        end
        if (size(obj.Lambda,1) ~= 1)
          error('NEF:nefGammaRV:Lambda','This is scalar version only')
        end

        if (size(obj.Lambda,2) ~= size(obj.C,2))
          error('NEF:nefGammaRV','Size of C nad Lambda must be same.')
        end
        try obj = storeSizes(obj,{C,Lambda});catch ME1,rethrow(ME1),end
      end

      obj.C = C;
      obj.Lambda = Lambda;

      obj.CNumeric= isnumeric(C);
      obj.LambdaNumeric= isnumeric(Lambda);

      obj.dimRV = size(obj.C,1);
      obj.isContinuous = 1;
    end % nefGammaRV constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ***EVALUATE PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varargin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,VARAGIN) returns value of the Gamma pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(1,[],[]);
      %      OBJ = NEFGammaRV(M,1);
      %      POINT = 1;
      %      VAL = EVALPDF(OBJ,POINT,10,[],[],[]);
      %
      %   See also DRAWSAMPLE.

      [C,Lambda] = evaParameters(obj,varargin{:});
      val = gampdf(point,C,Lambda);
    end % funtion evalPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %*** DRAW SAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varargin)
      % ***DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples 
      %   of the Gamma random variable OBJ. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(0.9,1,[]);
      %      OBJ = NEFGammaRV(M,1);
      %      COUNT = 5;utter
      %      VAL = DRAWSAMPLE(OBJ,COUNT,10,1,[],[]);
      %
      %   See also EVALPDF.

      [C,Lambda] = evalParameters(obj,varargin{:});
      val = gamrnd(C,Lambda,1,count);
    end % funtion drawSample

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE MEAN (1st central moment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalMean(obj,varargin)
      % ***EVALMEAN Evaluates mean of the random variable.
      %
      %   VAL = EVALMEAN(OBJ,VARARGIN) returns numerical value of the
      %   mean given state, input and time specified in VARARGIN.
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   See also EVALVARIANCE.
      [C,Lambda] = evaParameters(obj,varargin{:});
      val = C.*Lambda;
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
      [C,Lambda] = evaParameters(obj,varargin{:});
      val = C.* (Lambda.^2);
    end % funtion evalVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefGammaRV object G{C,Lambda}\n');
      C = obj.C;
      Lambda = obj.Lambda;
    end % function disp
  end %methods
  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALPARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [C,Lambda] = evalParameters(obj,varargin)
      % EVALPARAMETERS Evaluates parameters (including necessary checks)
      if obj.CNumeric
        C = obj.C;
      else
        C = evaluate(obj.C,varargin{:});
      end
      if obj.LambdaNumeric
        Lambda = obj.Lambda;
      else
        Lambda = evaluate(obj.Lambda,varargin{:});
      end
    end % function evalParameters
  end %methods private
end % classdef
