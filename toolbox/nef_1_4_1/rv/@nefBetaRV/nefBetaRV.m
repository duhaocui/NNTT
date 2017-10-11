classdef nefBetaRV < nefRV
  %file @nefBetaRV/nefBetaRV.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia


  properties (SetAccess = 'protected') % protected properties 
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVBeta PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%
    Alpha  % mean
    Beta   % covariance matrix
    %%%%%%%%%%%%%%%%%%%%%%%
    % RVGAUSSIAN ATTRIBUTES
    %%%%%%%%%%%%%%%%%%%%%%%
    AlphaNumeric % is numeric
    BetaNumeric %is numeric
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefBetaRV(Alpha,Beta,varargin)
      % NEFBeta Creates nefBetaRV object.
      %
      %   OBJ = NEFBetaRV(Alpha,Beta) creates a NEFBetaRV object OBJ 
      %   representing a Beta random variable with  Alpha and 
      %   Beta parameters
      %
      %   See also NEFRV,NEFGAMMARV, NEFUNIFORMRV,NEFGAUSSIANRV, NEFGAUSSIANSUMRV, NEFPOINTMASSRV.

      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFBETARV';
      p.addRequired('Alpha',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addRequired('Beta',@(x) isnumeric(x) || isa(x,'nefFunction'));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(Alpha,Beta, varargin{:});
      obj.check = p.Results.check;

      if obj.check
        if size(obj.Alpha,1) ~= 1
          error('NEF:nefBetaRV:Alpha','This is scalar version only')
        end
        if size(obj.Beta,1) ~= 1
          error('NEF:nefBetaRV:Beta','This is scalar version only')
        end
        if size(obj.Beta,2) ~= size(obj.Alpha,2)
          error('NEF:nefBetaRV','Size of Alpha nad Beta must be same.')
        end
        try obj = storeSizes(obj,{Alpha,Beta});catch ME1,rethrow(ME1),end
      end

      obj.AlphaNumeric = isnumeric(Alpha);
      obj.BetaNumeric = isnumeric(Beta);

      obj.Alpha = Alpha;
      obj.Beta = Beta;
      obj.dimRV = size(obj.Alpha,1);
      obj.isContinuous = 1;
    end % nefBetaRV constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ***EVALUATE PDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalPDF(obj,point,varargin)
      % EVALPDF Evaluates probability density function at a point.
      %
      %   VAL = EVALPDF(OBJ,POINT,VARAGIN) returns value of the Beta pdf of the random
      %   variable OBJ at a point given by POINT. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(1,[],[]);
      %      OBJ = NEFBetaRV(M,1);
      %      POINT = 1;
      %      VAL = EVALPDF(OBJ,POINT,10,[],[],[]);
      %
      %   See also DRAWSAMPLE.


      [Alpha,Beta] = evaParameters(obj,varargin{:});
      val = BETAPDF(point,Alpha,Beta);
    end % funtion evalPDF

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %*** DRAW SAMPLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = drawSample(obj,count,varargin)
      % ***DRAWSAMPLE Draws samples of the random variable.
      %
      %   VAL = DRAWSAMPLE(OBJ,COUNT,VARAGIN) returns COUNT samples 
      %   of the Beta random variable OBJ. VARARGIN contains values of state, input,
      %   and time that may be necessary for evaluation of the random variable parameters
      %   Note that the parameters must not depend on a noise and thus the noise must be 
      %   specified as empty.
      %
      %   Example
      %      M = NEFLINEARFUNCTION(0.9,1,[]);
      %      OBJ = NEFBetaRV(M,1);
      %      COUNT = 5;
      %      VAL = DRAWSAMPLE(OBJ,COUNT,10,1,[],[]);
      %
      %   See also EVALPDF.

      [Alpha,Beta] = evaParameters(obj,varargin{:});
      val = betarnd(Alpha,Beta,1,count);
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
      [Alpha,Beta] = evaParameters(obj,varargin{:});
      val = Alpha./(Alpha+Beta);
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
      [Alpha,Beta] = evaParameters(obj,varargin{:});
      val = Alpha.*Beta./((Alpha + Beta).^2.*(Alpha+Beta+1));
    end % funtion evalVariance

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefBetaRV object B{Alpha,Beta}\n');
      Alpha = obj.Alpha;
      Beta = obj.Beta;
    end % function disp
  end %methods
  methods (Access = 'private')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALPARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Alpha,Beta] = evalParameters(obj,varargin)
      % EVALPARAMETERS Evaluates parameters (including necessary checks)
      if obj.AlphaNumeric
        Alpha = obj.Alpha;
      else
        Alpha = evaluate(obj.Alpha,varargin{:});
      end
      if obj.BetaNumeric
        Beta = obj.Beta;
      else
        Beta = evaluate(obj.Beta,varargin{:});
      end
    end % function evalParameters
  end %methods private
end % classdef
