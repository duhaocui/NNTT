classdef nefConstFunction < nefFunction
  %file @nefConstFunction/nefConstFunction.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % CONSTANT
    %%%%%%%%%%%%%%%%%%%%%
    K@double = [];
  end % propeties

  methods
    function [obj] = nefConstFunction(k,varargin)
      % NEFCONSTFUNCTION Creates nefConstFunction object
      %
      %   OBJ = NEFCONSTFUNCTION(K) creates a NEFCONSTFUNCTION object OBJ
      %   representing a constant function
      %       f(X,U,XI,T) = K
      %   K can be either a matrix or a vector
      %
      %   See also NEFFUNCTION, NEFLINEFUNCTION, NEFHANDLEFUNCTION.

      p = inputParser;
      p.FunctionName = 'NEFCONSTFUNCTION';
      p.addRequired('k',@isnumeric);
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(k,varargin{:});
      obj.check = p.Results.check;

      obj.K = k;
      obj.dimFunction = size(k);
      obj.isLinear = 1;
      obj.isAdditive = 1;
      obj.isContiuous = 1;
      obj.isVectorized = 1;
      obj.isDifferentiable1State = 1;
      obj.isDifferentiable2State = 1;
      obj.isDifferentiable1Noise = 1;
      obj.isDifferentiable2Noise = 1;
      obj.isConst = 1;
    end % nefConstFunction constructor

    function [val] = evaluate(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=K, the function always returns K.
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end
      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      val = repmat(obj.K,[1 1 DataCount]);
    end % function evaluate

    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME.
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefConstFunction:FirstStateDerUndefined','First derivative along State undefined')
    end % function evalDiff1State

    function [val] = evalDiff2State(obj,State, Input,Noise,Time)
      % EVALDIFF2STATE Evaluating second derivative with respect to state
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefConstFunction:SecondStateDerUndefined','Second derivative along State undefined')
    end % function evalDiff2State

    function [val] = evalDiff1Noise(obj,State, Input,Noise,Time)
      % EVALDIFF1NOISE Evaluating first derivative with respect to noise
      %
      %   VAL = EVALDIFF1NOISEE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF2NOISE.
      error('NEF:nefConstFunction:FirstNoiseDerUndefined','First derivative along Noise undefined')
    end % function evalDiff1Noise

    function [val] = evalDiff2Noise(obj,State, Input,Noise,Time)
      % EVALDIFF2NOISE Evaluating second derivative with respect to noise
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE.
      error('NEF:nefConstFunction:SecondNoiseDerUndefined','Second derivative along Noise undefined')
    end % function evalDiff2Noise

    function [val] = size(obj,varargin)
      % SIZE Returns size of the function.
      %
      %   VAL = SIZE(OBJ,VARARGIN) Size of K is returned. The input and output 
      %   arguments can be used in accord with MATLAB function SIZE for a matrix
      %
      %   See also.
      val = size(obj.K,varargin{:});
    end % function size

    function disp(obj)
      disp('     A nefConstFunction object F(X,U,XI,T) = K');
    end % function disp
  end % methods
end % classdef
