classdef nefLinFunction < nefFunction
  %file @nefLinFunction/nefLinFunction.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % MATRICES
    %%%%%%%%%%%%%%%%%%%%%
    F@double = [];
    G@double = [];
    H@double = [];
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NEFLINEFUNCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefLinFunction(f,g,h,varargin)
      % NEFLINEFUNCTION Creates nefLinFunction object
      %
      %   OBJ = NEFLINEFUNCTION(f,g,h) creates a NEFLINEFUNCTION object OBJ
      %   representing a linear function
      %       f(X,U,XI,T) = F*X + G*U + H*XI
      %   F,G,H must be matrices of proper dimensions or empty matrices
      %
      %   See also NEFFUNCTION, NEFCONSTFUNCTION, NEFHANDLEFUNCTION.
      p = inputParser;
      p.FunctionName = 'NEFLINEFUNCTION';
      p.addRequired('f',@(x)isnumeric(x) || isempty(x));
      p.addRequired('g',@(x)isnumeric(x) || isempty(x));
      p.addRequired('h',@(x)isnumeric(x) || isempty(x));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(f,g,h, varargin{:});
      obj.check = p.Results.check;


      % CHECKING DIMENSIONS OF THE MATRICES
      nef = size(f,1);
      ng = size(g,1);
      nh = size(h,1);
      % IF SIZE(*,1) is ZERO, THE FUNCTION
      % DOES NOT DEPEND ON THE CORRESPONDING VARIABLE

      % CHECK NON-ZERO DIMS FOR CONSISTENCY
      dims = [];
      if obj.check && nef > 0, dims(end+1) = nef;end
      if obj.check && ng > 0, dims(end+1) = ng;end
      if obj.check && nh > 0, dims(end+1) = nh;end

      % AT LEAST ONE MATRIX MUST BE SPECIFIED
      if obj.check && isempty(dims)
        error('NEF:nefLinFunction:NoMatrixSpecified','At least one matrix must be specified')
      end

      % CHECK CONSISTENCY
      if obj.check && any(mean(dims) ~= dims)
        error('NEF:nefLinFunction:MatrixDim','Inconsistent matrix dimensions')
      end

      obj.F = f;
      obj.G = g;
      obj.H = h;
      obj.dimState = size(f,2);
      obj.dimInput = size(g,2);
      obj.dimNoise = size(h,2);
      obj.isLinear = 1;
      obj.isAdditive = 1;
      obj.isContiuous = 1;
      obj.isDifferentiable1State = 1;
      obj.isDifferentiable2State = 1;
      obj.isDifferentiable1Noise = 1;
      obj.isDifferentiable2Noise = 1;
      obj.isDiff1NoiseConst = 1;
      % dimFunction is given by the maximum of nef, ng, nh (for the case that some of them are zero)
      obj.dimFunction = [max([nef ng nh]),1];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evaluate(obj,State,Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. VAL =  f(X,U,XI,T) = F*X + G*U + H*XI.
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end
      %%%%%%%%%%%%%%%%%%%%%%%%
      % ACTUAL EVALUATION
      %%%%%%%%%%%%%%%%%%%%%%%%
      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      val = zeros(obj.dimFunction(1),DataCount);
      if obj.dimState > 0,val = val + obj.F*State;end
      if obj.dimInput > 0,val = val + obj.G*Input;end
      if obj.dimNoise > 0,val = val + obj.H*Noise;end

    end % function evaluate
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FUNCTION EVALUATION OVER ALL STATES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evaluateOverAllStates(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point over all states provided.
      %
      %   VAL = EVALUATEOVERALLSTATES(OBJ,STATE,INPUT,NOISE,TIME) returns values of
      %   the function OBJ at the points given by STATE (all values), INPUT, NOISE and
      %   TIME.
      %
      %   See also EVALUATE EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        if max([size(Input,2) size(Noise,2) size(Time,2)]) > 1;
          error('NEF:nefLinFunction:TooManyVariables','Only single values of Input, Noise and Time Allowed');
        end
      end
      DataCount = size(State,2);
      val = obj.F*State;
      if obj.dimInput > 0,val = bsxfun(@plus,val,obj.G*Input);end
      if obj.dimNoise > 0,val = bsxfun(@plus,val,obj.H*Noise);end
    end % function evaluateOverAllStates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE 1ST DERIVATIVE
    % WITH RESPECT TO STATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME.
      %   As f(X,U,XI,T)= F*X + G*U + H*XI, the derivative is always F.
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end
      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);

      val = obj.F(:,:,ones(1,DataCount));
    end % function evalDiff1State

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE 2ND DERIVATIVE
    % WITH RESPECT TO STATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalDiff2State(obj,State, Input,Noise,Time)
      % EVALDIFF2STATE Evaluating second derivative with respect to state
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)= F*X + G*U + H*XI, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      val = zeros([obj.dimState obj.dimState obj.dimFunction]);
    end % function evalDiff2State

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE 1ST DERIVATIVE
    % WITH RESPECT TO NOISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalDiff1Noise(obj,State, Input,Noise,Time)
      % EVALDIFF1NOISE Evaluating first derivative with respect to noise
      %
      %   VAL = EVALDIFF1NOISEE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)= F*X + G*U + H*XI, the derivative is always H.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      val = obj.H(:,:,ones(1,DataCount));
    end % function evalDiff1Noise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE 2ND DERIVATIVE
    % WITH RESPECT TO NOISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalDiff2Noise(obj,State, Input,Noise,Time)
      % EVALDIFF2NOISE Evaluating second derivative with respect to noise
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %   As f(X,U,XI,T)= F*X + G*U + H*XI, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      val = zeros([obj.dimNoise obj.dimNoise obj.dimFunction DataCount]);
    end % function evalDiff2Noise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = size(obj,varargin)
      % SIZE Returns size of the function.
      %
      %   VAL = SIZE(OBJ,VARARGIN) Size of function value is returned. The input and output
      %   arguments can be used in accord with MATLAB function SIZE for a matrix
      %
      %   See also.
      val = size(zeros(obj.dimFunction),varargin{:});
    end % function size

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      disp('     A nefLinFunction object F(X,U,XI,T) = F*X + G*U + H*XI');
    end % function disp

  end % methods
end % classdef
