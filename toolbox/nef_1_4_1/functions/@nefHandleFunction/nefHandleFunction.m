classdef nefHandleFunction < nefFunction
  %file @nefHandleFunction/nefHandleFunction.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % FUNCTION HANDLE
    %%%%%%%%%%%%%%%%%%%%%
    fHandle@function_handle;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROPERTIES LIST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PropertiesName= {'Diff1State','Diff2State','Diff1Noise','Diff2Noise'};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ASSIGNABLE PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Diff1State %@function_handle;
    Diff2State %@function_handle;
    Diff1Noise %@function_handle;
    Diff2Noise %@function_handle;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INDICATOR OF VECTORIZED FUNCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isVectorized@double = 0;
  end % propeties

  methods 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefHandleFunction(fun,sizes,varargin)
      % NEFHANDLEFUNCTION Creates nefHandleFunction object
      %
      %   OBJ = NEFHANDLEFUNCTION(FUN,SIZES) creates a NEFHANDLEFUNCTION object OBJ
      %   where the function is given by function_handle FUN and dimensions 
      %   of the state, input, noise and time are given in array SIZES. The function 
      %   FUN must have 4 arguments (STATE, INPUT, NOISE, TIME) even if they are not utilized.
      %
      %   The parameters FUN and SIZES can be followed by parameter/value pairs 
      %   to specify 1st or 2nd derivatives of the function with respect to either 
      %   state or noise. For example, NEFHANDLEFUNCTION(@(x,u,w,t) x+w,[1 0 1 0],'Diff1State',@(x,u,w,t) 1)
      %   will create function f(x,u,w,t) = x+w with the first derivative with respect to state 
      %   df(x,u,w,t)/dx = 1
      %
      %   Example:
      %   obj = nefHandleFunction(@(x,u,w,t) 0.95*x+sin(t*u)+w,[1 1 1 1])
      %
      %   See also NEFFUNCTION, NEFLINEARFUNCTION, NEFCONSTFUNCTION.

      p = inputParser;
      p.FunctionName = 'NEFHANDLEFUNCTION';
      p.addRequired('fun',@(x) isa(x, 'function_handle') && nargin(x) == 4);
      p.addRequired('sizes',@(x) isvector(x) && length(x) == 4 && all(x>=0) && all(mod(x,1)==0));
      p.addParamValue('Diff1State',[],@(x) isa(x, 'function_handle') && nargin(x) == 4);
      p.addParamValue('Diff2State',[],@(x) isa(x, 'function_handle') && nargin(x) == 4);
      p.addParamValue('Diff1Noise',[],@(x) isa(x, 'function_handle') && nargin(x) == 4);
      p.addParamValue('Diff2Noise',[],@(x) isa(x, 'function_handle') && nargin(x) == 4);
      p.addParamValue('isAdditive',0,@(x) isnumeric(x));
      p.addParamValue('isLinear',0,@(x) isnumeric(x));
      p.addParamValue('isVectorized',0,@(x) isnumeric(x));
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(fun,sizes, varargin{:});
      obj.check = p.Results.check;


      obj.fHandle = fun;

      obj.isAdditive = p.Results.isAdditive;
      obj.isLinear = p.Results.isLinear;
      obj.isVectorized = p.Results.isVectorized;

      obj.Diff1State = p.Results.Diff1State;
      obj.Diff2State = p.Results.Diff2State;
      obj.Diff1Noise = p.Results.Diff1Noise;
      obj.Diff2Noise = p.Results.Diff2Noise;

      obj.dimState = sizes(1);
      obj.dimInput = sizes(2);
      obj.dimNoise = sizes(3);
      obj.dimTime = sizes(4);
      % guessing!!! size
      % use a dummy zero variable for guessing size if its dimension is
      % zero
      obj.dimFunction = size(obj.fHandle(ones(max(obj.dimState,1),1),ones(max(obj.dimInput,1),1),ones(max(obj.dimNoise,1),1),ones(max(obj.dimTime,1),1)));
    end % nefHandleFunction constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FUNCTION EVALUATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evaluate(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME.
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      if obj.isVectorized
        val = obj.fHandle(State,Input,Noise,Time);
      else
        % to speed up the computation, the empty variables are expanded to zeros (i.e. can be referenced in the cycle)
        [State,Input,Noise,Time,DataCount] = nefHandleFunction.expandEmptyVariables(State,Input,Noise,Time);

        if obj.dimFunction(2) == 1 % vector function
          val = zeros(obj.dimFunction(1),DataCount);
          for i = 1:DataCount
            val(:,i) = obj.fHandle(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
          end
        else % matrix function
          val = zeros([obj.dimFunction DataCount]);
          for i = 1:DataCount
            val(:,:,i) = obj.fHandle(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
          end
        end
      end
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
          error('NEF:nefHandleFunction:TooManyVariables','Only single values of Input, Noise and Time Allowed');
        end
      end
      if obj.dimFunction(2) == 1 % vector function
        DataCount = size(State,2);
        if obj.isVectorized
          if obj.dimInput > 0, Input = Input(:,ones(1,DataCount));, end
          if obj.dimNoise > 0, Noise = Noise(:,ones(1,DataCount));, end
          if obj.dimTime > 0, time = Time(:,ones(1,DataCount));, end
          val = obj.fHandle(State,Input,Noise,Time);
        else
          val = zeros([obj.dimFunction(1) DataCount]);
          for i = 1:DataCount
            val(:,i) = obj.fHandle(State(:,i),Input,Noise,Time);
          end
        end
      else % matrix function
        DataCount = size(State,2);
        val = zeros([obj.dimFunction DataCount]);
        for i = 1:DataCount
          val(:,:,i) = obj.fHandle(State(:,i),Input,Noise,Time);
        end
      end
    end % function evaluateOverAllStates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATE !ST DERIVATIVE
    % WITH RESPECT TO STATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME.
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check
        if isempty(obj.Diff1State)
          error('NEF:nefHandleFunction:Diff1StateUndef','First derivative with respect to state undefined.')
        end
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      if obj.dimFunction(2) == 1 % vector function
        [State,Input,Noise,Time,DataCount] = nefHandleFunction.expandEmptyVariables(State,Input,Noise,Time);
        val = zeros([obj.dimFunction(1) obj.dimState DataCount]);
        for i = 1:DataCount
          val(:,:,i) = obj.Diff1State(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
        end
      else % matrix function
        error('NEF:nefHandleFunction:MFDiff1StateUndef','First derivative with respect to state undefined for Matrix Functions.')
      end
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
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if obj.check && isempty(obj.Diff2State)
        error('NEF:nefHandleFunction:Diff2StateUndef','Second derivative with respect to state undefined.')
      end
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      if obj.dimFunction(2) == 1 % vector function
        [State,Input,Noise,Time,DataCount] = nefHandleFunction.expandEmptyVariables(State,Input,Noise,Time);
        val = zeros([obj.dimFunction(1) obj.dimState obj.dimState DataCount]);
        for i = 1:DataCount
          val(:,:,:,i) = obj.Diff2State(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
        end
      else % matrix function
        error('NEF:nefHandleFunction:MFDiff2StateUndef','Second derivative with respect to state undefined for Matrix Functions.')
      end

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
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF2NOISE.
      if obj.check && isempty(obj.Diff1Noise)
        error('NEF:nefHandleFunction:Diff1NoiseUndef','First derivative with respect to noise undefined.')
      end
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      if obj.dimFunction(2) == 1 % vector function
        [State,Input,Noise,Time,DataCount] = nefHandleFunction.expandEmptyVariables(State,Input,Noise,Time);
        val = zeros([obj.dimFunction(1) obj.dimNoise DataCount]);
        for i = 1:DataCount
          val(:,:,i) = obj.Diff1Noise(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
        end
      else % matrix function
        error('NEF:nefHandleFunction:MFDiff1NoiseUndef','First derivative with respect to noise undefined for Matrix Functions.')
      end
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
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE.
      if obj.check && isempty(obj.Diff2Noise)
        error('NEF:nefHandleFunction:Diff2NoiseUndef','Second derivative with respect to noise undefined.')
      end
      if obj.check
        checkDimensions(obj,State,Input,Noise,Time);
        checkDataCount(obj,State,Input,Noise,Time);
      end

      if obj.dimFunction(2) == 1 % vector function
        [State,Input,Noise,Time,DataCount] = nefHandleFunction.expandEmptyVariables(State,Input,Noise,Time);
        val = zeros([obj.dimFunction(1) obj.dimNoise obj.dimNoise DataCount]);
        for i = 1:DataCount
          val(:,:,:,i) = obj.Diff2Noise(State(:,i),Input(:,i),Noise(:,i),Time(:,i));
        end
      else % matrix function
        error('NEF:nefHandleFunction:MFDiff2NoiseUndef','Second derivative with respect to noise undefined for Matrix Functions.')
      end
    end % function evalDiff2Noise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE OF THE FUNCTION
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
      %disp('     A nefHandleFunction object');
      disp(['     ',func2str(obj.fHandle)]);
    end % function disp
  end % methods

  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expand empty variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [State,Input,Noise,Time,DataCount] = expandEmptyVariables(State,Input,Noise,Time);
      DataCount = max([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      if isempty(State), State = zeros(1,DataCount);end
      if isempty(Input), Input = zeros(1,DataCount);end
      if isempty(Noise), Noise = zeros(1,DataCount);end
      if isempty(Time), Time = zeros(1,DataCount);end
    end % function expandEmptyVariables

  end % static methods
end % classdef
