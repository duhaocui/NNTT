classdef nefFunction < handle
  %file @nefFunction/nefFunction.m
  
  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIMENSIONS OF THE FUNCTION PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dimState@double = 0;
    dimInput@double = 0;
    dimNoise@double = 0;
    dimTime@double = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIMENSIONS OF THE FUNCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dimFunction@double = 0;
    %%%%%%%%%%%%%%%%%%%%%
    % FUNCTION PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%
    isAdditive@double = 0;
    isContiuous@double = 0;
    isLinear@double = 0;
    isDifferentiable1State@double = 1;
    isDifferentiable2State@double = 1;
    isDifferentiable1Noise@double = 1;
    isDifferentiable2Noise@double = 1;
    isConst@double = 0;
    isDiff1NoiseConst@double = 0;
    %%%%%%%%%%%%%%%%%%%%%
    % PERFORM CHECKS
    %%%%%%%%%%%%%%%%%%%%%
    check@double = 1;
  end % propeties

  methods
    function [obj] = nefFunction()
      % NEFFUNCTION Creates nefFunction object
      %
      %   OBJ = NEFFUNCTION() creates a NEFFUNCTION object OBJ
      %
      %   See also NEFCONSTFUNCTION, NEFLINEARFUNCTION, NEFHANDLEFUNCTION.
    end % nefFunction constructor

    function [val] = evaluate(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefFunction:EvaluateUndefined','Function evaluation undefined')
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
          error('NEF:nefFunction:TooManyVariables','Only single values of Input, Noise and Time Allowed');
        end
      end
      DataCount = size(State,2);
      if obj.dimFunction(2) == 1 % vector function
        val = zeros([obj.dimFunction(1) DataCount]);
        for i = 1:DataCount
          val(:,i) = evaluate(obj,State(:,i),Input,Noise,Time);
        end
      else % matrix function
        val = zeros([obj.dimFunction DataCount]);
        for i = 1:DataCount
          val(:,:,i) = evaluate(obj,State(:,i),Input,Noise,Time);
        end
      end
    end % function evaluateOverAllStates

    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefFunction:FirstStateDerUndefined','First derivative along State undefined')
    end % function evalDiff1State

    function [val] = evalDiff2State(obj,State, Input,Noise,Time)
      % EVALDIFF2STATE Evaluating second derivative with respect to state
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefFunction:SecondStateDerUndefined','Second derivative along State undefined')
    end % function evalDiff2State

    function [val] = evalDiff1Noise(obj,State, Input,Noise,Time)
      % EVALDIFF1NOISE Evaluating first derivative with respect to noise
      %
      %   VAL = EVALDIFF1NOISEE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF2NOISE.
      error('NEF:nefFunction:FirstNoiseDerUndefined','First derivative along Noise undefined')
    end % function evalDiff1Noise

    function [val] = evalDiff2Noise(obj,State, Input,Noise,Time)
      % EVALDIFF2NOISE Evaluating second derivative with respect to noise
      %
      %   VAL = EVALDIFF2STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the second derivative of the function OBJ with respect to
      %   noise at the point given by STATE, INPUT, NOISE and TIME
      %
      %   See also EVALUATE, EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE.
      error('NEF:nefFunction:SecondNoiseDerUndefined','Second derivative along Noise undefined')
    end % function evalDiff2Noise

    function [val] = size()
    end % function size


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECKDIMENSIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkDimensions(obj,State,Input,Noise,Time)
      if obj.check 
        if (obj.dimState > 0) && (size(State,1) ~= obj.dimState)
          error('NEF:nefFunction:StateDim','Wrong dimension of state')
        end
        if (obj.dimInput > 0) && (size(Input,1) ~= obj.dimInput)
          error('NEF:nefFunction:InputDim','Wrong dimension of input')
        end
        if (obj.dimNoise > 0) && (size(Noise,1) ~= obj.dimNoise)
          error('NEF:nefFunction:NoiseDim','Wrong dimension of noise')
        end
        if (obj.dimTime > 0) && (size(Time,1) ~= obj.dimTime)
          error('NEF:nefFunction:TimeDim','Wrong dimension of time')
        end
      end
    end % function checkDimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECKDATACOUNT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkDataCount(obj,State,Input,Noise,Time)
      % CHECKDATACOUNT Checks data count of variables provided
      %

      % get dataCount to be checked
      dataCount = ([size(State,2) size(Input,2) size(Noise,2) size(Time,2)]);
      if any(nonzeros(dataCount)~=dataCount(1)) % if all non-empty variables do not have same dataCount
        error('NEF:nefFunction:InconsStateInputNoise','Inconsistent dataCount of state, input, noise and time')
      end
    end %function checkDataCount

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disp(obj)
      fprintf('A nefFunction object\n');
    end % function disp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISABLECHECKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disableChecks(obj)
      obj.check = 0;
    end % function disableChecks
  end % methods
end % classdef
