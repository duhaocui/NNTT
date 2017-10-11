classdef nefEqSystem < nefSystem
  %file @nefEqSystem/nefEqSystem.m
  % NEFEQSYSTEM - defines a NEFEQSYSTEM OBJECT of a specified DIMENSION, STATE and
  % NOISE
  %
  % PROPERTIES:
  %   F    - state function
  %   H    - measurememt function
  %   W    - state noise
  %   V    - measurememt noise
  %   X0   - initial state
  %
  % METHODS:
  %   NEFEQSYSTEM - class constructor
  %   SIMULATE - simulates the system for specified number of steps
  %   INIT     - changes initial state
  %

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % nefEqSystem PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    f % state function
    h % measurement function
    w % state noise
    v % measurement noise
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % nefEqSystem LIKELIHOOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    log_likelihood_handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % nefEqSystem TRANSPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    log_transitionPDF_handle
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NEFEQSYSTEM CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefEqSystem(f,h,w,v,x0,varargin)
      % NEFSYSTEM Creates a nefEqSystem object.
      %
      %   S = NEFEQSYSTEM(F,H,W,V,X0) creates a NEFEQSYSTEM object S given 
      %   by neffunctions F and H, random variables W, V, state X0.
      %
      %   See also NEFPDFSYSTEM.
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFEQSYSTEM';
      p.addRequired('f',@(x) isa(x, 'nefFunction'));
      p.addRequired('h',@(x) isa(x, 'nefFunction'));
      p.addRequired('w',@(x) isa(x, 'nefRV'));
      p.addRequired('v',@(x) isa(x, 'nefRV'));
      p.addRequired('x0',@(x) isa(x, 'nefRV'));
      p.addParamValue('logLikelihood',[],@(x) isa(x, 'function_handle') && (nargin(x) == 4));
      p.addParamValue('logTransitionPDF',[],@(x) isa(x, 'function_handle') && (nargin(x) == 4));
      p.parse(f,h,w,v,x0,varargin{:});

      obj@nefSystem(p.Unmatched) % call nefSystem constructor

      % check functions f and h
      if obj.check 
        if (size(f,2) ~= 1)
          error('NEF:nefEqSystem:InconsistentFDimension','F must return a column vector')
        end
        if (size(h,2) ~= 1)
          error('NEF:nefEqSystem:InconsistentHDimension','H must return  a column vector')
        end

        if (size(f,1) ~= f.dimState)
          error('NEF:nefEqSystem:InconsistentFStateDimension','Dimension of STATE in F must be same as dimension of F ')
        end
        if (f.dimState ~= h.dimState)
          error('NEF:nefEqSystem:InconsistentStateDimension','Dimension of STATE in F and H must be same')
        end
        % permit control either only in trPDF or measPDF
        if ((f.dimInput > 0) && (h.dimInput > 0))
          error('NEF:nefEqSystem:TwoInputs','The control input is not permited to be present in both transient and measurement function')
        end

        if (f.dimNoise ~= f.dimState)
          error('NEF:nefEqSystem:InconsistentFNoiseAndFStateDimension','Dimension of NOISE and STATE in F must be same')
        end

        if (h.dimNoise ~= size(h,1))
          error('NEF:nefEqSystem:InconsistentDimension','Dimension of NOISE in H must be same as dimension of H')
        end

        if (w.dimRV ~= f.dimNoise)
          error('NEF:nefEqSystem:InconsistentDimension','Dimension of W must be same as dimension of NOISE in F')
        end

        if (v.dimRV ~= h.dimNoise)
          error('NEF:nefEqSystem:InconsistentDimension','Dimension of V must be same as dimension of NOISE in H')
        end

        if (size(x0,1) ~= f.dimState)
          error('NEF:nefEqSystem:InconsistentDimension','Dimension of X0 must be same as dimension of STATE in F')
        end

        % check initial condition
        % x0 has to be independent of state, input and time
        if (x0.dimState > 0)
          error('NEF:nefEqSystem:X0State','x0 must not be parametrized by state')
        end
        if (x0.dimInput > 0)
          error('NEF:nefEqSystem:x0Input','x0 must not be parametrized by input')
        end
        if (x0.dimTime > 0)
          error('NEF:nefEqSystem:x0Time','x0 must not be parametrized by time')
        end

        % check measurement noise w
        % w has to be independent of state, input and time
        if (w.dimState > 0)
          error('NEF:nefEqSystem:WState','W must not be parametrized by state')
        end

        % check measurement noise v
        % v has to be independent of state, input and time
        if (v.dimState > 0)
          error('NEF:nefEqSystem:WState','V must not be parametrized by state')
        end
      end % if obj.check

      obj.time = 0;
      obj.isLinear = (f.isLinear & h.isLinear);
      obj.isAdditive = (f.isAdditive & h.isAdditive);
      obj.isGaussian = ((isa(w,'nefGaussianRV')) & (isa(v,'nefGaussianRV')) & (isa(x0,'nefGaussianRV')));
      obj.f = f;
      obj.h = h;
      obj.w = w;
      obj.v = v;
      obj.x0 = x0;
      obj.log_likelihood_handle = p.Results.logLikelihood;
      obj.log_transitionPDF_handle = p.Results.logTransitionPDF;
      obj.dimState = max([f.dimState,h.dimState]);
      obj.dimMeasurement = size(h,1);
      obj.dimInput = max([f.dimInput,h.dimInput]);
      obj.dimTime = max([f.dimTime,h.dimTime]);

      obj.state = [];

    end % nefEqSystem constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [z,x,obj] = simulate(obj,steps,Input,varargin)
      % SIMULATE Simulates the system for given number of steps.
      %
      %   [Z,X] = SIMULATE(OBJ,STEPS,INPUT) simulates the system OBJ for
      %   STEPS time instants with input given by INPUT. If the system OBJ
      %   does not have specified the previous state value then the initial
      %   state is generated from the initial state PDF. The state and
      %   measurement noises are drawn from the corresponding PDFs.
      %   Measurement is stored in Z and state in X.
      %
      %   [Z,X,OBJ] = SIMULATE(OBJ,STEPS,INPUT,VARARGIN) simulates the
      %   system OBJ for STEPS time instants with input given by INPUT. The
      %   VARARGIN specifies additional parameters using  standard
      %   parameter-value MATLAB mechanism, i.e by specifying tuples of
      %   quantityName and quantity. The valid values of the quantityName
      %   are: stateNoise, measurementNoise and initialState. The
      %   corresponding quantity values determine the state noise,
      %   measurement noise and fixed initial state. Measurement is stored
      %   in Z and state in X. The OBJ output parameter is necessary in
      %   case the user wants to continue with the simulation in later
      %   point (so the last state value is preserved).
      %
      %   Example
      %    [z,x,model] = simulate(model,20,u,'stateNoise',w,'measurementNoise',v)
      %    [z,x,model] = simulate(model,1,[],'initialState',x0)
      %
      %   See also.

      p = inputParser;
      p.FunctionName = 'SIMULATE';
      p.addRequired('steps',@(x) isnumeric(x) | (x > 0));
      p.addRequired('Input',@(x) isnumeric(x) | isempty(x));
      p.addParamValue('w',[],@(x) isfloat(x));
      p.addParamValue('stateNoise',[],@(x) isfloat(x));      
      p.addParamValue('v',[],@(x) isfloat(x));
      p.addParamValue('measurementNoise',[],@(x) isfloat(x));
      p.addParamValue('x0',[],@(x) isfloat(x));
      p.addParamValue('initialState',[],@(x) isfloat(x));
      p.parse(steps,Input,varargin{:});

      steps = p.Results.steps;
      Input = p.Results.Input;
    
      uw = p.Results.stateNoise; % state noise supplied by user
      uv = p.Results.measurementNoise; % measurement noise supplied by user
      ux0 = p.Results.initialState; % initial state supplied by user
      
      if ~isempty(p.Results.w)
          warning('NEF:nefEqSystem:deprecatedSpecification','nefEqSystem: Deprecated specification of state noise parameter - use ''stateNoise'' instead.')
          if ~isempty(uw)
              warning('NEF:nefEqSystem:doubleSpecification','nefEqSystem: State noise specified twice. Only the ''stateNoise'' value is used.')
          else
              uw = p.Results.w; % state noise supplied by user
          end
      end
      
      if ~isempty(p.Results.v)
          warning('NEF:nefEqSystem:deprecatedSpecification','nefEqSystem: Deprecated specification of measurement noise parameter - use ''measurementNoise'' instead.')
          if ~isempty(uv)
              warning('NEF:nefEqSystem:doubleSpecification','nefEqSystem: Measurement noise specified twice. Only the ''measurementNoise'' value is used.')
          else
              uv = p.Results.v; % measurement noise supplied by user
          end
      end
      
      if isempty(Input)
        INPUT = zeros(1,steps);
      else
        if obj.check && (size(Input,2) ~= steps)
          error('NEF:nefEqSystem:InconsInputSteps','Inconsistent INPUT and STEPS')
        end
        INPUT = Input;
      end

      if obj.check
        if ~isempty(uw)
          if (size(uw,2) ~= steps)
            error('NEF:nefEqSystem:InconsWSteps','Inconsistent W and STEPS')
          end
          if (size(uw,1) ~= obj.dimState)
            error('NEF:nefEqSystem:InconsWState','Inconsistent size of W and STATE')
          end
        end

        if ~isempty(uv)
          if (size(uv,2) ~= steps)
            error('NEF:nefEqSystem:InconsVSteps','Inconsistent V and STEPS')
          end
          if (size(uv,1) ~= obj.dimMeasurement)
            error('NEF:nefEqSystem:InconsVState','Inconsistent size of V and MEASUREMENT')
          end
        end
        
        if ~isempty(ux0)
          if (size(ux0,1) ~= obj.dimState)
            error('NEF:nefEqSystem:InconsX0State','Inconsistent size of X0 and STATE')
          end
        end
      end % if obj.check

      x = zeros(obj.dimState,steps);
      z = zeros(obj.dimMeasurement,steps);
      for i = 1: steps
          % evaluate the state
          if isempty(obj.state) % no previous state present
              % determine the initial state - either it is given by the user
              % or generated from pdf
              if isempty(ux0)
                  x(:,i) = drawSample(obj.x0,1);
              else
                  x(:,i) = ux0(:,i);
              end
          else % we are alredy at later time instants
              if isempty(uw)
                  w = drawSample(obj.w, 1, obj.state, INPUT(:,i), [], obj.time);
              else
                  w = uw(:,i);
              end
              x(:,i) = evaluate(obj.f, obj.state, INPUT(:,i), w, obj.time);
          end          
          obj.state = x(:,i);

          % evaluate the measurement
          if isempty(uv)
              v = drawSample(obj.v, 1, x(:,i), INPUT(:,i), [], obj.time);
          else
              v = uv(:,i);
          end
          z(:,i) = evaluate(obj.h,x(:,i), INPUT(:,i), v, obj.time);
          
        % increment the time
        obj.time = obj.time + 1;
      end
    end % function simulate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOGLIKELIHOOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = logLikelihood(obj,Measurement,State,Input,Time)
      % LOGLIKELIHOOD Evaluates loglikelihood of the system.
      %
      %   [val] = LOGLIKELIHOOD(OBJ,MEASUREMENT,STATE,INPUT,TIME) evaluates system OBJ loglikelihood
      %   for input given by INPUT, state by STATE, measurement by
      %   MEASUREMENT. and time by TIME. The values of
      %   loglikelihoods is stored in VAL.
      %
      %   Example
      %
      %   See also.
      if obj.check && isempty(obj.log_likelihood_handle)
        error('NEF:nefEqSystem:UnspecifiedLogLikelihood','LogLikelihood is not specified (see nefEqSystem constructor help)')
      else
        val = obj.log_likelihood_handle(Measurement,State,Input,Time);
      end % if
    end % funtion logLikelihood

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOGTRANSITIONPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = logTransitionPDF(obj,newState,State,Input,Time)
      % LOGTRANSITIONPDF Evaluates logarithm of transitionPDF of the system.
      %
      %   [val] = LOGTRANSPDF(OBJ,NEWSTATE,STATE,INPUT,TIME) evaluates natural 
      %   logarithm of system OBJ transitionPDF
      %   for input given by INPUT, state by STATE, measurement by
      %   MEASUREMENT. and time by TIME. The values of
      %   likelihoods is stored in VAL.
      %
      %   Example
      %
      %   See also.
      if obj.check && isempty(obj.log_transitionPDF_handle)
        error('NEF:nefEqSystem:UnspecifiedlogTransitionPDF','logTransitionPDF is not specified (see nefEqSystem constructor help)')
      else
        val = obj.log_transitionPDF_handle(newState,State,Input,Time);
      end % if
    end % funtion logTransitionPDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISABLECHECKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disableChecks(obj)
      obj.check = 0;

      disableChecks(obj.f);
      disableChecks(obj.h);
      disableChecks(obj.w);
      disableChecks(obj.v);
    end % function disableChecks

  end %methods
end % the class
