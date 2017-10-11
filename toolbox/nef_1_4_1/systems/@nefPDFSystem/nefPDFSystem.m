classdef nefPDFSystem < nefSystem
  %file @nefPDFSystem/nefPDFSystem.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TWO RV'S DEFINING THE SYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x; % state
    z; % measurement
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefPDFSystem(trPDF,measPDF,initPDF)
      p = inputParser;
      p.KeepUnmatched = true;
      p.FunctionName = 'NEFPDFSYSTEM';
      p.addRequired('trPDF',@(x) isa(x, 'nefRV'));
      p.addRequired('measPDF',@(x) isa(x, 'nefRV'));
      p.addRequired('initPDF',@(x) isa(x, 'nefRV'));
      p.parse(trPDF,measPDF,initPDF);

      obj@nefSystem(p.Unmatched) % call nefSystem constructor

      if obj.check
        % initPDF has to be independent of state, input and time
        if (initPDF.dimState > 0)
          error('NEF:nefPDFSystem:InitPDFState','INITPDF must not be parametrized by state')
        end
        if (initPDF.dimInput > 0)
          error('NEF:nefPDFSystem:InitPDFInput','INITPDF must not be parametrized by input')
        end
        if (initPDF.dimTime > 0)
          error('NEF:nefPDFSystem:InitPDFTime','INITPDF must not be parametrized by time')
        end

        % dim of initPDF = dim of trPDF
        if (size(initPDF,1) ~= trPDF.dimRV)
          error('NEF:nefPDFSystem:InconsTrPDFInitPDF','Inconsistent TRPDF and INITPDF')
        end

        % dim of trPDF = dim of State in trPDF
        if (trPDF.dimState > 0)
          if trPDF.dimRV ~= trPDF.dimState
            error('NEF:nefPDFSystem:InvalidTrPDF','TRPDF is not valid transient pdf')
          end
        end

        % dim of trPDF = dim of state in measPDF
        if (measPDF.dimState > 0)
          if trPDF.dimRV ~= measPDF.dimState
            error('NEF:nefPDFSystem:InconsTrPDFMeasPDF','Inconsistent TRPDF and MEASPDF')
          end
        end
        
        % permit control either only in trPDF or measPDF
        if ((trPDF.dimInput > 0) && (measPDF.dimInput >0))
            error('NEF:nefPDFSystem:TwoInputs','The control input is not permited to be present in both TRPDF and MEASPDF')
        end
      end % if obj.check

      obj.x = trPDF;
      obj.z = measPDF;
      obj.x0 = initPDF;
      obj.dimState = trPDF.dimState;
      obj.dimInput = max([trPDF.dimInput,measPDF.dimInput]);
      obj.dimTime = max([trPDF.dimTime,measPDF.dimTime]);
      obj.dimMeasurement = size(measPDF,1);
      obj.state= [];
    end % nefPDFSystem constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE SYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [z,x] = simulate(obj,steps,Input,varargin)  
      % SIMULATE Simulates the system for given number of steps.
      %
      %   [Z,X] = SIMULATE(OBJ,STEPS,INPUT) simulates the system OBJ for
      %   STEPS time instants with input given by INPUT. If the system OBJ
      %   does not have specified the previous state value then the initial
      %   state is generated from the initial state PDF. Measurement is
      %   stored in Z and state in X.
      %
      %   [Z,X,OBJ] = SIMULATE(OBJ,STEPS,INPUT,VARARGIN) simulates the
      %   system OBJ for STEPS time instants with input given by INPUT. The
      %   VARARGIN specifies additional parameters using  standard
      %   parameter-value MATLAB mechanism, i.e by specifying tuples of
      %   quantityName and quantity. The currently only valid value of the
      %   quantityName is initialState. The corresponding quantity value of
      %   the  fixed initial state. Measurement is stored in Z and state in
      %   X. The OBJ output parameter is necessary in case the user wants
      %   to continue with the simulation in later point (so the last state
      %   value is preserved).
      %
      %   Example
      %    [z,x] = simulate(model,20,u)
      %    [z,x,model] = simulate(model,20,u,'initialState',x0)
      %
      
      p = inputParser;
      p.FunctionName = 'SIMULATE';
      p.addRequired('steps',@(x) isnumeric(x) | (x > 0));
      p.addRequired('Input',@(x) isnumeric(x) | isempty(x));
      p.addParamValue('initialState',[],@(x) isfloat(x));
      p.parse(steps,Input,varargin{:});

      steps = p.Results.steps;
      Input = p.Results.Input;
    
      ux0 = p.Results.initialState; % initial state supplied by user
      
      if obj.check && (steps <= 0)
        error('NEF:nefPDFSystem:NegativeSteps','STEPS must be greater than zero')
      end
      if isempty(Input)
        INPUT = zeros(1,steps);
      else
        if obj.check && (size(Input,2) ~= steps)
          error('NEF:nefPDFSystem:InconsInputSteps','Inconsistent INPUT and STEPS')
        end
        INPUT = Input;
      end
      
      if ~isempty(ux0)
          if (size(ux0,1) ~= obj.dimState)
              error('NEF:nefPDFSystem:InconsX0State','Inconsistent size of X0 and STATE')
          end
      end
      
      x = zeros(obj.dimState,steps);
      z = zeros(obj.dimMeasurement,steps);
      for k = 1:steps
        % evaluate new state
        if isempty(obj.state) % no previous state present
            % determine the initial state - either it is given by the user
            % or generated from pdf
            if isempty(ux0)
                x(:,k) = drawSample(obj.x0,1,[],[],[],obj.time);
            else
                x(:,k) = ux0(:,k);
            end
        else % we are alredy at later time instants
            x(:,k) = drawSample(obj.x,1,obj.state,INPUT(:,k),[],obj.time);
        end
        obj.state = x(:,k);
        
        % current state and corresponding measurement
        z(:,k) = drawSample(obj.z,1,x(:,k),INPUT(:,k),[],obj.time);
        % increment the time
        obj.time = obj.time + 1;
      end
    end % function simulate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISABLECHECKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function disableChecks(obj)
      obj.check = 0;

      disableChecks(obj.x);
      disableChecks(obj.z);
    end % function disableChecks
  end %methods
end % the class
