classdef nefSystem < handle
  %file @nefSystem/nefSystem.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia
  
  properties (SetAccess = 'protected') % protected properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATE AND TIME OF SYSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    state = [];
    time = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OTHER SYSTEM PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isLinear = 0;
    isGaussian = 0;
    isAdditive = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % X0 SPECIFICATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIMENSIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dimState@double = 0;
    dimMeasurement@double = 0;
    dimInput@double = 0;
    dimTime@double = 0;
    %%%%%%%%%%%%%%%%%%%%%
    % PERFORM CHECKS
    %%%%%%%%%%%%%%%%%%%%%
    check@double = 1;
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefSystem(varargin)
      p = inputParser;
      p.FunctionName = 'NEFSYSTEM';
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse( varargin{:});
      obj.check = p.Results.check;
    end % nefSystem constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE YSTEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [z,x] = simulate(obj,steps,Input)
    end % function simulate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = init(obj,x0)
      % INIT Changes initial state of the system and sets time to zero.
      %
      %   VAL = INIT(OBJ,X0) changes initial state of the system OBJ
      %   to XO and sets time to zero. X0 may be either a column
      %   vector of proper dimenstion or a nefRV object of proper
      %   dimensions. If X0 is a nefRV object, the initial state is
      %   set to a sample of the random variable.
      %
      %   VAL = INIT(OBJ) reinitializes the system OBJ by drawing the
      %   initial state according to the parameter X0 of the object. 
      %   fThe time is set to zero.
      %
      %
      %   See also SIMULATE.

      % INIT is method for changing initial state of NEFSYSTEM OBJECT
      if nargin == 2
        if ~isa(x0,'nefRV') && ~isfloat(newx0)
          error('NEF:nefSystem:X0','X0 must be a nefRV class or a numerical vector');
        end
        if (obj.dimState ~= size(x0,1))
          error('NEF:nefSystem:X0Dimension','Wrong dimension of X0');
        end
        if isfloat(x0)
          obj.state = x0;
        else
          obj.state = drawSample(x0,1);
        end
      else % no new x0 specified
        obj.state = drawSample(obj.x0,1);
      end
      obj.time = 0;
    end % funtion init
  end %methods
end % the class
