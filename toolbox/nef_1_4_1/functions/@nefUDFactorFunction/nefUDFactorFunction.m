classdef nefUDFactorFunction < nefFunction
  %file @nefUDFactorFunction/nefUDFactorFunction.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % CONSTANT
    %%%%%%%%%%%%%%%%%%%%%
    U@double = [];
    D@double = [];
    K@double = [];
  end % propeties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    function [obj] = nefUDFactorFunction(varargin)
      % NEFCONSTFUNCTION Creates nefUDFactorFunction object
      %
      %   OBJ = NEFCONSTFUNCTION(K) creates a NEFCONSTFUNCTION object OBJ
      %   representing a constant function
      %       f(X,U,XI,T) = K = UDU'
      %   K must be a square positive semidefinite matrix
      %
      %   OBJ = NEFCONSTFUNCTION(U,D) creates a NEFCONSTFUNCTION object OBJ
      %   representing a constant function
      %       f(X,U,XI,T) = UDU' = K
      %   U must be a square lower triangular matrix
      %   D is a vector of diagonal elements
      %
      %   See also NEFFUNCTION, NEFLINEFUNCTION, NEFHANDLEFUNCTION.

      p = inputParser;
      p.FunctionName = 'NEFUDFACTORFUNCTION';
      p.addRequired('X',@isnumeric);
      p.addOptional('Y',[], @isnumeric);
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(varargin{:});
      obj.check = p.Results.check;

      % if only one numerical parameter is given
      if isempty(p.Results.Y) 
        % check p.Results.X for positive semidefinitness
        K = p.Results.X;
        if obj.check && any(any(K ~= K'))
          error('NEF:nefUDFactorFunction:MatSymmetric','Matrix must be symmetric')
        end
        if obj.check && sum(eig(K) < 0)
          error('NEF:nefUDFactorFunction:MatPosSemidef','Matrix K must be positive semidefinite')
        end
        % store nonefactored parameter
        obj.K = p.Results.X;
        % make UD factor and store it
        [obj.U,obj.D] = nefUDFactorFunction.UDFactor(p.Results.X);
        % if two numerical parameters are given
      else
        U = p.Results.X;
        D = p.Results.Y;
        % check p.Results.X for triangularity
        if obj.check && any(any(U ~= triu(U)))
          error('NEF:nefUDFactorFunction:MatTriL','Matrix U must be upper triangular')
        end
        if obj.check && any(diag(U) ~= 1)
          error('NEF:nefUDFactorFunction:MatTriL1','Matrix U must have ones on diagonal')
        end
        % check p.Results.Y for size
        if obj.check && ~isvector(D)
          error('NEF:nefUDFactorFunction:MatDiag','D must be a vector')
        end
        % check consistency of p.Results.X and p.Results.Y
        [u1,u2] = size(U);
        [d1,d2] = size(D);
        if obj.check && u1 ~= u2
          error('NEF:nefUDFactorFunction:MatUSquare','Matrix U must be square')
        end
        if obj.check && d2 ~= 1
          error('NEF:nefUDFactorFunction:VectDColumn','Vector D must be a column')
        end
        if obj.check && u1 ~= d1
          error('NEF:nefUDFactorFunction:UDInconsistency','Matrix U and vector D are inconsistent')
        end
        % store the UD factor
        obj.U = U;
        obj.D = D;
        % unefactored matrix not yet computed
        obj.K = [];
      end
      obj.dimFunction = size(obj.U);
      obj.isLinear = 1;
      obj.isAdditive = 1;
      obj.isContiuous = 1;
      obj.isDifferentiable1State = 1;
      obj.isDifferentiable2State = 1;
      obj.isDifferentiable1Noise = 1;
      obj.isDifferentiable2Noise = 1;

    end % nefUDFactorFunction constructor


    function [val] = evaluate(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=K, the function always returns K.
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if isempty(obj.K)
        obj.K = obj.U*diag(obj.D)*obj.U';
      end
      val = obj.K;
    end % function evaluate

    function [val] = evaluateU(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATEU(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the U factor of the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=K, the function always returns U.
      %
      %   See also EVALUATE EVALUATED EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      val = obj.U;
    end % function evaluateU

    function [val] = evaluateD(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATED(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the D factor of the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=K, the function always returns D.
      %
      %   See also EVALUATE EVALUATEU EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      val = obj.D;
    end % function evaluateD

    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME.
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefUDFactorFunction:FirstStateDerUndefined','First derivative along State undefined')
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
      error('NEF:nefUDFactorFunction:SecondStateDerUndefined','Second derivative along State undefined')
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
      error('NEF:nefUDFactorFunction:FirstNoiseDerUndefined','First derivative along Noise undefined')
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
      error('NEF:nefUDFactorFunction:SecondNoiseDerUndefined','Second derivative along Noise undefined')
    end % function evalDiff2Noise

    function [val] = size(obj,varargin)
      % SIZE Returns size of the function.
      %
      %   VAL = SIZE(OBJ,VARARGIN) Size of K is returned. The input and output 
      %   arguments can be used in accord with MATLAB function SIZE for a matrix
      %
      %   See also.
      val = size(obj.U,varargin{:});
    end % function size

    function disp(obj)
      disp('     A nefUDFactorFunction object F(X,U,XI,T) = UDU');
    end % function disp
  end %methods

  methods(Static)
    function [U,D]= UDFactor(M)
      % Algorithm for decomposition of a symmetric positive-define 
      % matrix M. M is decomposited into products M=UDU' such that 
      % U is unit upper triangular and D is diagonal.
      % Input symetric positive-define matrix M
      % Outputs matrix U vector of diagonal elements D


      m=size(M,1);
      D = zeros(m,1);
      U = zeros(m,m);

      for j=m:-1:1
        for i=j:-1:1
          sigma=M(i,j);
          for k=j+1:m
            sigma=sigma-U(i,k)*D(k)*U(j,k);
          end
          if i==j
            D(j)=sigma; 
            U(j,j)=1; 
          else
            U(i,j)=sigma/D(j); 
          end
        end
      end
    end %UDFactor
  end % static methods
end % classdef
