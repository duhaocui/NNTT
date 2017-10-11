classdef nefCholFactorFunction < nefFunction
  %file @nefCholFactorFunction/nefCholFactorFunction.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia

  properties (SetAccess = 'protected')
    %%%%%%%%%%%%%%%%%%%%%
    % CONSTANT
    %%%%%%%%%%%%%%%%%%%%%
    S@double = [];
    P@double = [];
  end % propeties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    function [obj] = nefCholFactorFunction(varargin)
      % NEFCONSTFUNCTION Creates nefCholFactorFunction object
      %
      %   OBJ = NEFCONSTFUNCTION(P) creates a NEFCONSTFUNCTION object OBJ
      %   representing a constant function
      %       f(X,U,XI,T) = P = SS'
      %   P must be a square positive semidefinite matrix
      %
      %   OBJ = NEFCONSTFUNCTION(S) creates a NEFCONSTFUNCTION object OBJ
      %   representing a constant function
      %       f(X,U,XI,T) = SS' = P
      %   S triangular matrix     

      matrixTypeValues = {'cov','fact'};

      p = inputParser;
      p.FunctionName = 'NEFCHOLFACTORFUNCTION';
      p.addRequired('X',@isnumeric);
      p.addParamValue('matrixType',@(x) any(strcmpi(x,matrixTypeValues)))
      p.addParamValue('check',1,@(x)x==0 || x==1);
      p.parse(varargin{:});
      obj.check = p.Results.check;


      switch p.Results.matrixType
        case {'cov'}
          % check p.Results.X for positive semidefinitness
          P = p.Results.X;
          if obj.check && any(any(P ~= P'))
            error('NEF:nefCholFactorFunction:MatSymmetric','Matrix P must be symmetric.')
          end
          if obj.check && sum(eig(P) < 0)
            error('NEF:nefCholFactorFunction:MatPosSemidef','Matrix P must be positive semidefinite.')
          end
          % store nonefactored parameter
          obj.P = p.Results.X;
          % make S factor and store it
          [obj.S] = nefCholFactorFunction.CholFactor(p.Results.X);
        case {'fact'}
          S = p.Results.X;
          % check consistency of p.Results.X
          [s1,s2] = size(S);
          if obj.check && s1 ~= s2
            error('NEF:nefCholFactorFunction:MatUSquare','Matrix S must be square')
          end
          % store the CHOL factor
          obj.S = S;
          % unfactored matrix P not yet computed
          obj.P = [];
        otherwise
          disp('NEF:nefCholFactorFunction: Allowable matrix types are "cov" (for P -> S) and "fact" (for S -> P) only.')
      end % switch


      obj.dimFunction = size(obj.S);
      obj.isLinear = 1;
      obj.isAdditive = 1;
      obj.isContiuous = 1;
      obj.isDifferentiable1State = 1;
      obj.isDifferentiable2State = 1;
      obj.isDifferentiable1Noise = 1;
      obj.isDifferentiable2Noise = 1;

    end % nefCholFactorFunction constructor



    function [val] = evaluate(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATE(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=P, the function always returns P.
      %
      %   See also EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      if isempty(obj.P)
        obj.P = obj.S*obj.S';
      end
      val = obj.P;
    end % function evaluate

    function [val] = evaluateS(obj,State, Input,Noise,Time)
      % EVALUATE Evaluating function at specified point.
      %
      %   VAL = EVALUATEU(OBJ,STATE,INPUT,NOISE,TIME) returns value of
      %   the S factor of the function OBJ at the point given by STATE, INPUT, NOISE and
      %   TIME. As f(X,U,XI,T)=P, the function always returns S.
      %
      %   See also EVALUATE EVALUATED EVALDIFF1STATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      val = obj.S;
    end % function evaluateU

    function [val] = evalDiff1State(obj,State, Input,Noise,Time)
      % EVALDIFF1STATE Evaluating first derivative with respect to state
      %
      %   VAL = EVALDIFF1STATE(OBJ,STATE,INPUT,NOISE,TIME) returns value
      %   of the first derivative of the function OBJ with respect to
      %   state at the point given by STATE, INPUT, NOISE and TIME.
      %   As f(X,U,XI,T)=K, the derivative is always zero.
      %
      %   See also EVALUATE, EVALDIFF2STATE, EVALDIFF1NOISE, EVALDIFF2NOISE.
      error('NEF:nefCholFactorFunction:FirstStateDerUndefined','First derivative along State undefined')
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
      error('NEF:nefCholFactorFunction:SecondStateDerUndefined','Second derivative along State undefined')
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
      error('NEF:nefCholFactorFunction:FirstNoiseDerUndefined','First derivative along Noise undefined')
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
      error('NEF:nefCholFactorFunction:SecondNoiseDerUndefined','Second derivative along Noise undefined')
    end % function evalDiff2Noise

    function [val] = size(obj,varargin)
      % SIZE Returns size of the function.
      %
      %   VAL = SIZE(OBJ,VARARGIN) Size of K is returned. The input and output 
      %   arguments can be used in accord with MATLAB function SIZE for a matrix
      %
      %   See also.
      val = size(obj.S,varargin{:});
    end % function size

    function disp(obj)
      disp('     A nefCholFactorFunction object F(X,U,XI,T) = SS');
    end % function disp
  end % methods

  methods(Static)
    function [S]= CholFactor(P)
      %
      % Algorithm for decomposition of a symmetric positive-semidefinite 
      % matrix M. M is decomposited into products M=SS'
      % If P is positive definite, S is lower triangular
      % Input symmetric positive-semidefinite matrix M
      % Outputs matrix S
      D = det(P);

      if D == 0
        %[v,d] = eig(P);
        %S = v*diag(sqrt(diag(d)))*inv(v);
        [L,PerMat] = nefCholFactorFunction.modChol(P);
        S = (L*PerMat)';
      elseif D > 0
        S = chol(P)';
      else
        warning('NEF:nefCholFactorFunction:UnstableSysy','Possible unstability of system')
        error('NEF:nefCholFactorFunction:NegDefMat','Cannot find Cholesky factor of a negative definite matrix')
      end
    end %CholFactor

    function [L,P] = modChol(A)
      % modified Cholesky decomposition for positive semidefinite matrices
      % with pivoting
      [n,m] = size(A);
      L = A;
      idx = 1:n;

      P = eye(n); % permutation matrix
      r = zeros(1,n);
      for i = 1:n
        d = diag(L);
        [maxd,im] = max(d(i:n));
        r(i) = i+im-1;

        if (maxd == 0) % all done
          if norm(L(i+1:end,i),i)
            error('Matrix must be positive semidefinite')
          end
          break
        elseif (maxd < 0)
          error('Matrix must be positive semidefinite')
        end

        %Iperm = eye(n);
        if i ~= r(i) %interchange ith and r(i)th rows and columns
          L(:,[i r(i)]) = L(:,[r(i),i]);
          L([i r(i)],:) = L([r(i),i],:);
          idx([i, r(i)]) = idx([r(i), i]);
          %Iperm = Iperm(idx,:);
        end
        L(i,i) = sqrt(L(i,i));
        L(i,i+1:end) = L(i,i+1:end)/L(i,i);
        %    for j = i+1:n
        %      l(j,i) = (A(j,i) - sum(l(j,1:i-1).*l(i,1:i-1)))/l(i,i);
        %    end
        j = i+1:n;
        L(j,:) = L(j,:) - L(i,j)'*L(i,:);
        %L(j,j) = L(j,j) - sum(L(i,j).*L(i,j));
      end
      L = triu(L);
      P = eye(n);
      P = P(:,idx);
    end % function modChol
  end % static methods 
end % classdef
