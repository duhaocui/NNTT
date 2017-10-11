classdef nefIdentification
  % file @nefIdentificaton/nefIdentificaton.m

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  % Identification and Decision Making Research Group, Department of Cybernetics,
  %  University of West Bohemia
    
  properties (SetAccess = 'protected')
  end 
  
  methods(Static)
  function [L] = gainCalculation(f,h,w,v) 
        % Calculating of kalman gain L
        %
        % [L] = gainCalculation(f,h,w,v)
        %
        % f - must be object of nefLinFunction, state equation
        % h - must be object of nefLinFunction, measurement equation
        % w - objekt nefGaussianRV, state noises 
        % v - objekt nefGaussianRV, measurement noises 
        if(~isa(f,'nefLinFunction'))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'f is not objekt of nefLinFunction'])   
        end
        if(~isa(h,'nefLinFunction'))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'h is not objekt of nefLinFunction'])   
        end
        if(~isa(w,'nefGaussianRV'))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'w is not objekt of nefGaussianRV'])   
        end
        if(~isa(v,'nefGaussianRV'))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'c is not objekt of nefGaussianRV'])   
        end
        if(~prod(size(w)==size(f)))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'Q must be matrix relevant dimension'])   
        end
        if(~prod(size(v)==size(h)))
           error(['NEF:nefIdentificaton:calculationL - ',...
           'R must be matrix relevant dimension'])   
        end
        P = dare(f.F', h.F', f.H*w.Var*f.H', h.H*v.Var*h.H');
        L = P*h.F'*inv(h.F*P*h.F'+h.H*v.Var*h.H');     
    end
  end
end