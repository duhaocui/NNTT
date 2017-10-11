classdef queue
%file @queue/queue.m
% QUEUE - defines a queue of a specified length
%
% PROPERTIES:
%   LEN - length of the queue.
%   POS - current position in queue.
%   Q   - a cell array of queue elements.
%
% METHODS:
%   QUEUE - class constructor.
%   END   - returns an object at actual position in QUEUE.
%   SUBSREF - implements a special subscripted reference.
%   SUBSASGN - implements a special subscripted assignment.
%

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

properties (SetAccess = 'protected')
        %%%%%%%%%%%%%%%%%%%%%
        % LENGTH OF QUEUE
        %%%%%%%%%%%%%%%%%%%%%
        len@double = 0;
        %%%%%%%%%%%%%%%%%%%%%
        % POINTER OF QEUE
        %%%%%%%%%%%%%%%%%%%%%
        pos@double = 0;
end % protected properties
properties
        %%%%%%%%%%%%%%%%%%%%%
        % queue
        %%%%%%%%%%%%%%%%%%%%%
        q@cell = {};
end %properties
methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUEUE CONSTRUCTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = queue(len)
        % QUEUE Creates a queue of a specified length.
        %
        %   OBJ = QUEUE(LEN) creates a QUEUE object OBJ of length LEN.
        %
        %   See also.
        obj.q = cell(1,len);
        obj.len = len;
        end % queue constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [val] = end(obj,k,n)
        % Returns an object at actual position in QUEUE.
        val = obj.pos;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSREF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function val = subsref(obj,s)
        % Implements a special subscripted reference.
        switch s.type
        case '()'
            val = obj.q{mod(s.subs{:}-1,obj.len)+1};
        end
        end % function  subsref

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSASGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function  [obj] = subsasgn(obj,s,element)
        % Implements a special subscripted assignment.
        switch s.type
        case '()'
            if s.subs{:} == obj.pos+1 % if index is end+1,
                obj.pos = mod(s.subs{:}-1,obj.len)+1; % increase pos
            end
            obj.q{mod(s.subs{:}-1,obj.len)+1} = element; % assign the element
        end
        end % function subsasgn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function disp(obj)
        % DISP Display object.
        obj.q
    end % function disp

end % methods
end % classdef
