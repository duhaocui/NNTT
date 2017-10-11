classdef nefEstimator < handle
  %file @nefEstimator/nefEstimator.m
  % nefEstimator Properties:
  %   nefEstimator_optsDef - definition of the nefPF options
  %   sys - object defining system
  %   estQueue - internal queue for estimates
  %   time - current time instant
  %   x0 - initial condition
  %   opts - value of the nefPF options
  % nefEstimator Methods:
  %   nefEstimator - class constructor
  %   fixedLagSmoothing - Runs fixed lag smoothing process.
  %   prediction - Runs prediction process.
  %   filtering - Runs filtering process.
  %   estimate - Runs estimation process (filtering, prediction or fixedLagSmoothing).
  %   disp - display class info
  %   processOptions - Processes all options.
  %   setupOption - Setup option OPTS to given or default values.
  %   showOptions - Displays value of all options.
  %   showOption - Displays value of option given by index IDX.
  %   getOption - Gets value of option OPT.
  %   setOption - Sets option OPT to value VAL.

  % NEF version 1.4.1
  % Copyright (c) 2006 - 2017 NFT developement team,
  %              Identification and Decision Making Research Group, Department of Cybernetics,
  %              University of West Bohemia
  
  properties (Constant)
    % Definition of nefEstimator options
    % field #1: importance w.r.t visibility - show 0=never, 1=sometimes, 2=always
    % field #2: parent ('' for root or 'option.:value' for child)
    % field #3: option name
    % field #4: option default value
    % field #5: option description
    % field #6: option test
    nefEstimator_optsDef = {1,'','taskType','filtering','type of the estimation task',@(x) any(strcmpi(x,{'prediction','filtering','fixedLagSmoothing','fixedPointSmoothing','fixedIntervalSmoothing'}));
    1,'','taskPar',0,'parameter of the estimation task in effect',@(x) isscalar(x) && (mod(x,1)==0);
    1,'','dataProcessing','sequential','type of data processing',@(x) any(strcmpi(x,{'sequential','batch'}));
    0,'','disableChecks',0,'determines whether the checks in objects should be in effect',@(x)x==0 || x==1;
    0,'','verbose',0,'determines verbosity',@(x)x==0 || x==1};
  end % constant properties

  properties (SetAccess = 'protected') % protected properties
    sys; % System which state is estimated
    estQueue; % Buffer for results to remember
    time = 0; % Specification of time
    x0; % Initial condition for estiamtion
    opts=struct('Visibility',{},'Parent',{},'Name',{},'Default',{},'Current',{},'Description',{}); % Storage for estimator options
  end % properties
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = nefEstimator(system,varargin)
    % NEFESTIMATOR Creates NEFESTIMATOR object.
    %
    %   OBJ = NEFESTIMATOR(system,varargin) creates a NEFESTIMATOR object OBJ 
    %   representing particle filter for model SYSTEM.
    %
    %   the following user-defined parameters can be changed 
    %   via the standard Parameter-Value MATLAB mechanism
    %  PARAMETER              DEFAULT VAL  DESCRIPTION                  VALID VALUES
    %  =============================================================================
    %  taskType               filtering    type of task to be run       'prediction','filtering','fixedLagSmoothing','fixedPointSmoothing','fixedIntervalSmoothing'
    %  dataProcessing         sequential   way to process input data    'sequential','batch'
    %  disableChecks          0            disables checvks in the object 1,0
    %  verbose                0            verbosity                     0,1
    %
    %   See also ESTIMATE.

      p = inputParser;
      p.FunctionName = 'NEFESTIMATOR';
      p.addRequired('system',@(x) isa(x, 'nefSystem'));
      p.addParamValue('x0',[],@(x) isa(x,'nefRV'));
      for  i = 1:size(nefEstimator.nefEstimator_optsDef,1)
        p.addParamValue(nefEstimator.nefEstimator_optsDef{i,3},[],nefEstimator.nefEstimator_optsDef{i,6});
      end
      p.parse(system,varargin{:});


      obj.sys = p.Results.system;
      obj.time = 0;

      processOptions(obj,nefEstimator.nefEstimator_optsDef,p.Results,'showMessages',0)
      % chain - like disabling checks
      if getOption(obj,'disableChecks')
        disableChecks(obj.sys);
      end

      % setup x0
      % if initial condition undefined use initial condition of the system
      if isempty(p.Results.x0)
        obj.x0 = obj.sys.x0;
        if getOption(obj,'verbose')
          fprintf('%s: Initial condition (\"x0\") unspecified, using initial condition of the model: %s\n',class(obj),class(obj.sys.x0));
        end
      else
        if size(p.Results.x0,1) ~= obj.sys.dimState
          error('NEF:NEFESTIMATOR:WrongInitStateDim','INITSTATE has wrong dimension')
        end
        obj.x0 = p.Results.x0;
      end
      % setup estQueue
      switch getOption(obj,'taskType')
        case {'filtering','prediction'}
          obj.estQueue = queue(1);
        case {'fixedlagsmoothing'}
          obj.estQueue = queue(abs(getOption(obj,'taskPar')));
        case {'fixedpointsmoothing'}
          % TODO add fixedPointSmoothing
        case {'fixedintervalsmoothing'}
          % TODO add fixedIntervalSmoothing
      end
    end % nefEstimator constructor

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FUNCTION ESTIMATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val,varargout] = estimate(obj,Measurement,Input)
    % ESTIMATE Runs estimation process (filtering, prediction or fixedLagSmoothing).
    %
    %   VAL = ESTIMATE(OBJ,MEASUREMENT, INPUT) runs the estimation process 
    %   given the MEASUREMENT and INPUT. The type of process is specified on
    %   object creation. IF second output argument is specified, and taskType is 'filtering'
    %   also the predictive estimates are returned.
    %
    %   Example
    %
    %   See also.

      nout = max(nargout,1)-1;
      switch getOption(obj,'taskType')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'filtering'
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if nout == 0 % only filtering estimate is requested
            [val] = filtering(obj,Input,Measurement);
          else % filtering and predictive estimates are requested
            [val,tmp] = filtering(obj,Input,Measurement);
            varargout(1) = {tmp};
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'prediction'
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          [val] = prediction(obj,Input,Measurement);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'fixedlagsmoothing'
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          [val] = fixedLagSmoothing(obj,Input,Measurement);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'fixedpointsmoothing'
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % TODO
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'fixedintervalsmoothing'
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % TODO
      end % case
    end % function estimate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FILTERING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val,varargout] = filtering(obj,Input,Measurement)
    % FILTERING Runs filtering process.

      % expand empty Input to proper length
      if isempty(Input)
        Input = zeros(1,size(Measurement,2));
      end
      % prepare empty cell array for results
      if strcmp(getOption(obj,'dataProcessing'),'sequential')
        val = cell(1,size(Measurement,2));
      else
        val = cell(1,1);
      end
      % if also prediction estimate is requested by the second output argument
      nout = max(nargout,1)-1;
      if nout == 1
        if strcmp(getOption(obj,'dataProcessing'),'sequential')
          valPred = cell(1,size(Measurement,2));
        else
          valPred = cell(1,1);
        end
      end

      % obtain prediction estimate from queue
      predEstimate = obj.estQueue(end);
      for k = 1 : size(Measurement,2)
        % return prediction estimate if requested
        if nout == 1
          if strcmp(getOption(obj,'dataProcessing'),'sequential')
            valPred{k} = predEstimate.RV;
          else % batch dataProcessing
            if k == size(Measurement,2)
              valPred{1} = predEstimate.RV;
            end
          end
        end
        % measurement update -> p(x_time|z^time)
        filtEstimate = measurementUpdate(obj,predEstimate,Input(:,k),Measurement(:,k),obj.time);
        % time update -> p(x_time+1|z^time)
        predEstimate = timeUpdate(obj,filtEstimate,Input(:,k),obj.time);
        % return filtering estimate:
        % at each time instant for sequential dataProcessing
        % at last time instant for batch dataProcessing
        if strcmp(getOption(obj,'dataProcessing'),'sequential')
          val{k} = filtEstimate.RV;
        else % batch dataProcessing
          if k == size(Measurement,2)
            val{1} = filtEstimate.RV;
          end
        end
        % increase time
        obj.time = obj.time + 1;
      end % for
      % add the last prediction into the queue
      obj.estQueue(end+1) = predEstimate;
      % return prediction estimate if requested
      if nout == 1
        varargout(1) = {valPred};
      end
    end % function filtering

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREDICTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = prediction(obj,Input,Measurement)
    % PREDICTION Runs prediction process.
  
      % expand empty Input to proper length
      if isempty(Input)
        Input = zeros(1,size(Measurement,2)+getOption(obj,'taskPar')-1);
      end
      % prepare empty cell array for results
      % for time = 0 , return also prediction p(x_lag|z^-1)
      if strcmp(getOption(obj,'dataProcessing'),'sequential')
        if obj.time == 0
          val = cell(1,size(Measurement,2)+1);
        else
          val = cell(1,size(Measurement,2));
        end
      else % batch dataProcessing
        val = cell(1,1);
      end
      % obtain prediction estimate from queue
      predEstimate = obj.estQueue(end);
      %----------------------------------------------------------
      j = 0; % pointer to val
      % for time = 0 , return also prediction p(x_lag|z^-1)
      if obj.time == 0
        % compute and return prediction estimate p(x_lag|z^-1)
        % at each time instant for sequential dataProcessing
        % at last time instant for batch dataProcessing (in this case (k=0) if no measurement is given)
        if strcmp(getOption(obj,'dataProcessing'),'sequential')|| size(Measurement,2) == 0
          pred = predEstimate;
          for i = 1:getOption(obj,'taskPar')-1
            pred = timeUpdate(obj,pred,Input(:,i),obj.time+i-1);
          end
          j = j+1;
          val{j} = pred.RV;
        end
      end
      %----------------------------------------------------------
      for k = 1:size(Measurement,2)
        % measurementUpdate -> p(x_time|z^time)
        filtEstimate = measurementUpdate(obj,predEstimate,Input(:,k),Measurement(:,k),obj.time);
        % time Update -> p(x_time+1|z^time)
        predEstimate = timeUpdate(obj,filtEstimate,Input(:,k),obj.time);
        %----------------------------------------------------------
        % compute and return prediction estimate:
        % at each time instant for sequential dataProcessing or at last time instant  for batch dataProcessing
        if strcmp(getOption(obj,'dataProcessing'),'sequential')|| (k == size(Measurement,2))
          % compute prediction estimate p(x_time+lag|z^time)
          pred = predEstimate;
          for i = 1:getOption(obj,'taskPar')-1
            pred = timeUpdate(obj,pred,Input(:,i+k),obj.time+i+k-1);
          end
          j = j+1;
          val{j} = pred.RV;
        end
        %----------------------------------------------------------
        % increase time
        obj.time = obj.time + 1;
      end
      % add the last prediction into the queue
      obj.estQueue(end+1) = predEstimate.RV;
    end % function prediction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FLSMOOTHING Rauch Tung Striebel smoothing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = fixedLagSmoothing(obj,Input,Measurement)
    % FIXEDLAGSMOOTHING Runs fixed lag smoothing process.

      % expand empty Input to proper length
      if isempty(Input)
        Input = zeros(1,size(Measurement,2));
      end
      val = cell(1,size(Measurement,2)-max(0,abs(getOption(obj,'taskPar'))-obj.time));
      % INSERT NEW DATA TO QUEUE
      % obtain prediction estimate from queue
      p = obj.estQueue(end);
      predEstimate = p.predEstimate;
      j = 0;
      for kM = 1:size(Measurement,2)
        %==============================================================
        % measurement update -> p(x_k|z^k)
        [filtEstimate] = measurementUpdate(obj,predEstimate,Input(:,kM),Measurement(:,kM),obj.time);
        % time update -> p(x_k+1|z^k)
        [predEstimate] = timeUpdate(obj,filtEstimate,Input(:,kM),obj.time,'smoothingPurpose',1);
        %==============================================================
        % if enough data avaliable
        if obj.time >= abs(getOption(obj,'taskPar'))
          % compute and return the fixed lag estimate
          % at each time instant for sequential dataProcessing or at last time instant for batch dataProcessing
          if strcmp(getOption(obj,'dataProcessing'),'sequential')|| (kM == size(Measurement,2))
            % initial condition for smoothing
            smoothEstimate = filtEstimate;
            % going back in time
            for i = 0 : abs(getOption(obj,'taskPar'))-1
              p = obj.estQueue(end-i);   % obtain data
              [smoothEstimate]  = smoothUpdate(obj,smoothEstimate,p.filtEstimate,p.predEstimate);
            end
            j = j+1;
            val{j} = smoothEstimate.RV;
          end % if
        end % if
        % update queue
        p.filtEstimate = filtEstimate;
        p.predEstimate = predEstimate;
        p.time = obj.time;
        obj.estQueue(end+1) = p;
        % increase time
        obj.time = obj.time + 1;
      end % for
    end % function fixedLagSmoothing


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROCESS OPTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function processOptions(obj,optsDef,pResults,varargin)
    % PROCESSOPTIONS Processes all options.

      p = inputParser;
      p.FunctionName = 'PROCESSOPTIONS';
      p.addParamValue('showMessages',1,@(x)x==0 || x==1);
      p.parse(varargin{:});


      opts.Visibility = {optsDef{:,1}};
      opts.Parent = {optsDef{:,2}};
      opts.Name = {optsDef{:,3}};
      opts.Default = {optsDef{:,4}};
      opts.Current = {optsDef{:,4}};
      opts.Description = {optsDef{:,5}};

      for i = 1:length(opts.Parent)
        if strcmpi(opts.Parent{i},'') % if it is a root option
          [opts] = setupOption(obj,opts,pResults,i,p.Results.showMessages);
          % now, check for children (only for string values
          if ischar(opts.Current{i})
            % if a child exists
            children = strcmpi(opts.Parent,[opts.Name{i},':',opts.Current{i}]);
            if any(children)
              for j = 1:length(children)
                if children(j)
                  [opts] = setupOption(obj,opts,pResults,j,p.Results.showMessages);
                end
              end
            end
          end
        end
      end
      % appending options to the current ones
      if isempty(obj.opts)
        obj.opts = opts;
      else
        obj.opts.Visibility = {obj.opts.Visibility{:} opts.Visibility{:}};
        obj.opts.Parent = {obj.opts.Parent{:} opts.Parent{:}};
        obj.opts.Name = {obj.opts.Name{:} opts.Name{:}};
        obj.opts.Default = {obj.opts.Default{:} opts.Default{:}};
        obj.opts.Current = {obj.opts.Current{:} opts.Current{:}};
        obj.opts.Description = {obj.opts.Description{:} opts.Description{:}};
      end
    end % function processOptions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETUP OPTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [opts] = setupOption(obj,opts,pResults,idx,showMessages)
    % SETUPOPTION Setup option OPTS to given or default values.

      optSpecified = isfield(pResults,opts.Name{idx});
      % if the option is not specified
      if (optSpecified && isempty(pResults.(opts.Name{idx}))) || ~optSpecified
        if isa(opts.Default{idx},'function_handle')
          opts.Current{idx} = opts.Default{idx}(obj);
        elseif ischar(opts.Default{idx}) % convert strings to lowercase
          opts.Current{idx} = lower(opts.Default{idx});
        else
          opts.Current{idx} = opts.Default{idx};
        end
        if showMessages && opts.Visibility{idx} >= 1 
          tmpmsg1 = [class(obj) ': Option %s (\"%s\") unspecified, using default value: '];
          % using Current instead of Default as Default may be a function
          switch class(opts.Current{idx})
            case 'char'
              tmpmsg2 = sprintf('\"%s\"',opts.Current{idx});
            case {'single','double'}              
              if floor(opts.Current{idx}) == opts.Current{idx}
                tmpmsg2 = sprintf('\"%u\"',opts.Current{idx});
              else
                tmpmsg2 = sprintf('\"%f\"',opts.Current{idx});
              end
            otherwise
              error('unknown type')
          end
          fprintf([tmpmsg1,tmpmsg2,'\n'],opts.Description{idx},opts.Name{idx});
        end
      else % take the user-specified value
        if ischar(pResults.(opts.Name{idx})) % convert strings to lowercase
          opts.Current{idx} = lower(pResults.(opts.Name{idx}));
        else
          opts.Current{idx} = pResults.(opts.Name{idx});
        end
      end
    end % function setupOptions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SHOW OPTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function showOptions(obj)
    % SHOWOPTIONS Displays value of all options.

      for i = 1:length(obj.opts.Parent)
        if strcmpi(obj.opts.Parent{i},'') % if it is a root option
          showOption(obj,i);
          % now, check for children (only for string values
          if ischar(obj.opts.Current{i})
            % if a child exists
            children = strcmpi(obj.opts.Parent,[obj.opts.Name{i},':',obj.opts.Current{i}]);
            if any(children)
              for j = 1:length(children)
                if children(j)
                  showOption(obj,j);
                end
              end
            end
          end
        end
      end
    end % function showOptions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SHOWOPTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function showOption(obj,idx)
    % SHOWOPTION Displays value of option given by index IDX.

      if obj.opts.Visibility{idx} >= 1 
        %tmpmsg1 = '  Option %s (\"%s\") is:';
        tmpmsg1 = ['  ',class(obj),'.%s = '];
        switch class(obj.opts.Current{idx})
          case 'char'
            tmpmsg2 = sprintf('\"%s\"',obj.opts.Current{idx});
          case {'single','double'}
            if all(floor(obj.opts.Current{idx}) == obj.opts.Current{idx})
              tmpmsg2 = sprintf('\"%u\"',obj.opts.Current{idx});
            else
              tmpmsg2 = sprintf('\"%f\"',obj.opts.Current{idx});
            end
          otherwise
            error('unknown type')
        end
        fprintf([tmpmsg1,tmpmsg2,'\n'],obj.opts.Name{idx});
      end
    end % function showOption

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETOPTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val] = getOption(obj,opt)
    % GETOPTION Gets value of option OPT.

      val = obj.opts.Current{strcmpi(obj.opts.Name,opt)};
    end % function getOption

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETOPTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setOption(obj,opt,val)
    % SETOPTION Sets option OPT to value VAL.

      obj.opts.Current{strcmpi(obj.opts.Name,opt)} = val;
    end % function setOption

  end %methods

  %methods (Abstract = logical(1))
  %        function [predEst] = timeUpdate(obj,prevPred)
  %        function [filtEst] = measurementUpdate(obj,prevEst,z)
  %        function [smoothEst] = smoothingUpdate(obj,prevSmooth,prevPred,prevFilt,Input)
  %end % abstract methods

end % classdef
