function pdfTimePlot(rvs, varargin)
% PDFTIMEPLOT(M,V,'VALCOUNT',X,'RANGE',[A B],'TRAJECTORY',x)
% PLOTS A TIME DEVELOPMENT OF PDF OF A GAUSSIAN SCALAR RANDOM VARIABLE
%
%   Required Parameters
%   rvs - cell array of random variable (nefRV) objects
%
%   Optional parameters
%   'dimension' 'valCount' 'range' 'trajectory' 'fill'
%   'dimension' - select dimension of rvs to plot
%   'valCount' - number of point of each PDF (default: 1000)
%   'range' -  1x2 vector of plot range (default: computed from means and variances)
%   'trajectory' - 1xn vector of trajectory to be plotted
%   'fill' - scalar determining whether the plot is to be filled (fill~=0) or not (fill=0)

% NEF version 1.4.1
% Copyright (c) 2006 - 2017 NFT developement team,
%              Identification and Decision Making Research Group, Department of Cybernetics,
%              University of West Bohemia

p = inputParser;
p.FunctionName = 'pdfTimePlot';
p.addRequired('rvs', @(x)iscell(x));
p.addParamValue('index',1,@(x)isscalar(x) & x>0);
p.addParamValue('valCount',1000, @(x)isnumeric(x) && x>1);
p.addParamValue('range',[],@(x)isvector(x) & length(x)==2);
p.addParamValue('trajectory',[], @(x)isvector(x));
p.addParamValue('fill',1, @(x)isscalar(x));
p.parse(rvs, varargin{:});

TI = length(p.Results.rvs); % time instants

for i = 1:length(p.Results.rvs)
  if ~isa(p.Results.rvs{i},'nefRV')
    error('RVS must be cell array of random variables')
  end
end

if ~isempty(p.Results.trajectory) & (TI ~= length(p.Results.trajectory))
  error('RVS and TRAJECTORY must be of same size')
end

% assign variables
rvs = p.Results.rvs;

idx = p.Results.index;

valCount = p.Results.valCount;

if isempty(p.Results.range)
  beta = 3; % multiple of deviation
  for ti = 1:TI
    M = evalMean(rvs{ti});
    m(ti) = M(idx);
    V = evalVariance(rvs{ti});
    v(ti) = V(idx,idx);
  end
  lo = min(m-beta*sqrt(v));
  hi = max(m+beta*sqrt(v));
else
  if p.Results.range(1) ==p.Results.range(2)
    error('Values in range are equal')
  end
  lo = min(p.Results.range);
  hi = max(p.Results.range);
end

trajectory = p.Results.trajectory;


AZIM = 46;
ELEV = 70;

delta = (hi-lo)/(valCount-1);
points = repmat([lo:delta:hi]',1,TI);
pdfs = zeros(valCount,TI);

for ti = 1:TI
  pdfs(:,ti) = evalMarginalPDF(rvs{ti},points(:,ti)',idx);
end
% prepare figure
hold on;
grid on;
view(AZIM, ELEV);

plotcolor = 'black';
fillcolor = [0.5 1 0.2];
for ti = 1:TI;
  if p.Results.fill
    fill3(points(:,ti),repmat(ti,valCount,1), pdfs(:,ti) , fillcolor);
  else
    plot3(points(:,ti),repmat(ti,valCount,1), pdfs(:,ti) ,plotcolor, 'LineWidth',1)
  end
end
plot3(m,[1:TI],zeros(1,TI),'gs','LineWidth',2,...
  'MarkerEdgeColor','k',...
  'MarkerFaceColor','g',...
  'MarkerSize',10)
if ~isempty(trajectory)
  plot3(trajectory,[1:TI],zeros(1,TI),'--rs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
end
hold off
zlabel('pdf')
ylabel('time')
