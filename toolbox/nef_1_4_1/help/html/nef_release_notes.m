%% Nonlinear Estimation Framework Release Notes
%
%% Version 1.3.0, 14-Oct-2013
% * implemented Stochastic Integration Filter (nefSIF)
% * the whole toolbox was more optimized for speed bringing about up 2x speedup
% * corrected issue with new MATLAB version regarding cell array with one zero dimension
% * it is possible to use square-root version of the local filters as local filters in nefGSM
% * needs at least MATLAB R2009b
%
%% Version 1.2.0, 19-Mar-2010
% * added two new estimators nefSDD2 (numerically stable version of nefDD2)
% and nefEnKF (ensemble Kalman filter)
% * nefUKF, nefSUKF, nefDD1, nefDD2, nefSDD1, nefSDD2 support adaptive
% choice of parameter
% * nefGSM and nefItKalman can use all the local filters as base
% * added help into class files for all the estimators
% * added first version of documentation to the documentation browser
% * corrected scripts addNEFPath and rmNEFPath to work flawlessly on all
% platforms supported by MATLAB

%% 
%% Version 1.1.0, 19-Mar-2010
% * added performance evaluators (for example of their use see example13)
% * estimator constructors reports use of default values
%
%% Version 1.0.2 19-Mar-2010
%
% * removed cause of possible MATLAB crash after calling just clear all between two examples
% * |nefGSM| estimator can now work with |nefPDFSystem|
%
%% Version 1.0.1, 12-Mar-2010
%
% * stabilized |nefPF| algorithm (scaling weights for proper normalization)
% * added |evalLogPDF| to |nefGaussianRV|
% * replaced the |likelihood| and |transPDF| by |logLikelihood| and |logTransPDF| in |nefEqSystem|
% * added some checks of the noises to |nefEqSystem|
% * removed critical bug in PF preventing resampling from being executed in some cases
% * added optimalized code for additive systems in all derivative free filters
% * removed bug in |nefUDKalman| when the prediction was not properly updated for all measurement elements
% * |nefUD| completely rewritten and renamed to |nefUDKalman|
% * added |nefSKalman| implementing Square-root Kalman filter
% * |nefItFilter| and |nefGSM| can use any of the local filters
%
%% Version 1.0.0, 14-Jul-2009
% 
% * initial public release
% * stabilized API of the new framework
% 
% Copyright 2007-2010 NFT developement team, Identification and Decision Making Research Group, Department of Cybernetics, University of West Bohemia   