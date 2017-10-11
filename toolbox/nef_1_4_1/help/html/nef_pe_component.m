%% Performace evaluators component
%

%%
% As was mentioned in the estimation component section all the NEF
% estimators provide the estimates in the form of a conditional pdf of the
% state. A common task is to measure estimation error and compare
% performance of several estimators against the true value of the state.
% Obtaining such a performance index requires a procedure consisting of 
% 
% # collecting data from Monte Carlo (MC) simulations, 
% # extracting appropriate indicators from the conditional distribution of
% the state provided by individual estimators,
% # evaluating the performance index.
% 
% Within the NEF, this process is facilitated by the performance evaluator
% component. It contains methods for initialization, data processing and
% evaluation of the performance index which is selected during initialization
% as well as the number of estimators to be compared and the expected
% number of MC simulations. The performance index can be evaluated at any
% time, not necessarily after completing all the MC simulations. Currently,
% the performance evaluation component provides performance indices given
% in the following Table.

%%
% 
% <html>
% <table cellspacing="0" class="body" cellpadding="4" border="2" width="90%">
% 	<caption>Performance indices implemented in the NEF performance
% evaluation component</caption>
%   <thead>
%     <tr>
%       <th colspan="2" width="100%" align="center" bgcolor="#B2B2B2" valign="top">ABSOLUTE ERROR MEASURES</th>
%     </tr>
%   </thead>
%   <tbody>
%     <tr>
%       <td>MSEM</td>
%       <td>mean squared error matrix</td>
%     </tr>
%     <tr>
%       <td>RMSE</td>
%       <td>root mean squared error</td>
%     </tr>
%     <tr>
%       <td>AEE</td>
%       <td>average Euclidean error</td>
%     </tr>
%     <tr>
%       <td>HAE</td>
%       <td>harmonic average error</td>
%     </tr>
%     <tr>
%       <td>GAE</td>
%       <td>geometric average error</td>
%     </tr>
%     <tr>
%       <td>MEDE</td>
%       <td>median error</td>
%     </tr>
%     <tr>
%       <td>MODE</td>
%       <td>mode error</td>
%     </tr>
%   </tbody>
%   <thead align="center" valign="top">
%     <tr>
%       <th colspan="2" width="100%" align="center" bgcolor="#B2B2B2" valign="top">RELATIVE ERROR MEASURES</th>
%     </tr>
%   </thead>
%   <tbody>
%     <tr>
%       <td>RMSRE</td>
%       <td>root mean squared relative error</td>
%     </tr>
%     <tr>
%       <td>ARE</td>
%       <td>average Euclidean relative error</td>
%     </tr>
%     <tr>
%       <td>BEEQ</td>
%       <td>Bayesian estimation error quotient</td>
%     </tr>
%     <tr>
%       <td>EMER</td>
%       <td>estimation error relative to measurement
% error</td>
%     </tr>
%   </tbody>
%   <thead align="center" valign="top">
%     <tr>
%       <th colspan="2" width="100%" align="center" bgcolor="#B2B2B2" valign="top">PERFORMANCE MEASURES</th>
%     </tr>
%   </thead>
%   <tbody>
%     <tr>
%       <td>NCI</td>
%       <td>non-credibility index</td>
%     </tr>
%     <tr>
%       <td>ANEES</td>
%       <td>average normalized estimation error
% squared</td>
%     </tr>
%   </tbody>
% </table>
% </html>
% 

%%
% *Example of performance evaluator use*
%
% This example will demonstrate the use of |nefPerformaceEvaluator| class
% for determination of the root mean square error (RMSE). The RMSE will be
% compared for unscented Kalman filter and particle filter described by
% objects |UKF| and |PF|
PF = nefPF(model,'sampleSize',1000);
UKF = nefUKF(model);

%%
% Thus te RMSE wil be evaluated for two filters and 100 Monte Carlo runs
% will be processed
nFilters = 2;
MCRuns = 100;

%%
% The RMSE performace evaluator can be created as instance of the
% |nefPerformaceEvaluator| class in the following way
RMSE_PE = nefPerformanceEvaluator(model,timeHorizon,MCRuns,nFilters,'method','RMSE');

%%
% Now, the Monte Carlosimulations will be performed with processing the
% results obtained from individual estimators in performance evaluators
for i = 1:MCRuns
    [val_PF] = estimate(PF,z,u);
    [val_UKF] = estimate(UKF,z,u);
    for k = 1:timeHorizon
        Data.state{1,k} = [x(:,k); u(:,k)];
        estData{1,1,k} = val_PF{k};
        estData{2,1,k} = val_UKF{k};
    end
    processData(RMSE_PE,Data,estData);
end

%%
% Finally, the performance indices are obtained using command
rmse = performanceValue(RMSE_PE);