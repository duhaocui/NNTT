%% Estimators component
%
% The estimation component of the NEF provides implementation of a number
% of estimators used to solve the estimation problem where the main task is
% to find a conditional probability density function of the state
% conditioned by available measurements. The estimate of the state $x_k$ is
% given by the posterior pdf $p(x_k|z^\ell,u^\ell)$, where $z^\ell$ is the
% sequence of measurements up to time instant $\ell$, i.e., $z^\ell=[z_0^T,
% z_1^T,\dots, z_\ell^T]^T$.
%
% The estimation problem itself can be divided according to
% the relation of $k$ and $\ell$ into the following three special cases.
%
% * If $\ell=k$, the problem is called filtering.
% * If $\ell<k$, the problem is called prediction.
% * If $\ell>k$, the problem is called smoothing.
%
% Solution to all of the three problems is provided by the Bayesian
% functional relations (BFR) that are analytically tractable only for a few
% special models. Nevertheless, most of the state estimators can be seen as
% an approximative solution to the BFR and therefore the NEF estimation
% component follows the BFR idea.
%

%% Implemented estimators
%
% The NEF provides many well known local and global filters (see the table
% below). Most of them provide besides "canonival" implementation also
% numerically stable version. Notable is also the fact that the Gaussian
% sum filter can utilize any of the local filters (i.e. it is not
% restricted to use of the extended Kalman filter). 
%
% <html>
% <table cellspacing="0" class="body" cellpadding="4" border="2">
%   <caption>Estimators implemented in the NEF estimation
% component</caption>
%   <thead>
%     <tr bgcolor="#B2B2B2" valign="top">
%       <th>NEF class</th>
%       <th>estimators</th>
%     </tr>
%   </thead>
%   <tbody>
%     <tr>
%       <td>nefKalman,
% nefSKalman
% nefUDKalman</td>
%       <td>(extended) Kalman filter (standard,
% square-root and UD versions)</td>
%     </tr>
%     <tr>
%       <td>nefDD1,
% nefSDD1,
% nefDD2, nefSDD2</td>
%       <td>central difference Kalman filter, divided
% difference filter (1st and 2nd
% order) (standard and square-root
% version)</td>
%     </tr>
%     <tr>
%       <td>nefUKF,
% nefSUKF</td>
%       <td>unscented Kalman filter (standard
% and square-root version), cubature
% Kalman filter</td>
%     </tr>
%     <tr>
%       <td>nefItKalman</td>
%       <td>iterated Kalman filter based on any
% of the above local filter</td>
%     </tr>
%     <tr>
%       <td>nefGSM</td>
%       <td>Gaussian sum filter based on any of
% the above local filter</td>
%     </tr>
%     <tr>
%       <td>nefEnKF</td>
%       <td>ensemble Kalman filter</td>
%     </tr>
%     <tr>
%       <td>nefPF</td>
%       <td>bootstrap filter, generic particle filter,
% auxiliary particle filter, unscented
% particle filter</td>
%     </tr>
%     <tr>
%       <td>nefSIF</td>
%       <td>stochastic integration filter
%       </td>
%     </tr>
%   </tbody>
% </table> 
% </html>

%%
% All the estimator implement filtering and prediction task and most of
% them provide also the ability to perform the smoothing.
%
% <html>
% <table cellspacing="0" class="body" cellpadding="4" border="2" width="80%">
%   <caption>Estimation tasks supported by individual estimators</caption>
%   <thead>
%     <tr bgcolor="#B2B2B2" valign="top">
%       <th width="40%">Estimator</th>
%       <th width="20%">filtering</th>
%       <th width="20%">prediction</th>
%       <th width="20%">smoothing</th>
%     </tr>
%   </thead>
%   <tbody>
%     <tr>
%       <td>nefKalman</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefSKalman</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefUDKalman</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td></td>
%     </tr>
%     <tr>
%       <td>nefItKalman</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefDD1</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefSDD1</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefDD2</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefSDD2</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefUKF</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefSUKF</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%     </tr>
%     <tr>
%       <td>nefGSM</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td></td>
%     </tr>
%     <tr>
%       <td>nefEnKalman</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td></td>
%     </tr>
%     <tr>
%       <td>nefPF</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td></td>
%     </tr>
%     <tr>
%       <td>nefSIF</td>
%       <td align="center">x</td>
%       <td align="center">x</td>
%       <td></td>
%     </tr>
%   </tbody>
% </table>
% </html>
% 

%% Estimation experiment setup
% 
% The estimation experiment setup is quite simple. As first it is neccessary
% to create instance of any of the estimator class. The class constructor
% always expects object describing the model as first input argument.
% Aditional argumets which are used to change the default behaviour of the
% estimators age appended as touples of parameter name and patameter value.
% Upon object creation the class contructor reposts the current initial
% settings of the estimator.
%
% For example an (extended) Kalman smoother can be created issuing the
% following command
 ekf_smoother = nefKalman(model,'taskType','fixedLagSmoothing','taskPar',2);
 
 %%
 % The estimation process can be then carried out in to way. Either
 % emplying the |estimate| method which is an universal interface for any
 % estimation task or using the individual methods |timeUpdate|,
 % |measurementUpdate| and |smoothUpdate|. 
 %
 % The first method is the prefered one as it executes the update method in
 % the right order and passes to then the appropriate data automatically.
 % The second is useful e.g. in case of tight coupling with an controller
 % where it is necesary to "insert" the evaluation of the control law
 % between the |measurementUpdate| and |timeUpdate| tasks. 
 %
 % The following line represents typical command for execution of
 % estimation experiment 
 est_pdfs = estimate(ekf_smoother,z,u);
 
 %%
 % The results of the estimation proces are stored as cell arrays of pdf's.
 % The point estimates can be then determined using appropriate method of
 % the random variable class. For example to determine the mean value of
 % the smoothing pdf at time instant |k| one uses the method of the random
 % variable |evalMean| in following way
 xMean = evalMean(est_pdfs{k});