function [mu, Sigma] = KalmanFilterUpdate(mu, Sigma, u, z, A, B, C, Q, R)

%%
% 
% PURPOSE
% --------------
% A single update at time t of the traditional state space model:
% 
% z(t) = C*x(t) + delta(t), delta(t) ~ N(0,Q_(t))             Measurement equation
% x(t) = A*x(t-1) + B*u(t) + eps(t), eps(t) ~ N(0,R_(t))      State equation
% where x(t) is the n-dim state, u(t) is the n-dim control, z(t) is the k-dim observed data.
% 
%
% CALL
% --------------
% [mu, Sigma] = KalmanFilterUpdate(mu, Sigma, u, z, A, B, C, Q, R)
%        
% INPUTS
% --------------
% mu                n-by-1          mu(t-1), the posterior mean of the state at time t-1                 
% Sigma             n-by-n          Sigma(t-1), the posterior covariance matrix of the state at time t-1
% u                 n-by-1          u(t), the controls at time t
% z                 k-by-1          z(t), the measurments at time t
% A, B, C, Q, R                     Model parameters, see the state-space model above.
%
% OUTPUTS
% ---------------
% mu                n-by-1          mu(t), the posterior mean of the state at time t                 
% Sigma             n-by-n          Sigma(t), the posterior covariance matrix of the state at time t
%
% AUTHOR
% ---------------
% Mattias Villani, Linkoping University. e-mail: mattias.villani@gmail.com
%
% VERSION DATING
% ---------------
% FIRST     2015-07-29
% CURRENT   2015-07-29
%
% REFERENCES
% ---------------
% Thrun, Burgard and Fox (2006). Probabilistic Robotics, Algorithm Kalman_filter in Table 3.1.

%% Prelims
n = length(mu);

%% Prediction step - just moving the state on step forward without taking new measurement into account
muBar = A*mu + B*u;
SigmaBar = A*Sigma*A' + R;

%% Measurement update - updating the N(muBar, SigmaBar) prior with the new data point
K = SigmaBar*C' / (C*SigmaBar*C' + Q); % Kalman Gain
mu = muBar + K*(z - C*muBar);
Sigma = (eye(n)- K*C)*SigmaBar;

