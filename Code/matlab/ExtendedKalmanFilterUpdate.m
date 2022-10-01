function [mu, Sigma] = ExtendedKalmanFilterUpdate(mu, Sigma, u, z, gFunc, gParam, hFunc, hParam, G, H, Q, R)

%%
% 
% PURPOSE
% --------------
% A single update at time t of the traditional state space model:
% 
% z(t) = h(x(t)) + delta(t), delta(t) ~ N(0,Q_(t))             Measurement equation
% x(t) = g(u(t),x(t-1)) + eps(t), eps(t) ~ N(0,R_(t))          State equation
% where x(t) is the n-dim state, u(t) is the n-dim control, z(t) is the k-dim observed data.
% 
%
% CALL
% --------------
% [mu, Sigma] = KalmanFilterUpdate(mu, Sigma, u, z, gFunc, gParam, hFunc, hParam, Q, R)
%        
% INPUTS
% --------------
% mu                n-by-1          mu(t-1), the posterior mean of the state at time t-1                 
% Sigma             n-by-n          Sigma(t-1), the posterior covariance matrix of the state at time t-1
% u                 n-by-1          u(t), the controls at time t
% z                 k-by-1          z(t), the measurments at time t
% gFunc             function        State transition function
% gParam            cell array      Other arguments to gFunc
% hFunc             function        Measurement function
% hParam            cell array      Other arguments to hFunc
% Q, R                              Model parameters, see the state-space model above.
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
[muBar, G] = feval(gFunc, u, mu, gParam); % Computing g() and the matrix of derivatives at mu
SigmaBar = G*Sigma*G' + R;

%% Measurement update - updating the N(muBar, SigmaBar) prior with the new data point
[hBar, H] = feval(hFunc, muBar, hParam); % Evaluating the measurement function and derivatives at muBar
K = SigmaBar*H' / (H*SigmaBar*H' + Q); % Kalman Gain
mu = muBar + K*(z - hBar);
Sigma = (eye(n)- K*H)*SigmaBar;

