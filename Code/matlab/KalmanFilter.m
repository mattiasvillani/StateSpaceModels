function [MU, SIGMA] = KalmanFilter(U, Z, A, B, C, Q, R, mu0, Sigma0)

%%
% 
% PURPOSE
% --------------
% Running a full Kalman filtering round for the traditional state space model:
% 
% z(t) = C*x(t) + delta(t), delta(t) ~ N(0,Q_(t))             Measurement equation
% x(t) = A*x(t-1) + B*u(t) + eps(t), eps(t) ~ N(0,R_(t))      State equation
% where x(t) is the n-dim state, u(t) is the n-dim control, z(t) is the k-dim observed data.
% 
%
% CALL
% --------------
% [mu, Sigma] = KalmanFilter(mu, Sigma, u, z, A, B, C, Q, R)
%        
% INPUTS
% --------------
% U                 T-by-m          U(t,:) is u(t)'             
% Z                 T-by-k          Z(t,:) is z(t)'
% A, B, C, Q, R                     Model parameters, see the state-space model above.
% mu0               n-by-1          Initial mean of state
% Sigma0            n-by-n          Initial covariance of state
%
% OUTPUTS
% ---------------
% MU                T-by-n          MU(t,:) is the posterior mean of the state at time t                 
% SIGMA             n-by-n-by-T     SIGMA(:,:,t) is the posterior covariance matrix of the state at time t
%
% AUTHOR
% ---------------
% Mattias Villani, Link?ping University. e-mail: mattias.villani@gmail.com
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
T = size(Z,1); 
n = length(mu0);

MU = zeros(T,n);
SIGMA = zeros(n,n,T);

% The Kalman iterations
mu = mu0;
Sigma = Sigma0;
for t = 1:T
    [mu, Sigma] = KalmanFilterUpdate(mu, Sigma, U(t,:)', Z(t,:)', A, B, C, Q, R);
    MU(t,:) = mu;
    SIGMA(:,:,t) = Sigma;
end
