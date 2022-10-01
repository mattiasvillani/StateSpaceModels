function [mu, Sigma, N, Fx] = EKFSLAMUpdate(mu, Sigma, u, z, N, Q, R, Fx)

%%
% 
% PURPOSE
% --------------
% A single update at time t of the Extended Kalman Filter SLAM with unknown correspondances
% 
%
% CALL
% --------------
% [mu, Sigma] = EKFSLAMUpdate(mu, Sigma, u, z, Q, R)
%        
% INPUTS
% --------------
% mu                n-by-1          mu(t-1), the posterior mean of the state at time t-1                 
% Sigma             n-by-n          Sigma(t-1), the posterior covariance matrix of the state at time t-1
% u                 n-by-1          u(t), the controls at time t
% Z                 k-by-3          The measured features (r,phi,s) at time t. Z(k,:) is the range-direction-signature 
%                                   triple for the k:th feature
% N                 Scalar          Number of sighted landmarks at time t-1
% Q                 k-by-k          Covariance matrix of the observation noise, Q = diag(sigmaR,sigmaPhi, sigmaS) for
%                                   range sensor model (r = range, Phi = direction, S = feature signature).
% R                 n-by-n          Covariance matrix of state innovations
% Fx                3-by-N          Matrix with ones and zeros to reduce the state vector to only the three pose vars.
%
% OUTPUTS
% ---------------
% mu                n-by-1          mu(t), the posterior mean of the state at time t                 
% Sigma             n-by-n          Sigma(t), the posterior covariance matrix of the state at time t
% N                 scalar          Number of sighted landmarks at time t
%
% AUTHOR
% ---------------
% Mattias Villani, Linkoping University. 
% mattias.villani@gmail.com
% http://mattiasvillani.com
%
% VERSION DATING
% ---------------
% FIRST     2016-12-30
% CURRENT   2016-12-30
%
% REFERENCES
% ---------------
% Thrun, Burgard and Fox (2006). Probabilistic Robotics, Algorithm Kalman_filter in Table 3.1.

%% Prelims
n = length(mu);
Nt = N;
deltaT = 1; % TODO: make this an input
theta = mu(3);
k = size(Z,1); % Number of observed features

%% Prediction step - just moving the state on step forward without taking new measurement into account
muBar = mu + Fx'*[(u(1)/u(2))*(-sin(theta) + sin(theta + u(2)*deltaT)) ;
                  (u(1)/u(2))*(cos(theta) - cos(theta + u(2)*deltaT))  ; 
                  u(2)*deltaT];

G = eye(3+3*Nt) + Fx'*[zeros(3,2) [(u(1)/u(2))*(-cos(theta) + cos(theta + u(2)*deltaT)); (u(1)/u(2))*(-sin(theta) + sin(theta + u(2)*deltaT)) ; 0] ]*Fx;
M = diag([alpha(1)*u(1)^2 + alpha(2)*u(2)^2 , alpha(3)*u(1)^2 + alpha(4)*u(2)^2]); % Motion noise
V = [(-sin(theta) + sin(theta + u(2)*deltaT))/u(2)    u(1)*(sin(theta) - sin(theta + u(2)*deltaT))/u(2)^2 + u(1)*(cos(theta + u(2)*deltaT)*deltaT))/u(2) ;
     (cos(theta) - cos(theta + u(2)*deltaT))/u(2)    -u(1)*(cos(theta) - cos(theta + u(2)*deltaT))/u(2)^2 + u(1)*(sin(theta + u(2)*deltaT)*deltaT))/u(2);
     0                                                deltaT];
R = V*M*V'; % Motion noise in state space
SigmaBar = G*Sigma*G' + Fx'*R*Fx;

%% Measurement update - updating the N(muBar, SigmaBar) prior with the new data point
for i = 1:k
    z = Z(k,:);
    muBarNew
end

[hBar, H] = feval(hFunc, muBar, hParam); % Evaluating the measurement function and derivatives at muBar
K = SigmaBar*H' / (H*SigmaBar*H' + Q); % Kalman Gain
mu = muBar + K*(z - hBar);
Sigma = (eye(n)- K*H)*SigmaBar;

