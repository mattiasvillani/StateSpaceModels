function [Z, X] = SimStateSpace(T, U, A, B, C, Q, R, mu0, Sigma0)

%%
% 
% PURPOSE
% --------------
% Simulate data from the traditional state space model:
% 
% z(t) = C*x(t) + delta(t), delta(t) ~ N(0,Q_(t))             Measurement equation
% x(t) = A*x(t-1) + B*u(t) + eps(t), eps(t) ~ N(0,R_(t))      State equation
% where x(t) is the n-dim state, u(t) is the n-dim control, z(t) is the k-dim observed data.
% 
%
% CALL
% --------------
% Z = SimStateSpace(U, A, B, C, Q, R, mu0, Sigma0)
%        
% INPUTS
% --------------
% U, A, B, C, Q, R, mu0, Sigma0
%
% OUTPUTS
% ---------------
% Z                 T-by-k          Measurement data. Z(t,:) is z(t)'
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
% Thrun, Burgard and Fox (2006). Probabilistic Robotics.

%% Prelims
n = size(A,1);
k = size(C,1);
X = zeros(T,n);
Z = zeros(T,k);
CholR = chol(R)';
CholQ = chol(Q)';

%% Simulating the data
xLag = randmn(mu0, Sigma0); % Initial state
for t = 1:T
    x = A*xLag + B*U(t,:)' + CholR*randn(n,1);
    xLag = x;
    X(t,:) = x';
    Z(t,:) = (C*x + CholQ*randn(k,1))';
end