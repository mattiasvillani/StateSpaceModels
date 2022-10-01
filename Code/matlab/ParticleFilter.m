function Psamples = ParticleFilter(U, Z, P0, TransitionSim, TransitionArgs, MeasurementLogDensity, MeasurementArgs)

%%
% 
% PURPOSE
% --------------
% Running a full particle filtering round.
%
% CALL
% --------------
% Psamples = ParticleFilter(U, Z, P0, TransSim, TransArgs, MeasureDensity, MeasureArgs)
%        
% INPUTS
% --------------
% U                 T-by-m          U(t,:) is u(t)'             
% Z                 T-by-k          Z(t,:) is z(t)'
% P0                M-by-n          Inital particles at time t = 0
% TransitionSim     string          Name of function that simulates p(x(t) \ u(t), x(t-1)). 
% TransitionArgs    cell array      Additional arguments needed for TransSim
% MeasurementLogDensity string          Name of function that evalutes p(z(t) \ x(t)). 
% MeasurementArgs       cell array      Additional arguments needed for MeasureDensity
%
% OUTPUTS
% ---------------
% Psamples          M-by-n-by-T     Psamples(m,:,t) is m:th posterior particle for x(t)                 
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
[M,n] = size(P0);

% The Kalman iterations
Psamples = zeros(M,n,T);
P = P0;
for t = 1:T
    P = ParticleFilterUpdate(P, U(t,:)', Z(t,:)', TransitionSim, TransitionArgs, MeasurementLogDensity, MeasurementArgs);
    Psamples(:,:,t) = P;
end
