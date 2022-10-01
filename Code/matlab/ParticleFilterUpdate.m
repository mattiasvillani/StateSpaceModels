function P = ParticleFilterUpdate(P, u, z, TransitionSim, TransitionArgs, MeasurementLogDensity, MeasurementArgs)

%%
% 
% PURPOSE
% --------------
% A single particle filter update at time t for the model.
% 
% CALL
% --------------
% [mu, Sigma] = KalmanFilterUpdate(mu, Sigma, u, z, A, B, C, Q, R)
%        
% INPUTS
% --------------
% P                         M-by-n          Set of M particles at time t-1                
% u                         n-by-1          u(t), the controls at time t
% z                         k-by-1          z(t), the measurments at time t
% TransitionSim             string          Name of function that simulates p(x(t) \ u(t), x(t-1)). 
% TransitionArgs            cell array      Additional arguments needed for TransSim
% MeasurementLogDensity     string          Name of function that evalutes p(z(t) \ x(t)). 
% MeasurementArgs           cell array      Additional arguments needed for MeasureDensity

% OUTPUTS
% ---------------
% P                         M-by-n          Set of M particles at time t           
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
% Thrun, Burgard and Fox (2006). Probabilistic Robotics, Algorithm Particle_filter in Table 4.3.

%% Prelims
M = size(P,1);
Wprop = zeros(M,1);

for m = 1:M
    P(m,:) = feval(TransitionSim, u, P(m,:)', TransitionArgs{:});  % Proposal draws for time t. Propagating particles forward.
    Wprop(m) = feval(MeasurementLogDensity, z, P(m,:)', MeasurementArgs{:});
end
Wprop = exp(Wprop);
Wprop = Wprop/sum(Wprop);

SelectedObs = randsample(1:M, M, true, Wprop);
P = P(SelectedObs,:);


