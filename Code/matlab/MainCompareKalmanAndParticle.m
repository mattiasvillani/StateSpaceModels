%% Main file that simulates some state space data and the runs the Kalman filter followed by the particle filter
% Author: Mattias Villani, Linkoping University.

%% Set up the linear Gaussian state space model
%
% z(t) = C*x(t) + delta(t), delta(t) ~ N(0,Q_(t))             Measurement equation
% x(t) = A*x(t-1) + B*u(t) + eps(t), eps(t) ~ N(0,R_(t))      State equation
% where x(t) is the n-dim state, u(t) is the n-dim control, z(t) is the k-dim observed data.

A = 0.9;
B = 1;
C = 1;
Q = 10;
R = 1;
TransitionSim = 'SimNonLinStateAddErr';             % File name for state transition simulatior: p(x_t | x_t-1, u_t)
TransitionArgs = {'LinearStateEq',{A,B},R};         % Additional arguments for the TransitionSim file
MeasurementLogDensity = 'DensNonLinMeasureAddErr';  % File name for measurement density: p(z_t | x_t)
MeasurementArgs = {'LinearMeasurementEq',{C},Q};    % Additional arguments for the MeasurementLogDensity file
mu0 = 0;                                            % Mean of initial state x_0 ~ N(mu0,Sigma0^2)
Sigma0 = 10;                                        % Covariance of initial state x_0 ~ N(mu0,Sigma0^2)
T = 100;                                            % Number of time periods for the simulated data
U = abs(randn(T,1));                                % Control variables are just random numbers


%% Settings
M = 1000;                        % The number of particles
IntervalProb = 0.95;             % Interval coverage in plots
fontSize = 14;
twoFigs = false;

%% Prelims
n = size(A,1);


%% Simulate state space data
[Z, X] = SimStateSpace(T, U, A, B, C, Q, R, mu0, Sigma0);


%% Kalman filter

% Run the Kalman filter
[MU, SIGMA] = KalmanFilter(U, Z, A, B, C, Q, R, mu0, Sigma0);

% Compute the posterior standard deviations of the inferred states
n = size(A,1);
PostStdState = zeros(T,n);
for t = 1:T
    PostStdState(t,:) = sqrt(diag(SIGMA(:,:,t)))';
end


if twoFigs
    figure('name','KalmanFilterWithData')
else
    figure('name','KalmanFilter')
    subplot(2,1,1)
end

plot(X, 'linewidth', 2)
hold on
plot(MU(:,1), 'linewidth', 2)
plot(Z, 'linewidth', 2)


xlabel('time')
legend('True state','Posterior mean','Observed data', 'location','northwest')
set(gca,'fontsize',fontSize)


if twoFigs
    figure('name','KalmanFilterStateUncertainty')
else
    subplot(2,1,2)
end

lowerBand = MU(:,1) - 1.96*PostStdState(:,1);
upperBand = MU(:,1) + 1.96*PostStdState(:,1);

plot(X, 'linewidth', 2)
hold on
colorband = 0.8*[1 1 1];
patchHandle = patch([(1:T) fliplr(1:T)],[lowerBand' fliplr(upperBand')], colorband);
set(patchHandle,'faceLighting','phong','facealpha',0.5,'edgecolor',min([1.05*colorband],[1 1 1]),'edgealpha',0.0)
plot(MU(:,1), 'linewidth', 2)
xlabel('time')
legend('True state', '95% posterior intervals', 'Posterior mean (Kalman)', 'location','northwest')
set(gca,'fontsize',fontSize)






%% Particle filter

% Run the particle filter
P0 = randn(M,n);
Psamples = ParticleFilter(U, Z, P0, TransitionSim, TransitionArgs, MeasurementLogDensity, MeasurementArgs);

% Compute the posterior mean, covariance and standard deviations of the inferred states
if size(Psamples,2) == 1 % Only one state variable. Need to handle it separately
    MU = squeeze(mean(Psamples,1)); % MU is T-by-1
else
    MU = squeeze(mean(Psamples,1))'; % MU is T-by-n
end

SIGMA = zeros(n,n,T);
for t = 1:T
    SIGMA(:,:,t) = cov(Psamples(:,:,t));
end
PostStdState = zeros(T,n);
for t = 1:T
    PostStdState(t,:) = sqrt(diag(SIGMA(:,:,t)))';
end


if twoFigs
    figure('name','ParticleFilterWithData')
else
    figure('name','ParticleFilter')
    subplot(2,1,1)
end

plot(X, 'linewidth', 2)
hold on
plot(MU(:,1), 'linewidth', 2)
plot(Z, 'linewidth', 2)


xlabel('time')
legend('True state','Posterior mean','Observed data', 'location','northwest')
set(gca,'fontsize',fontSize)


if twoFigs
    figure('name','KalmanFilterStateUncertainty')
else
    subplot(2,1,2)
end

lowerBand = MU(:,1) - 1.96*PostStdState(:,1);
upperBand = MU(:,1) + 1.96*PostStdState(:,1);

plot(X, 'linewidth', 2)
hold on
colorband = 0.8*[1 1 1];
patchHandle = patch([(1:T) fliplr(1:T)],[lowerBand' fliplr(upperBand')], colorband);
set(patchHandle,'faceLighting','phong','facealpha',0.5,'edgecolor',min([1.05*colorband],[1 1 1]),'edgealpha',0.0)
plot(MU(:,1), 'linewidth', 2)
xlabel('time')
legend('True state', '95% posterior intervals', 'Posterior mean (Particle)', 'location','northwest')
set(gca,'fontsize',fontSize)

