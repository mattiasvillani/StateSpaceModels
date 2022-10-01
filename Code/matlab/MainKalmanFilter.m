% Main file that simulates som state space data and the runs the Kalman filter

fontSize = 14;
twoFigs = true;

% Set up the model
A = 0.9;
B = 1;
C = 1;
Q = 1; % Measurement error variance
R = 0.5; % State innovation variance
mu0 = 0;
Sigma0 = 10;
T = 100;
U = abs(1*randn(T,1));
IntervalProb = 0.95;

% Simulate state space data
[Z, X] = SimStateSpace(T, U, A, B, C, Q, R, mu0, Sigma0);

% Run the Kalman filter
[MU, SIGMA] = KalmanFilter(U, Z, A, B, C, Q, R, mu0, Sigma0);
n = size(A,1);
PostStdState = zeros(T,n);
for t = 1:T
    PostStdState(t,:) = sqrt(diag(SIGMA(:,:,t)))';
end

% Plot the data and inferences

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
legend('True state', '95% posterior intervals', 'Posterior mean', 'location','northwest')
set(gca,'fontsize',fontSize)


