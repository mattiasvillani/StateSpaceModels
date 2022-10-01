%% Main file for the Model Velocity model for robot motion
%
% AUTHOR
% ---------------
% Mattias Villani, Linkoping University. e-mail: mattias.villani@gmail.com
%
% VERSION DATING
% ---------------
% FIRST     2016-04-14
% CURRENT   2016-04-14
%
% REFERENCES
% ---------------
% Thrun, Burgard and Fox (2006). Probabilistic Robotics.

%% User settings
deltaT = 1;
param = [1 0.01 0.1 0.1 0.01 0.01];
T = 100;            % Number of simulated time steps
U = 0.2*ones(T,2);  % The control actions
plotPath = 1;       % If = 1, robot path is plotted.
pauseTime = 0.1;    % Pause between each time step (for plotting)

%% Simulating the robot's path
if plotPath 
    hFig = figure('name','Simulated path');
    axis square
    set(gca, 'xlim',[-1 1], 'ylim',[-1 1],'xticklabels','','yticklabels','')
    hold on
    box on
end
X = zeros(T,3);
X(1,:) = 0;         % Starting at (x,y) = 0 and position along y-axis (theta=0).
for t = 2:T
    X(t,:) = SampleMotionModelVelocity(U(t,:), X(t-1,:), param, deltaT);
    if plotPath
        PlotOrientPos(X(t,:),0.03)
        line([X(t-1,1) X(t,1)],[X(t-1,2) X(t,2)],'color','g')
        drawnow
        pause(pauseTime)
    end
end

