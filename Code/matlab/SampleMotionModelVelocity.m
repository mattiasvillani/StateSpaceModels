function xt = SampleMotionModelVelocity(ut,xLag,param,deltaT)


%%
% 
% PURPOSE
% --------------
% Simulate a time step with the velocity motion model (2D case)
% 
% CALL
% --------------
% xt = SampleMotionModelVelocity(u,xLag,param)
%        
% INPUTS
% --------------
% ut                        1-by-2          ut = (vt,wt), translational and rotational velocities at time t               
% xLag                      1-by-3          xLag = (x,y,theta), current pose (coordinates x,y and direction theta)
% param                     1-by-6          Parameters (alpha_1, ,..., alpha_6) of the motion velocity model
% deltaT                    scalar          The time step (default = 1).
% OUTPUTS
% ---------------
% xt                        1-by-3          Simulated pose (xTilde,yTilde,ThetaTilde) at time t          
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
% Thrun, Burgard and Fox (2006). Probabilistic Robotics
% Algorithm sample_motion_model_velocity in Table 5.3.

if nargin == 3
    deltaT = 1; 
end

v = ut(1); % translational velocity
w = ut(2); % rotational velocity
x = xLag(1);
y = xLag(2);
theta = xLag(3);

vHat = randn*sqrt((param(1)*v^2 + param(2)*w^2)); % Random intial translation
wHat = randn*sqrt((param(3)*v^2 + param(4)*w^2)); % Random rotation
gammaHat = randn*sqrt((param(5)*v^2 + param(6)*w^2)); % Additional rotation to break degeneracy

xTilde = x - (vHat/wHat)*sin(theta) + (vHat/wHat)*sin(theta + wHat*deltaT);
yTilde = y + (vHat/wHat)*cos(theta) - (vHat/wHat)*cos(theta + wHat*deltaT);
thetaTilde = theta + wHat*deltaT + gammaHat*deltaT;

xt = [xTilde yTilde thetaTilde];
