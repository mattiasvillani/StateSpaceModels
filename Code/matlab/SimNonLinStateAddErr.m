function x = SimNonLinStateAddErr(u, xLag, gFunc, gArgs, R)

%%
% 
% PURPOSE
% --------------
% A simulate from the nonlinear state equation with additive Gaussian erros:
% 
% x(t) = g(u(t),x(t-1)) + eps(t), eps(t) ~ N(0,R), and g is a function.
% 
% CALL
% --------------
% x = SimNonLinStateAddErr(u, xLag, gFunc, gArgs, R)
%        
% INPUTS
% --------------               
% u                 m-by-1          u(t), the controls at time t
% xLag              n-by-1          x(t-1)
% gFunc             string          Name of g-function 
% gArgs             cell array      Additional arguments needed for gFunc
% R                 n-by-n          System innovation covariance matrix
%
% OUTPUTS
% ---------------
% x                 n-by-1          Simulation of state x(t)         
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

n = size(R,1);
x = feval(gFunc, u, xLag, gArgs{:}) + chol(R)'*randn(n,1);