function LogDens = DensNonLinMeasureAddErr(z, x, hFunc, hArgs, Q)

%%
% 
% PURPOSE
% --------------
% Evalute the measurement (log) density from the nonlinear measurement equation with additive Gaussian erros:
% 
% z(t) = h(u(t),x(t-1)) + eps(t), eps(t) ~ N(0,Q), and g is a function.
% 
% CALL
% --------------
% LogDens = DensNonLinMeasureAddErr(z, x, hFunc, hArgs, Q)
%        
% INPUTS
% --------------      
% z                 k-by-1          z(t), the observed data
% x                 n-by-1          x(t), the state
% hFunc             string          Name of h-function 
% hArgs             cell array      Additional arguments needed for hFunc
% Q                 k-by-k          Measurement errors covariance matrix
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
% REFERENCES
% ---------------
% Thrun, Burgard and Fox (2006). Probabilistic Robotics.

mu = feval(hFunc, x, hArgs{:});
LogDens = log(mvnpdf(z,mu,Q)); % TODO: replace with logdensity function
