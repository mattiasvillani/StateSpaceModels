""" 
    kalmanfilter_update(μ, Σ, u, y, A, B, C, Q, R) 

A single Kalman filter update at time t of the state space model: 

yₜ = Cxₜ + δₜ,           δₜ ~ N(0,Qₜ)         Measurement equation
xₜ = Axₜ₋₁+ Buₜ + εₜ,    εₜ ~ N(0,Rₜ)         State equation

where
xₜ is the n-dim state
uₜ is the n-dim control
yₜ is the k-dim observed data. 

Reference: Thrun, Burgard and Fox (2006). Probabilistic Robotics, Algorithm Kalman_filter in Table 3.

"""
function kalmanfilter_update(μ, Σ, u, y, A, B, C, Q, R)

    # Prediction step - moving state forward without new measurement
    μ̄ = A*μ + B*u;
    Σ̄ = A*Σ*A' + R;

    # Measurement update - updating the N(μ̄, Σ̄) prior with the new data point
    K = Σ̄*C' / (C*Σ̄*C' + Q); # Kalman Gain
    μ = μ̄ + K*(y - C*μ̄);
    Σ = ( I(length(μ)) - K*C )*Σ̄;

    return μ, Σ
end

""" 
    kalmanfilter(U, Y, A, B, C, Q, R, μ₀, Σ₀) 

A full run of the Kalman filter update at time t of the state space model: 

yₜ = Cxₜ + δₜ,           δₜ ~ N(0,Qₜ)         Measurement equation
xₜ = Axₜ₋₁+ Buₜ + εₜ,    εₜ ~ N(0,Rₜ)         State equation

where
xₜ is the n-dim state
uₜ is the n-dim control
yₜ is the k-dim observed data. 

The observed data observations are the rows of the T×k matrix Y
The control signals are the rows of the T×m matrix U
μ₀ and Σ₀ are the mean and covariance of the initial state vector x₀

Reference: Thrun, Burgard and Fox (2006). Probabilistic Robotics, Algorithm Kalman_filter in Table 3.

"""

function kalmanfilter(U, Y, A, B, C, Q, R, μ₀, Σ₀)
    
    # Prelims
    T = size(Y,1)
    n = length(μ₀)
    
    μ_all = zeros(T,n)
    Σ_all = zeros(n,n,T)
    
    # The Kalman iterations
    μ = μ₀
    Σ = Σ₀
    for t = 1:T
        μ, Σ = kalmanfilter_update(μ, Σ, U[t,:]', Y[t,:]', A, B, C, Q, R)
        μ_all[t,:] = μ
        Σ_all[:,:,t] = Σ
    end

    return μ_all, Σ_all
end
    