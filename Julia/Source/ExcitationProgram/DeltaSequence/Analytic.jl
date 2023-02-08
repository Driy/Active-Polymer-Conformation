module Analytic

"""
gaussian_kernel(Δs, w)

Corresponds to a Gaussian distribution (here of active processes) with a characteristic width.
"""
@. function gaussian_kernel(Δs, w)
    return (π*w^2)^(-0.5) * exp( -(Δs/w)^2 )
end

"""
sequence_kernel(τ, Δs)

Indicates how time-delayed excitations travel along the backbone of the polymer. Δs measures the distance and τ the elapsed time. This is basically a diffusion kernel, which widens over time.
"""
@. function sequence_kernel(τ, Δs)
    return (4π*τ)^(-0.5) * exp( - Δs^2 / (4τ) )
end

"""
effective_correlation(A, τ, w0, s0)

Calculates an effective matrix of correlated excitations that arises from a time-delayed pulse (delay time 'τ') with amplitude 'A', which has width 'w0' and is centered around 's0'.
"""
@. function effective_correlation(A, τ, w0, s0)
    return function(s1, s2)
        return A * ( gaussian_kernel(s1-s0, w0) + gaussian_kernel(s2-s0, w0) ) * sequence_kernel(τ, s2-s1)
    end
end

end