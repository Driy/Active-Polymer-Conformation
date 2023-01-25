module DeltaSequence

@. function gaussian_kernel(Δs, w)
    return sqrt(π*w^2) * exp( -(Δs/w)^2 )
end

@. function sequence_kernel(τ, Δs)
    return (4π*τ)^(-0.5) * exp( - Δs^2 / (4τ) )
end

@. function effective_correlation(A, τ, w0, s0)
    return function(s1, s2)
        return A * ( gaussian_kernel(s1-s0, w0) + gaussian_kernel(s2-s0, w0) ) * sequence_kernel(τ, s2-s1)
    end
end

end