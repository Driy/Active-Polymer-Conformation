module Spectral

using ....Methods.FastFourier

"""
Takes the delay sequence and the excitation sequence in Fourier space, and gives the effective resulting correlation (that is, effective correlations between instantaneous excitations that would yield the same correlation between Rouse modes).
"""
function effective_correlation(delay_sequence, excitation_sequence, J::Function; fourier_type)
    @assert length(delay_sequence) == length(excitation_sequence)
    
    # initialize matrix which will be returned
    matrix = zeros(excitation_sequence[begin] |> size)
    
    # iterate over all excitations
    for (α, C) in zip(delay_sequence, excitation_sequence)
        # iterate over all matrix elements
        for id in CartesianIndices(matrix)
            let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
                # manipulate as needed in Fourier space
                tmp = exp( -α*J(q) ) * C'[id] + C[id] * exp( -α*J(k) )
                matrix[id] += tmp                
            end
        end        
    end
    return matrix
end

"""
Takes the delay sequence and the excitation sequence in Fourier space, and gives the correction term ('B') for the fluctuations.
"""
function delayed_correlation(Δt, delay_sequence, excitation_sequence, J::Function; fourier_type)
    @assert length(delay_sequence) == length(excitation_sequence)
    
    # initialize matrix which will be returned
    matrix = zeros(excitation_sequence[begin] |> size)
    
    # iterate over all excitations
    for (α, C) in zip(delay_sequence, excitation_sequence)
        # iterate over all matrix elements
        for id in CartesianIndices(matrix)
            let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
                # manipulate as needed in Fourier space
                tmp = C[id] * exp( -α*J(k) )
                tmp /= J(q) + J(k)
                tmp *= exp( (J(q) + J(k)) * min(Δt, α) ) - 1.0
                
                # check if we divided by zero and, if yes, then take the proper limit
                if !isfinite(tmp)
                    tmp = C[id] * min(Δt, α)
                end
                matrix[id] += tmp                
            end
        end        
    end
    return matrix
end

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, τ::Real, fourier_type)

This method is used to calculate the pairwise velocity correlation between different points, with lag time τ and measurement window δ. Needs to be divided by δ^2 in the end.

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function apply_correlation_correction!(matrix::AbstractMatrix, delay_sequence, excitation_sequence, J::Function, δ::Real, τ::Real; fourier_type)
    B₀ = delayed_correlation(τ, delay_sequence, excitation_sequence, J; fourier_type = fourier_type)
    B₊ = delayed_correlation(τ+δ, delay_sequence, excitation_sequence, J; fourier_type = fourier_type)
    B₋ = delayed_correlation(abs(δ-τ), delay_sequence, excitation_sequence, J; fourier_type = fourier_type)
    
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            if τ >= δ
                matrix[id] += 2exp( -τ*J(q) ) * B₀[id] - exp( -(τ+δ)*J(q) ) * B₊[id] - exp( -(τ-δ)*J(q) ) * B₋[id]
            else
                matrix[id] += 2exp( -τ*J(q) ) * B₀[id] - exp( -(τ+δ)*J(q) ) * B₊[id] - exp( -(δ-τ)*J(k) ) * B₋'[id]
            end
        end
    end
end

end